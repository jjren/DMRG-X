	Subroutine statecorrect(Davidwork,checksymm)
	use mpi
	use variables
	use BLAS95
	use f95_precision

	implicit none

	real(kind=8) :: Davidwork(nstate*ngoodstates),norm,checksymm
	integer :: i,j,k,error,m,total,exsit
	logical :: done,Ifexsit
	integer :: spinline(ngoodstates,2),C2line(ngoodstates,2),fulline(4,2)

	
	if(myid==0) then
		write(*,*) "enter statecorrect subroutine"
		
		spinline=10
		C2line=10
		fulline=10

		norm=0.0D0

		if(logic_spinreversal/=0) then
			do i=1,ngoodstates,1
			! the symmlink is himself
				if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1) .and. &
				abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(i,2)) then
					if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
					/=logic_spinreversal*((-1)**(mod(nelecs/2,2)))) then
						Davidwork(i:nstate*ngoodstates:ngoodstates)=0.0D0
						spinline(i,1)=i
						spinline(i,2)=0
					else
						spinline(i,1)=i
						spinline(i,2)=1
					end if
			! either of the L or R space symmlink is not himself
				else if(abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1) .or. &
				abs(symmlinkbig(symmlinkgood(i,2),1,2))/=symmlinkgood(i,2)) then
					done=.false.
					do j=1,ngoodstates,1
						if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(j,1) .and. &
						abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(j,2)) then
							done=.true.
							if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
							*logic_spinreversal*((-1)**(mod(nelecs/2,2)))==-1) then
								spinline(i,1)=j
								spinline(j,1)=i
								spinline(i,2)=-1
								spinline(j,2)=-1
							else if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
							*logic_spinreversal*((-1)**(mod(nelecs/2,2)))==1) then
								spinline(i,1)=j
								spinline(j,1)=i
								spinline(i,2)=1
								spinline(j,2)=1
							else
								write(*,*) "failed!"
								write(*,*) sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
								*logic_spinreversal*((-1)**(mod(nelecs/2,2)))
								stop
							end if
							exit
						end if
					end do
					if(done==.false.) then
						write(*,*) "-----------------------------------------------------------------------------"
						write(*,*) "in statecorrect did't find the state",i,symmlinkgood(i,1),symmlinkgood(i,2),"corrsponds",symmlinkbig(symmlinkgood(i,1),1,1),symmlinkbig(symmlinkgood(i,2),1,2)
						write(*,*) "-----------------------------------------------------------------------------"
						stop
					end if
				end if
			end do
		end if

! C2 symmetry part
! the rule is 
! if L space fai1 fai2, corresponding R space fai3,fai4
! C2(fai1*fai4)=fai3*fai2=-1^(nstate3*nstate2)*fai2*fai3
		if(logic_C2/=0 .and. nleft==nright) then
			do i=1,ngoodstates,1
				if(symmlinkgood(i,1)==symmlinkgood(i,2)) then
					if(logic_C2*(-1)**mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)==-1) then
						Davidwork(i:nstate*ngoodstates:ngoodstates)=0.0D0
						C2line(i,1)=i
						C2line(i,2)=0
					else
						C2line(i,1)=i
						C2line(i,2)=1
					end if
				else
					done=.false.
					do j=1,ngoodstates,1
					if(symmlinkgood(j,1)==symmlinkgood(i,2) .and. &
					symmlinkgood(j,2)==symmlinkgood(i,1)) then
						if(logic_C2*(-1)**mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)==-1) then
							C2line(i,1)=j
							C2line(j,1)=i
							C2line(i,2)=-1
							C2line(j,2)=-1
						else
							C2line(i,1)=j
							C2line(j,1)=i
							C2line(i,2)=1
							C2line(j,2)=1
						end if
						done=.true.
						exit
					end if
					end do
					if(done==.false.) then
						write(*,*) "---------------------------------------"
						write(*,*) "in Initial random guess C2 done=.false."
						write(*,*) "---------------------------------------"
						stop
					end if
				end if
			end do
		end if

! the core part! find the link diagram
		if((logic_C2/=0 .and. nleft==nright) .or. logic_spinreversal/=0) then
			do i=1,ngoodstates,1
				m=1
				total=1
				fulline(1,1)=i
				! index
				fulline(1,2)=1
				! phase

				do while(.true.)
				if(logic_C2/=0 .and. nleft==nright) then
					Ifexsit=.false.
					do j=1,total,1
						if(C2line(fulline(m,1),1)==fulline(j,1)) then
							Ifexsit=.true.
							exsit=j
							exit
						end if
					end do
					if(Ifexsit==.false.) then
						fulline(total+1,1)=C2line(fulline(m,1),1)
						fulline(total+1,2)=fulline(m,2)*C2line(fulline(m,1),2)
						total=total+1
					else if(C2line(fulline(m,1),2)*fulline(m,2)/=fulline(exsit,2)) then
						do j=1,total,1
							fulline(j,2)=0
						end do
					end if
				end if

				if(logic_spinreversal/=0) then
					Ifexsit=.false.
					do j=1,total,1
						if(spinline(fulline(m,1),1)==fulline(j,1)) then
							Ifexsit=.true.
							exsit=j
							exit
						end if
					end do
					if(Ifexsit==.false.) then
						fulline(total+1,1)=spinline(fulline(m,1),1)
						fulline(total+1,2)=fulline(m,2)*spinline(fulline(m,1),2)
						total=total+1
					else if(spinline(fulline(m,1),2)*fulline(m,2)/=fulline(exsit,2)) then
						do j=1,total,1
							fulline(j,2)=0
						end do
					end if
				end if
				m=m+1
				if(m>total) exit
				end do

				if(fulline(1,2)==0) then
					do k=1,total,1
						Davidwork(fulline(k,1):nstate*ngoodstates:ngoodstates)=0.0D0
					end do
				else
					do k=2,total,1
						Davidwork(fulline(k,1):nstate*ngoodstates:ngoodstates)=&
							Davidwork(fulline(1,1):nstate*ngoodstates:ngoodstates)*DBLE(fulline(k,2)/fulline(1,2))
					end do
				end if

			end do
		end if
!		write(*,*) "Davidwork"
!		write(*,*) Davidwork

		norm=dot(Davidwork(1:ngoodstates),Davidwork(1:ngoodstates))
		write(*,*) "statecorrect state1 norm=",norm
		checksymm=norm
		if(norm<1.0D-10) then
			write(*,*) "--------------------------"
			write(*,*) "norm is < 1.0D-10,caution!"
			write(*,*) "--------------------------"
		end if
		Davidwork(1:ngoodstates)=Davidwork(1:ngoodstates)/sqrt(norm)
		
! Gram-Schmit Orthogonalization
		if(nstate >= 2) then
		do i=2,nstate,1
			do j=1,i-1,1
! checkout if the symmetry is correct
				norm=dot(Davidwork((i-1)*ngoodstates+1:i*ngoodstates),&
				Davidwork((i-1)*ngoodstates+1:i*ngoodstates))
				write(*,*) "statecorrect state",i,"norm=",norm

				norm=dot(Davidwork((i-1)*ngoodstates+1:i*ngoodstates),&
				Davidwork((j-1)*ngoodstates+1:j*16*ngoodstates))
				Davidwork((i-1)*ngoodstates+1:i*ngoodstates)=&
					Davidwork((i-1)*ngoodstates+1:i*ngoodstates)-&
					norm*Davidwork((j-1)*ngoodstates+1:j*ngoodstates)
			end do
				norm=dot(Davidwork((i-1)*ngoodstates+1:i*ngoodstates),&
				Davidwork((i-1)*ngoodstates+1:i*ngoodstates))
				if(norm<1.0D-10) then
					write(*,*) "norm is < 1.0D-10,caution!"
				end if
				Davidwork((i-1)*ngoodstates+1:i*ngoodstates)=&
					Davidwork((i-1)*ngoodstates+1:i*ngoodstates)/sqrt(norm)
		end do
		end if

	end if

	return
	end subroutine
!===============================================================
