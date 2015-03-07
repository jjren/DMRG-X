MODULE InitialGuess
	use mpi
	use variables
	implicit none
! this module contains the subroutine that generate the initial guess of the MPS
! procedure

	contains
!===========================================================
	Subroutine Initialfinit(guesscoeff,direction)
	! when nstate=1 then we can use the last stored matrix to
	! contruct the Initial Guess
	! only used in the finit MPS and nstate=1
	use mpi
	use variables
	USE BLAS95
	USE F95_PRECISION

! ngoodstates is the number of states fullfill good quantum number
	implicit none
	real(kind=8) :: guesscoeff(ngoodstates)
	real(kind=8),allocatable :: leftu(:,:),rightv(:,:),singularvalue(:)&
	,LRcoeff(:,:)
	logical :: alive
	integer :: reclength
	integer :: error,i,j,m
	real(kind=8) :: norm
	character(len=1) :: direction
	logical :: done

	if(myid==0) then
		write(*,*) "enter Initialfinit subroutine"
		! two site dmrg
		if((nright+nleft+2)/=norbs) then
			write(*,*) "-----------------------------------"
			write(*,*) "two site dmrg nright+nleft+2/=norbs"
			write(*,*) "-----------------------------------"
			stop
		end if

		allocate(leftu(4*Lrealdim,subM),stat=error)
		if(error/=0) stop
		allocate(rightv(subM,4*Rrealdim),stat=error)
		if(error/=0) stop
		allocate(singularvalue(subM),stat=error)
		if(error/=0) stop

		reclength=2*subM*subM

		inquire(file="wavefunction.tmp",exist=alive)
		if(alive) then
			open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) "wavefunction.tmp doesn't exist"
			stop
		end if

		do i=1,4,1
			read(105,rec=4*nleft+i) leftu((i-1)*Lrealdim+1:i*Lrealdim,1:subM)
		end do
		do i=1,4,1
			read(105,rec=4*(nleft+1)+i) rightv(1:subM,i:4*Rrealdim:4)
		end do

		open(unit=106,file="singularvalue.tmp",status="old")
		read(106,*) singularvalue(1:subM)
! let the singularvalue'sum be 1
		norm=sum(singularvalue(1:subM))
		singularvalue=singularvalue/norm
		singularvalue=sqrt(singularvalue)

! be careful the finit initial guessvector
		if(direction=='l' .and. nleft>exactsite) then
			if(Lrealdim/=subM) then
				write(*,*) "-----------------------------------"
				write(*,*) "guessfinit, Lrealdim/=subM failed!",Lrealdim,subM
				write(*,*) "-----------------------------------"
				stop
			end if

			do i=1,subM,1
				leftu(i:4*subM:subM,:)=leftu(i:4*subM:subM,:)*singularvalue(i)
			end do

		else if(direction=='r' .and. nright>exactsite) then
			if(Rrealdim/=subM) then
				write(*,*) "-----------------------------------"
				write(*,*) "guessfinit, Rreadlim/=subM failed!"
				write(*,*) "-----------------------------------"
				stop
			end if

			do i=1,subM,1
				rightv(:,(i-1)*4+1:i*4)=rightv(:,(i-1)*4+1:i*4)*singularvalue(i)
			end do

		else if((direction=='l' .and. nleft==exactsite) .or. (direction=='r' .and. nright==exactsite)) then

			do i=1,subM,1
				rightv(i,:)=rightv(i,:)*singularvalue(i)
			end do
		end if
			
		! recombine the two site sigmaL sigmaR coefficient
		
		allocate(LRcoeff(4*Lrealdim,4*Rrealdim),stat=error)
		if(error/=0) stop
		call gemm(leftu,rightv,LRcoeff,'N','N',1.0D0,0.0D0)

		!if(logic_C2/=0) then
		!	if((direction=='l' .and. nleft==exactsite) .or. (direction=='r' .and. nright==exactsite)) then
		!		write(*,*)
		!	else if(direction=='l') then
		!		do i=subM+1,3*subM,1
		!			LRcoeff(i,:)=LRcoeff(i,:)*DBLE((-1)**quantabigR(:,1))
		!		end do
		!	else if(direction=='r') then
		!		do i=subM+1,3*subM,1
		!			LRcoeff(:,i)=LRcoeff(:,i)*DBLE((-1)**(quantabigR(:,1)-1))
		!		end do
		!	end if
		!end if

		m=1
		do i=1,4*Rrealdim,1
		do j=1,4*Lrealdim,1
			if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
				quantabigL(j,2)+quantabigR(i,2)==totalSz) then
				guesscoeff(m)=LRcoeff(j,i)
				m=m+1
			end if
		end do
		end do

		m=m-1

		if(m/=ngoodstates) then
			write(*,*) "----------------------------------------------"
			write(*,*) "guesscoeff good quantum states number wrong! failed!"
			write(*,*) "----------------------------------------------"
			stop
		end if

		if(logic_spinreversal/=0) then
			do i=1,ngoodstates,1
			! the symmlink is himself
				if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1) .and. &
				abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(i,2)) then
					if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
					/=logic_spinreversal*((-1)**(mod(nelecs/2,2)))) then
						guesscoeff(i)=0.0D0
					end if
			! either of the L or R space symmlink is not himself
				else if(abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1) .or. &
				abs(symmlinkbig(symmlinkgood(i,2),1,2))/=symmlinkgood(i,2)) then
					done=.false.
					do j=1,ngoodstates,1
						if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(j,1) .and. &
						abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(j,2)) then
							guesscoeff(j)=guesscoeff(i)&
							*DBLE(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2)))&
							*DBLE(logic_spinreversal)*((-1)**(mod(nelecs/2,2)))
							done=.true.
							exit
						end if
					end do
					if(done==.false.) then
						write(*,*) "-----------------------------------------------------------------------------"
						write(*,*) "did't find the state",i,symmlinkgood(i,1),symmlinkgood(i,2),"corrsponds",symmlinkbig(symmlinkgood(i,1),1,1),symmlinkbig(symmlinkgood(i,2),1,2)
						write(*,*) "-----------------------------------------------------------------------------"
						stop
					end if
				end if
			end do
		end if

! In the MPS form, no simple C2 symmetry
! C2 symmetry part
! the rule is 
! if L space fai1 fai2, corresponding R space fai3,fai4
! C2(fai1*fai4)=fai3*fai2=-1^(num3*num2)*fai2*fai3
	!	if(logic_C2/=0 .and. nleft==nright) then
	!		write(*,*) nleft,nright,"nleft"
	!		do i=1,ngoodstates,1
	!			if(symmlinkgood(i,1)==symmlinkgood(i,2)) then
	!				if(logic_C2*(-1)**mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)==-1) then
	!					guesscoeff(i)=0.0D0
	!				end if
	!			else
	!				done=.false.
	!				do j=1,ngoodstates,1
	!				if(symmlinkgood(j,1)==symmlinkgood(i,2) .and. &
	!				symmlinkgood(j,2)==symmlinkgood(i,1)) then
	!					guesscoeff(j)=guesscoeff(i)*&
	!					DBLE(logic_C2*(-1)**(mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)))
	!					done=.true.
	!					exit
	!				end if
	!				end do
	!				if(done==.false.) then
	!					write(*,*) "---------------------------------------"
	!					write(*,*) "in Initial finit guess C2 done=.false."
	!					write(*,*) "---------------------------------------"
	!					stop
	!				end if
	!			end if
	!		end do
	!	end if

		norm=dot(guesscoeff(1:ngoodstates),guesscoeff(1:ngoodstates))
		write(*,*) "finit guesscoeff norm",norm
		if(norm<1.0D-10) then
			write(*,*) "--------------------------"
			write(*,*) "initialfinit norm is < 1.0D-10,caution!"
			write(*,*) "--------------------------"
		end if
		guesscoeff(1:ngoodstates)=guesscoeff(1:ngoodstates)/sqrt(norm)


		close(105)
		close(106)


		deallocate(leftu)
		deallocate(rightv)
		deallocate(singularvalue)
		deallocate(LRcoeff)
	end if
	return
	end subroutine

!================================================================

! Initial Guess in the infinit DMRG and nstate/=1 DMRG
! every state fullfilled good quantum number have random weight

	Subroutine Initialrandomweight(guessvector,num)
	use mpi
	use variables
	use BLAS95
	use f95_precision

	implicit none

	integer :: num
	real(kind=8) :: guessvector(num*ngoodstates),randomx,norm
	integer :: i,j,k,error,m,total,exsit
	logical :: done,Ifexsit
	integer :: spinline(ngoodstates,2),C2line(ngoodstates,2),fulline(4,2)

	
	if(myid==0) then
		write(*,*) "enter Initialrandomweight subroutine"
		
		spinline=10
		C2line=10
		fulline=10

		norm=0.0D0
		call random_seed()

		do i=1,num*ngoodstates,1
			call random_number(randomx)
			guessvector(i)=randomx
		end do

! when L space Sz>0 then we need to make the symmetry pair using same coefficient
! when L space Sz=0 then we need to make the L+R space basis to have the specfic 
! spin parity. Then set others to zero
		if(logic_spinreversal/=0) then
		!	do i=1,ngoodstates,1
		!		if(quantabigL(symmlinkgood(i,1),2)>=0 .and. abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1)) then
		!			done=.false.
		!		do j=1,ngoodstates,1
		!			if(symmlinkgood(j,1)==abs(symmlinkbig(symmlinkgood(i,1),1,1)) &
		!				.and. symmlinkgood(j,2)==abs(symmlinkbig(symmlinkgood(i,2),1,2))) then
		!				guessvector(j:num*ngoodstates:ngoodstates)=guessvector(i:num*ngoodstates:ngoodstates)&
		!				*DBLE(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2)))&
		!				*DBLE(logic_spinreversal)
		!				done=.true.
		!				exit
		!			end if
		!		end do
		!			if(done/=.true.) then
		!				write(*,*) "-------------------------------------------------"
		!				write(*,*) "initialrandomweight spin reversal adapted failed!"
		!				write(*,*) "-------------------------------------------------"
		!				stop
		!			end if
		!		else if(quantabigL(symmlinkgood(i,1),2)==0 .and. &
		!		abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1)) then
		!			if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))/=logic_spinreversal) then
		!				guessvector(i:num*ngoodstates:ngoodstates)=0.0D0
		!			end if
		!		end if
		!	end do
			do i=1,ngoodstates,1
			! the symmlink is himself
				if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1) .and. &
				abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(i,2)) then
					if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
					/=logic_spinreversal*((-1)**(mod(nelecs/2,2)))) then
						guessvector(i:num*ngoodstates:ngoodstates)=0.0D0
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
						write(*,*) "in Initial did't find the state",i,symmlinkgood(i,1),symmlinkgood(i,2),"corrsponds",symmlinkbig(symmlinkgood(i,1),1,1),symmlinkbig(symmlinkgood(i,2),1,2)
						write(*,*) "-----------------------------------------------------------------------------"
						stop
					end if
				end if
			end do
		end if

! C2 symmetry part
! the rule is 
! if L space fai1 fai2, corresponding R space fai3,fai4
! C2(fai1*fai4)=fai3*fai2=-1^(num3*num2)*fai2*fai3
		if(logic_C2/=0 .and. nleft==nright) then
			do i=1,ngoodstates,1
				if(symmlinkgood(i,1)==symmlinkgood(i,2)) then
					if(logic_C2*(-1)**mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)==-1) then
						guessvector(i:num*ngoodstates:ngoodstates)=0.0D0
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
						guessvector(fulline(k,1):num*ngoodstates:ngoodstates)=0.0D0
					end do
				else
					do k=2,total,1
						guessvector(fulline(k,1):num*ngoodstates:ngoodstates)=&
							guessvector(fulline(1,1):num*ngoodstates:ngoodstates)*DBLE(fulline(k,2)/fulline(1,2))
					end do
				end if

			end do
		end if
	!	write(*,*) "guessvector"
	!	write(*,*) guessvector

		norm=dot(guessvector(1:ngoodstates),guessvector(1:ngoodstates))
		if(norm<1.0D-10) then
			write(*,*) "--------------------------"
			write(*,*) "norm is < 1.0D-10,caution!"
			write(*,*) "--------------------------"
		end if
		guessvector(1:ngoodstates)=guessvector(1:ngoodstates)/sqrt(norm)
		
! Gram-Schmit Orthogonalization
		if(num >= 2) then
		do i=2,num,1
			do j=1,i-1,1
				norm=dot(guessvector((i-1)*ngoodstates+1:i*ngoodstates),&
				guessvector((j-1)*ngoodstates+1:j*16*ngoodstates))
				guessvector((i-1)*ngoodstates+1:i*ngoodstates)=&
					guessvector((i-1)*ngoodstates+1:i*ngoodstates)-&
					norm*guessvector((j-1)*ngoodstates+1:j*ngoodstates)
			end do
				norm=dot(guessvector((i-1)*ngoodstates+1:i*ngoodstates),&
				guessvector((i-1)*ngoodstates+1:i*ngoodstates))
				if(norm<1.0D-10) then
					write(*,*) "norm is < 1.0D-10,caution!"
				end if
				guessvector((i-1)*ngoodstates+1:i*ngoodstates)=&
					guessvector((i-1)*ngoodstates+1:i*ngoodstates)/sqrt(norm)
		end do
		end if

	end if

	return
	end subroutine
!===============================================================



!================================================================

! unit vector that corresponding the smallest diagonal element!
! because we want to get the smallest eigenvalue
! but this method may not good in PPP model
! because initially the R space and L space does not have non-diagonal
! interation

	Subroutine Initialunitvector(HDIAG,guessvector,num)
	use mpi
	use variables
	implicit none

	integer :: num
	real(kind=8) :: guessvector(num*ngoodstates)
	real(kind=8) :: HDIAG(ngoodstates)
	integer :: i,j,k,l,minindex(num)
	real(kind=8) :: mindiag(num)

	
	if(myid==0) then
		write(*,*) "enter Initialunitvector subroutine"
! the first is the smallest number
		minindex=0
		mindiag=1.0D10
		do i=1,ngoodstates,1
			do k=1,num,1
				if(HDIAG(i)<mindiag(k)) then
					do l=num,k,-1
						if(l==k) then
							mindiag(l)=HDIAG(i)
							minindex(l)=i
						else
							mindiag(l)=mindiag(l-1)
							minindex(l)=minindex(l-1)
						end if
					end do
					exit
				end if
			end do
		end do

		guessvector=0.0D0
		do i=1,num,1
			if(minindex(i)==0) then
				write(*,*) "---------------------------------------------"
				write(*,*) " initialunivector good quantum number failed!"
				write(*,*) "---------------------------------------------"
				stop
			end if
		end do
		do i=1,num,1
			guessvector((i-1)*ngoodstates+minindex(i))=1.0D0
		end do
	end if

	return
	end subroutine
!===============================================================

end MODULE
