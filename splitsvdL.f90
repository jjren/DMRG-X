subroutine splitsvdL(singularvalue,leftu,statebegin,stateend,indexlp1)
! this subroutine is used to split the reduced density matrix
! to different subspace according to good quantum number
! and diagonalizaiton it to get the renormalized vector
	use mpi
	use variables
	USE blas95
	use lapack95
	USE F95_PRECISION

	implicit none

	real(kind=8) :: leftu(4*Lrealdim,subM),singularvalue(subM),discard
	integer :: statebegin,stateend,indexlp1
	real(kind=8),allocatable :: valuework(:),coeffwork(:,:),coeffbuffer(:,:)&
	,coeffresult(:,:)
	integer,allocatable :: valueindex(:),szzeroindex(:),quantabigLbuffer(:,:)
	integer :: i,error,j,k,l,m,n,info,i1,symm1,p,p1
	integer :: szzero,pair1,szzero1,szl0
	logical :: done,done2
	integer :: tmp
! szl0 means the number of sz>0 states
! szzero means the number of sz=0 state including spin parity=+-1
! pair1 means the number of sz>0 and sz=0 but the symmetry pair is not himself state
! szzero1 means the number sz=0 and spin reversal parity=1


	write(*,*) "enter in subroutine splitsvdL"

	allocate(coeffbuffer(4*subM,4*subM),stat=error)
	if(error/=0) stop
	allocate(coeffwork(4*subM,4*subM),stat=error)
	if(error/=0) stop
	allocate(coeffresult(4*subM,4*subM),stat=error)
	if(error/=0) stop
	
	coeffbuffer=0.0D0
	coeffresult=0.0D0

	write(*,*) "coeffIF"
	write(*,*) coeffIF(1:4*Lrealdim,1:4*Rrealdim,:)

! the L+sigmaL space reduced density matrix
	do i=statebegin,stateend,1
		call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),&
		coeffbuffer(1:4*Lrealdim,1:4*Lrealdim),'N','T',1.0D0,0.0D0)
		if(exscheme==1) then
			coeffwork=coeffwork+coeffbuffer*nweight(i)
		else
			coeffwork=coeffwork+coeffbuffer
		end if
	end do
	write(*,*) "coeffwork"
	write(*,*) coeffwork
! split the reduced density matrix to different good quantum
! number subspace
	
	allocate(valuework(4*subM),stat=error)
	if(error/=0) stop
	valuework=0.0D0
	allocate(valueindex(subM),stat=error)
	if(error/=0) stop
	allocate(quantabigLbuffer(4*subM,2),stat=error)
	if(error/=0) stop
	
	if(logic_spinreversal/=0) then
	allocate(szzeroindex(4*subM),stat=error)
	if(error/=0) stop
	szzeroindex=0
	end if
	
	coeffbuffer=coeffwork
	m=0
! m is the number of good quantum state so far
! j is Sz
! i is nelecs
	do j=nleft+1,-nleft-1,-1
		
		if(logic_spinreversal/=0 .and. j<0) then
			exit
		! using pair symmetry copy
		end if

! when Sz=0 we need to seperate spin reversal=1 or -1 and symmetry pair is not
! himself
		if(logic_spinreversal/=0 .and. j==0) then
			symm1=3
		else
			symm1=1
		end if
				

		do i1=1,symm1,1
		
		do i=0,2*(nleft+1),1

		
			n=m
			coeffwork=coeffbuffer
! n is the last total m
			if(symm1==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			!	write(*,*) m
			else if(symm1==3 .and. i1==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. abs(symmlinkbig(k,1,1))/=k) then
						done2=.true.
						do p1=n+1,m,1
							if(abs(symmlinkbig(k,1,1))/=szzeroindex(p1)) then
								done2=.true.
							else
								done2=.false.
								exit
							end if
						end do
						if(done2==.true.) then
							m=m+1
							szzeroindex(m)=k
							call swap(coeffwork(:,m),coeffwork(:,k))
							call swap(coeffwork(m,:),coeffwork(k,:))
						end if
					end if
				end do
			!	write(*,*) m
			else if(symm1==3 .and. i1==2) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==k) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			!	write(*,*) m
			else if(symm1==3 .and. i1==3) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==-k) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			!	write(*,*) m
			end if

			if(m/=n) then
				call syevd(coeffwork(n+1:m,n+1:m),valuework(n+1:m),'V','U',info)
				if(info/=0) then
					write(*,*) "left diagnolization failed!"
					stop
				end if
				coeffresult(n+1:m,n+1:m)=coeffwork(n+1:m,n+1:m)
				
			p=n+1
			if(symm1==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
						call copy(coeffresult(p,n+1:m),coeffresult(k,n+1:m))
						coeffresult(p,n+1:m)=0.0D0
						p=p+1
					end if
				end do
			else if(symm1==3 .and. i1==1) then
					do p1=1,4*Lrealdim
						if(szzeroindex(p1)/=0) then
							call copy(coeffresult(p,n+1:m),coeffresult(szzeroindex(p1),n+1:m))
							coeffresult(p,n+1:m)=0.0D0
							p=p+1
						end if
					end do
			else if(symm1==3 .and. i1==2) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==k) then
						call copy(coeffresult(p,n+1:m),coeffresult(k,n+1:m))
						coeffresult(p,n+1:m)=0.0D0
						p=p+1
					end if
				end do
			else if(symm1==3 .and. i1==3) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==-k) then
						call copy(coeffresult(p,n+1:m),coeffresult(k,n+1:m))
						coeffresult(p,n+1:m)=0.0D0
						p=p+1
					end if
				end do
			end if

				if(p-1/=m) then
					write(*,*) "p-1/=m-n failed!",p-1,m-n
					stop
				end if
				quantabigLbuffer(n+1:m,1)=i
				quantabigLbuffer(n+1:m,2)=j
			end if
		end do
		
		if(j>0) then
			szl0=m
		else if(j==0 .and. i1==1 .and. logic_spinreversal/=0) then
			pair1=m
		else if(j==0 .and. i1==2 .and. logic_spinreversal/=0) then
			szzero1=m-pair1
		else if(j==0 .and. i1==3 .and. logic_spinreversal/=0) then
			szzero=m-pair1
		end if

		end do

	end do
	
	if(logic_spinreversal/=0) then
	if((pair1*2+szzero)/=4*Lrealdim) then
		write(*,*) "(pair1*2+szzero)/=4*Lrealdim failed!","pair1=",pair1,"szzero",szzero
		stop
	end if
	end if
	


! copy the Sz>0 part and Sz=0(symmetry pair is not himself) to the symmetry pair
	if(logic_spinreversal/=0) then
		do i=1,4*Lrealdim,1
			if(quantabigL(i,2)>0) then
				coeffresult(abs(symmlinkbig(i,1,1)),szl0+szzero+1:2*szl0+szzero)=coeffresult(i,1:szl0)*DBLE(sign(1,symmlinkbig(i,1,1)))
			end if
		end do
			do p1=1,4*Lrealdim,1
				if(szzeroindex(p1)/=0) then
				coeffresult(abs(symmlinkbig(szzeroindex(p1),1,1)),pair1+szzero+szl0+1:4*Lrealdim)&
				=coeffresult(szzeroindex(p1),szl0+1:pair1)*&
					DBLE(sign(1,symmlinkbig(szzeroindex(p1),1,1)))
				end if
			end do
		quantabigLbuffer(pair1+szzero+1:4*Lrealdim,:)=quantabigLbuffer(1:pair1,:)
	end if

	singularvalue=0.0D0
	valueindex=0
	if(logic_spinreversal==0) then
		do i=1,4*Lrealdim,1
			do j=1,subM,1
				if(valuework(i)>singularvalue(j)) then
					valueindex(j+1:subM)=valueindex(j:subM-1)
					singularvalue(j+1:subM)=singularvalue(j:subM-1)
					valueindex(j)=i
					singularvalue(j)=valuework(i)
					exit
				end if
			end do
		end do
	else
		do i=1,szzero+pair1,1
			do j=1,subM,1
				if(valuework(i)>singularvalue(j)) then
					if(i<=pair1) then
						valueindex(j+2:subM)=valueindex(j:subM-2)
						valueindex(j)=i
						valueindex(j+1)=pair1+szzero+i
						singularvalue(j+2:subM)=singularvalue(j:subM-2)
						singularvalue(j)=valuework(i)
						singularvalue(j+1)=valuework(i)
					else
						valueindex(j+1:subM)=valueindex(j:subM-1)
						singularvalue(j+1:subM)=singularvalue(j:subM-1)
						valueindex(j)=i
						singularvalue(j)=valuework(i)
					end if
				end if
			end do
		end do

		if(abs(symmlinkbig(valueindex(subM),1,1))/=symmlinkbig(valueindex(subM),1,1) .or. &
		abs(symmlinkbig(valueindex(subM),1,1))/=symmlinkbig(valueindex(subM-1),1,1)) then
			do i=pair1+szzero,pair1+1,-1
			
			do j=1,subM,1
				if(valueindex(j)==i) then
					done=.false.
					exit
				else
					done=.true.
				end if
			end do
				if(done==.true.) then
				valueindex(subM)=i
				singularvalue(subM)=valuework(i)
				exit
				end if
			end do
		end if
	end if

	do i=1,subM,1
		leftu(:,i)=coeffresult(1:4*Lrealdim,valueindex(i))
		quantasmaL(i,:)=quantabigLbuffer(valueindex(i),:)
		if(logic_spinreversal/=0) then
			if(valueindex(i)<=pair1) then
				symmlinksma(i,1,1)=i+1
			else if(valueindex(i)>pair1 .and. valueindex(i)<=pair1+szzero) then
				if(valueindex(i)<=pair1+szzero1) then
					symmlinksma(i,1,1)=i
				else
					symmlinksma(i,1,1)=-i
				end if
			else 
				symmlinksma(i,1,1)=i-1
			end if
		end if
	end do
		write(*,*) valuework
		write(*,*) singularvalue
		write(*,*) valueindex
!   the total discard weight between site indexLp1,indexRm1
		discard=1.0D0-sum(singularvalue(1:subM))
		write(*,'(A20,I4,D12.5)') "totaldiscardL=",indexLp1,discard


deallocate(valuework)
deallocate(coeffwork)
deallocate(quantabigLbuffer)
deallocate(valueindex)
deallocate(coeffbuffer)
deallocate(coeffresult)
if(logic_spinreversal/=0) then
	deallocate(szzeroindex)
end if
!	tmp=0
!	do while(tmp==0) 
!		call sleep(2)
!	end do
return
end subroutine

