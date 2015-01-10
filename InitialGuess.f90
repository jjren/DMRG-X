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
	character(len=1) :: direction

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
		read(105,rec=4*(norbs-nright-1)+i) rightv(1:subM,i:4*Rrealdim:4)
		end do

		open(unit=106,file="singularvalue.tmp",status="old")
		read(106,*) singularvalue(1:subM) 
		singularvalue=sqrt(singularvalue)

! be careful the finit initial guessvector
		if((direction=='l' .and. nleft>exactsite) .or.  &
			(direction=='r' .and. nleft==norbs-exactsite-2)) then
			if(Lrealdim/=subM) then
				write(*,*) "-----------------------------------"
				write(*,*) "guessfinit, Lreadlim/=subM failed!"
				write(*,*) "-----------------------------------"
				stop
			end if
			do i=1,subM,1
				leftu(i:4*subM:subM,:)=leftu(i:4*subM:subM,:)*singularvalue(i)
			end do
		else if((direction=='r' .and. nright>exactsite) .or.  &
			(direction=='l' .and. nright==norbs-exactsite-2)) then
			if(Rrealdim/=subM) then
				write(*,*) "-----------------------------------"
				write(*,*) "guessfinit, Rreadlim/=subM failed!"
				write(*,*) "-----------------------------------"
				stop
			end if
			do i=1,subM,1
				rightv(:,i:4*subM:subM)=rightv(:,i:4*subM:subM)*singularvalue(i)
			end do
		end if

		! recombine the two site sigmaL sigmaR coefficient
		
		allocate(LRcoeff(4*Lrealdim,4*Rrealdim),stat=error)
		if(error/=0) stop
		call gemm(leftu,rightv,LRcoeff,'N','N',1.0D0,0.0D0)

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

	implicit none

	integer :: num
	real(kind=8) :: guessvector(num*ngoodstates),randomx,norm
	integer :: i,j,k
	logical :: done

	
	if(myid==0) then
		write(*,*) "enter Initialrandomweight subroutine"
		
		norm=0.0D0
		call random_seed()

		do i=1,num*ngoodstates,1
			call random_number(randomx)
			guessvector(i)=randomx
		end do

! when L space Sz>0 then we need to make the symmetry pair using some coefficient
! when L space Sz=0 then we need to make the L+R space basis to have the specfic 
! spin parity. Then set others to zero
		if(logic_spinreversal/=0) then
			do i=1,ngoodstates,1
				if(quantabigL(symmlinkgood(i,1),2)>=0 .and. abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1)) then
					done=.false.
				do j=1,ngoodstates,1
					if(symmlinkgood(j,1)==abs(symmlinkbig(symmlinkgood(i,1),1,1)) &
						.and. symmlinkgood(j,2)==abs(symmlinkbig(symmlinkgood(i,2),1,2))) then
						guessvector(j:num*ngoodstates:ngoodstates)=guessvector(i:num*ngoodstates:ngoodstates)&
						*DBLE(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2)))&
						*DBLE(logic_spinreversal)
						done=.true.
						exit
					end if
				end do
					if(done/=.true.) then
						write(*,*) "-------------------------------------------------"
						write(*,*) "initialrandomweight spin reversal adapted failed!"
						write(*,*) "-------------------------------------------------"
						stop
					end if
				else if(quantabigL(symmlinkgood(i,1),2)==0 .and. &
				abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1)) then
					if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))/=logic_spinreversal) then
						guessvector(i:num*ngoodstates:ngoodstates)=0.0D0
					end if
				end if
			end do
		end if
		
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
