MODULE InitialGuess
	use variables
	use communicate
	use kinds_mod
	use module_sparse
	implicit none
! this module contains the subroutine that generate the initial guess of the MPS
! procedure

	contains
!===========================================================

subroutine InitialStarter(direction,lvector,nvector,initialcoeff)
! get the intital vector
! here we only consider
! 1. nstate=1 finit
! 2. nstate>1 finit exscheme=1
! 3  infinit
! we can add exscheme=2 and later
	use mathlib
	use symmetry
	implicit none

	character(len=1) :: direction
	integer :: lvector,nvector
	real(kind=r8) :: initialcoeff(lvector*nvector),norm(nvector)
	real(kind=r8),allocatable :: nosymmguess(:)
	integer :: i,error
	
	if(nvector/=nstate) then
		call master_print_message(nvector,"nvector/=nstate")
		stop
	end if

	if(direction/='i' .and. logic_C2==0 .and. &
	formernelecs==nelecs .and. ifopenperturbation==.false.) then

		allocate(nosymmguess(ngoodstates*nvector),stat=error)
		if(error/=0) stop

		if(nvector==1) then
			call SingleInitialFinite(nosymmguess,ngoodstates,direction)
		else
			call MoreInitialFinite(nosymmguess,ngoodstates,direction)
		end if

		if(logic_spinreversal/=0) then
			do i=1,nvector,1
				call symmetrizestate(ngoodstates,nosymmguess((i-1)*ngoodstates+1:i*ngoodstates),&
				initialcoeff((i-1)*lvector+1:i*lvector),'s')
			end do
		else
			call copy(nosymmguess,initialcoeff)
		end if

		deallocate(nosymmguess)
		
		call GramSchmit(nvector,lvector,initialcoeff,norm)
		write(*,*) "the Initial guessvector norm",norm
	else
		call InitialRandom(initialcoeff,nvector,lvector)
		call GramSchmit(nvector,lvector,initialcoeff,norm)
	!	write(*,*) "the Initial guessvector norm",norm

	!	there are some problem(such as local minimun) when us univector
	!	especically in non interacting system
		call master_print_message("use the algorithm default as initialguess")
		! univector in davidson
		nvector=0  ! niv=0  
	end if
	

return

end subroutine InitialStarter


!===========================================================
	subroutine MoreInitialFinite(guesscoeff,lvector,direction)
	! the new Initial Guess can support nstate/=1
	! the algrithom can be read from the manual
	! the guesscoeff could be approximated as U+CB though UU+/=I

	USE BLAS95
	USE F95_PRECISION

	implicit none
	include "mkl_spblas.fi"
	integer :: lvector
	real(kind=r8) :: guesscoeff(lvector*nstate)
	real(kind=r8),allocatable :: leftu(:,:),rightv(:,:)&
	,LRcoeff(:,:,:),LRcoeff1(:,:,:)
	logical :: alive
	integer :: reclength
	integer :: subMbefore
	integer :: error,i,j,m,k
	character(len=1) :: direction
	integer :: job(8),info

	write(*,*) "enter MoreInitialFinit subroutine"
	! two site dmrg
	if((nright+nleft+2)/=norbs) then
		write(*,*) "-----------------------------------"
		write(*,*) "two site dmrg nright+nleft+2/=norbs"
		write(*,*) "-----------------------------------"
		stop
	end if

	allocate(leftu(4*subM,subM),stat=error)
	if(error/=0) stop
	allocate(rightv(subM,4*subM),stat=error)
	if(error/=0) stop

	reclength=2*subM*subM

	inquire(file="wavefunction.tmp",exist=alive)
	if(alive) then
		open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
	else
		write(*,*) "wavefunction.tmp doesn't exist"
		stop
	end if
	if(direction=='l' .and. nleft>exactsite) then
		if(nleft==exactsite+1) then   ! the last step Lrealdim
			subMbefore=4**exactsite
		else
			subMbefore=subM
		end if

		do i=1,4,1
			read(105,rec=4*(nleft-1)+i) leftu((i-1)*subMbefore+1:i*subMbefore,1:subM)
		end do

		do i=1,4,1
			read(105,rec=4*(nleft+1)+i) rightv(1:subM,i:4*Rrealdim:4)
		end do
	else if(direction=='r' .and. nright>exactsite) then
		
		do i=1,4,1
			read(105,rec=4*nleft+i) leftu((i-1)*Lrealdim+1:i*Lrealdim,1:subM)
		end do
		
		if(nright==exactsite+1) then ! the last step Rrealdim
			subMbefore=4**exactsite
		else
			subMbefore=subM
		end if
		
		do i=1,4,1
			read(105,rec=4*(nleft+2)+i) rightv(1:subM,i:4*subMbefore:4)
		end do
	end if

	allocate(LRcoeff(4*subM,4*subM,nstate),stat=error)
	if(error/=0) stop
	allocate(LRcoeff1(4*subM,4*subM,nstate),stat=error)
	if(error/=0) stop
	
	job(1)=1
	job(2)=1
	job(3)=1
	job(4)=2
	job(5)=0
	job(6)=1
	! transfer coeffIF to dense matrix
	do i=1,nstate,1
		call mkl_ddnscsr(job,4*subM,4*subM,LRcoeff1(:,:,i),4*subM,coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),info)
		if(info/=0) then
			call master_print_message(info,"MoreInitFinite info/=")
			stop
		end if
	end do

	if(direction=='l' .and. nleft>exactsite) then
		if(Lrealdim/=subM) then
			write(*,*) "-----------------------------------"
			write(*,*) "newguessfinit, Lrealdim/=subM failed!",Lrealdim,subM
			write(*,*) "-----------------------------------"
			stop
		end if

		do i=1,nstate,1
			call gemm(leftu(1:4*subMbefore,:),LRcoeff1(1:4*subMbefore,:,i),LRcoeff(1:subM,:,i),'T','N',1.0D0,0.0D0)   ! do U+*C
			LRcoeff1(:,:,i)=0.0D0
			do j=1,subM,1
				do k=1,4,1
				LRcoeff1(subM*(k-1)+1:subM*k,j,i)=LRcoeff(1:subM,(j-1)*4+k,i)   ! tranfer U+C to new form
				end do
			end do
			call gemm(LRcoeff1(:,1:subM,i),rightv,LRcoeff(:,1:4*Rrealdim,i),'N','N',1.0D0,0.0D0)  !U+C*V
		end do

	else if(direction=='r' .and. nright>exactsite) then
		if(Rrealdim/=subM) then
			write(*,*) "-----------------------------------"
			write(*,*) "newguessfinit, Rreadlim/=subM failed!"
			write(*,*) "-----------------------------------"
			stop
		end if
		do i=1,nstate,1
			call gemm(LRcoeff1(:,1:4*subMbefore,i),rightv(:,1:4*subMbefore),LRcoeff(:,1:subM,i),'N','T',1.0D0,0.0D0)
			LRcoeff1(:,:,i)=0.0D0
			do j=1,subM,1
				do k=1,4,1
					LRcoeff1(1:subM,(j-1)*4+k,i)=LRcoeff((k-1)*subM+1:k*subM,j,i)
				end do
			end do
			call gemm(leftu,LRcoeff1(1:subM,:,i),LRcoeff(1:4*Lrealdim,:,i),'N','N',1.0D0,0.0D0)
		end do

	else if((direction=='l' .and. nleft==exactsite) .or. (direction=='r' .and. nright==exactsite)) then
		LRcoeff=LRcoeff1  ! directly use the last step result
	end if

	m=0
	do k=1,nstate,1
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
			quantabigL(j,2)+quantabigR(i,2)==totalSz) then
			m=m+1
			guesscoeff(m)=LRcoeff(j,i,k)
		end if
	end do
	end do
	end do

	if(m/=lvector*nstate) then
		write(*,*) "----------------------------------------------"
		write(*,*) "guesscoeff good quantum states number wrong! failed!"
		write(*,*) "----------------------------------------------"
		stop
	end if


	deallocate(leftu)
	deallocate(rightv)
	deallocate(LRcoeff1)
	deallocate(LRcoeff)
	end subroutine MoreInitialFinite

!===========================================================

!===========================================================
	Subroutine SingleInitialFinite(nosymmguess,lvector,direction)
	! when nstate=1 then we can use the last stored matrix to
	! contruct the Initial Guess
	! only used in the finit MPS and nstate=1

	!out :: nosymmguess(lvector) :: the finit initial guess
	!in  :: lvector :: ngoodstates,the nnosymmstate
	!in  :: direction :: the process of DMRG
	
	USE BLAS95
	USE F95_PRECISION

	implicit none
	include "mkl_spblas.fi"
	integer :: lvector
	real(kind=r8) :: nosymmguess(lvector)
	character(len=1) :: direction
	real(kind=r8),allocatable :: leftu(:,:),rightv(:,:),singularvalue(:)&
	,LRcoeff(:,:)
	logical :: alive
	integer :: reclength
	integer :: error,i,j,m
	real(kind=r8) :: norm
	integer :: job(8),info
	
	call master_print_message("enter SingleInitialFinite subroutine")
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
	end if

		
	allocate(LRcoeff(4*Lrealdim,4*Rrealdim),stat=error)
	if(error/=0) stop
	
	if((direction=='l' .and. nleft==exactsite) .or. (direction=='r' .and. nright==exactsite)) then
		! direct use the last step result
		if(nstate==1) then
		!	LRcoeff=coeffIF(1:4*Lrealdim,1:4*Rrealdim,1)
			job(1)=1
			job(2)=1
			job(3)=1
			job(4)=2
			job(5)=0
			job(6)=1
			call mkl_ddnscsr(job,4*Lrealdim,4*Rrealdim,LRcoeff,4*Lrealdim,coeffIF(:,1),coeffIFcolindex(:,1),coeffIFrowindex(:,1),info)
			if(info/=0) then
				call master_print_message(info,"SingleInitFinit info/=0")
				stop
			end if
		else
			write(*,*) "================================="
			write(*,*) "garnet chan's specific algorithom"
			write(*,*) "================================="
			stop
		end if
	else
		! recombine the two site sigmaL sigmaR coefficient
		call gemm(leftu,rightv,LRcoeff,'N','N',1.0D0,0.0D0)
	end if

	
	m=0
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
			quantabigL(j,2)+quantabigR(i,2)==totalSz) then
			m=m+1
			nosymmguess(m)=LRcoeff(j,i)
		end if
	end do
	end do

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
return
end subroutine SingleInitialFinite
!================================================================

!================================================================

! Random Initial Guess in the infinit DMRG and nstate/=1 DMRG

Subroutine InitialRandom(guessvector,lvector,num)
!   in  :: lvector :: the dimension of the davidson Hamiltonian matrix
!   in  :: num     :: the number of target states
!   out :: guessvetor(lvector*num) :: the initial guessvetor workspace

	implicit none

	integer :: num,lvector
	real(kind=r8) :: guessvector(num*lvector),randomx
	integer :: i
	
	call master_print_message("enter InitialRandom subroutine")
	call random_seed()

	do i=1,num*lvector,1
		call random_number(randomx)
		guessvector(i)=randomx
	end do

return
end subroutine InitialRandom
!===============================================================



!================================================================

! unit vector that corresponding the smallest diagonal element!
! because we want to get the smallest eigenvalue
! but this method may not good in PPP model
! because initially the R space and L space does not have non-diagonal
! interation

	Subroutine Initialunitvector(HDIAG,guessvector,lvector,num)
	use mpi
	use variables
	implicit none

	integer :: num,lvector
	real(kind=r8) :: guessvector(num*lvector)
	real(kind=r8) :: HDIAG(lvector)
	integer :: i,j,k,l,minindex(num)
	real(kind=r8) :: mindiag(num)

	
	if(myid==0) then
		write(*,*) "enter Initialunitvector subroutine"
! the first is the smallest number
		minindex=0
		mindiag=1.0D10
		do i=1,lvector,1
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
			guessvector((i-1)*lvector+minindex(i))=1.0D0
		end do
	end if

	return
	end subroutine
!===============================================================

end MODULE
