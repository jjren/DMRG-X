subroutine Hamiltonian(direction)
! this is the core subroutine
! contruct the hamiltonian on the fly and muliply it with coefficient of
! the wavefunction and transfter it to the 0 process

	USE variables, only : nstate,sweeps,isweep
	use communicate
	use mpi

	implicit none

	integer :: &
		lim   ,  &       ! the expanding small space in davidson iteration
		ilow  ,  &       ! index of lowest eigenpair 
		ihigh ,  &       ! index of the highest eigenpair
		niv   ,  &       ! number of initial vector
		mblock,  &       ! number of vector to be targeted in each iteration
		maxiter          ! maxiter iterations
	real(kind=8) :: &
		crite ,   &      ! convergence eigenvalue
		critc ,   &      ! convergence coefficiency
		critr ,   &      ! convergence residual vector norm
		ortho            ! orthogonality threshold
	integer,allocatable :: iselec(:)   ! selected eigenpair not used
	integer :: error
	character(len=1) :: direction
	! i,l,r direction l is L space to R space sweep
	
	call master_print_message("enter hamiltonian subroutine")

	! default value
	ilow=1
	ihigh=nstate
	lim=nstate+20    ! this number 20 can be changed consider the efficiency
	niv=nstate
	mblock=nstate
	maxiter=800
	allocate(iselec(lim),stat=error)
	if(error/=0) stop
	iselec=-1

	! in the initial few sweeps
	if(1.0D-5*(1.0D-1)**isweep>1.1D-10) then
		crite=1.0D-5*(1.0D-1)**isweep
		critc=1.0D-5*(1.0D-1)**isweep
		critr=1.0D-5*(1.0D-1)**isweep
		ortho=1.0D-4*(1.0D-1)**isweep
	else
		crite=1.0D-9
		critc=1.0D-9
		critr=1.0D-9
		ortho=1.0D-8
	end if

	! in the last 3 sweeps
	if(isweep==sweeps .or. isweep==sweeps-1 .or. &
	isweep==sweeps-2) then
		crite=1.0D-10
		critc=1.0D-10
		critr=1.0D-10
		ortho=1.0D-9
	end if

	call Davidson_Wrapper(direction,lim,ilow,ihigh,iselec,niv,mblock,crite,critc,critr,ortho,maxiter)
	
	deallocate(iselec)
return

end subroutine Hamiltonian

