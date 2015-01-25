subroutine hamiltonian(direction)
! this is the core subroutine
! contruct the hamiltonian on the fly and muliply it with coefficient of
! the wavefunction and transfter it to the 0 process

	USE MPI
	USE variables

	implicit none
	integer :: lim,ilow,ihigh,niv,mblock,maxiter
	real(kind=8) :: crite,critc,critr,ortho
	integer,allocatable :: iselec(:)
	integer :: error
	character(len=1) :: direction
	
	if(myid==0) then
		write(*,*) "enter hamiltonian subroutine"
	end if

! default value
	ilow=1
	ihigh=nstate
	lim=nstate+20
	niv=nstate
	mblock=nstate
	
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


	if(isweep==sweeps .or. isweep==sweeps-1 .or. &
	isweep==sweeps-2) then
		crite=1.0D-10
		critc=1.0D-10
		critr=1.0D-10
		ortho=1.0D-9
	end if

	maxiter=200
	allocate(iselec(lim),stat=error)
	if(error/=0) stop
	iselec=-1
!-----------------------

	call davidson_wrapper(direction,lim,ilow,ihigh,iselec,niv,mblock,crite,critc,critr,ortho,maxiter)

	deallocate(iselec)
return

end subroutine hamiltonian







