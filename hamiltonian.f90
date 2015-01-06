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
! default value
	ilow=1
	ihigh=nstate
	lim=nstate+20
	niv=nstate
	mblock=nstate
	crite=1.0D-7
	critc=1.0D-7
	critr=1.0D-7
	ortho=1.0D-6
	maxiter=200
	allocate(iselec(lim),stat=error)
	if(error/=0) stop
	iselec=-1
!-----------------------

	call davidson_wrapper(direction,lim,ilow,ihigh,iselec,niv,mblock,crite,critc,critr,ortho,maxiter)

return

end subroutine hamiltonian







