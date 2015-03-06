subroutine C2_resort
! this program is to resort the orbital and integral to combine to
! fullfill C2 symmetry(make the orbital in different irrep)
	use mpi
	use variables

	implicit none
	integer :: i,error
	real(kind=8),allocatable :: tsort(:,:),hubbardUsort(:),pppVsort(:,:)
	integer(kind=4),allocatable :: bondlinksort(:,:)
	
	allocate(tsort(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(hubbardUsort(norbs),stat=error)
	allocate(hubbardUsort(norbs),stat=error)
	if(error/=0) stop
	do i=1,norbs,1
		



return
end subroutine
