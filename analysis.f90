subroutine analysis
! this subroutine is to anaylsis the wavefuntion and output the info

	use mpi
	use variables

	implicit none
	integer :: i
! do another loop to calculate all the operator

! parity number operator expectation
	do i=1,norbs,1
		
		

! bond order
