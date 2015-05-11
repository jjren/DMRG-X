module analysis
! this subroutine is to anaylsis the wavefuntion and output the analysis result

	use variables
	use module_sparse
	use module_para

	implicit none
	private 
	save

	integer :: orbidnew(norbs),operalocate(norbs,norbs)  
	! orbidnew including 0 process; all the process is the same 
	! operalocate is the one partical term location on the specific process

	contains
!===================================================================
!===================================================================

subroutine OneRDM
! this subroutine is to calculate the one partical density matrix
! a(i,sigam)^+ aj(j,sigma)
! this subroutine is called after converged
! k^2*M*2 storage
! in PPP model; only the bond order matrix is needed







! parity number operator expectation
	do i=1,norbs,1
		
		

! bond order
