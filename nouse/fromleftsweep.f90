subroutine fromleftsweep
! this subroutine is finit MPS procedure
! left space is the system
	USE mpi
	USE variables
	use communicate

	implicit none
	
	call master_print_message("enter in subroutine fromleftsweep")
	
	if(nleft/=exactsite) then
		call OnesiteMatrix(nleft+1)
		call System_Big('L')
		call System_Constructquanta('L')
		call Store_Operator('L')
	else
		call enviro_bigL
	end if
	
	if(nleft==norbs/2-1 .and. logic_C2/=0) then
		call C2_copy('l')
	else
		call enviro_bigR
	end if
	call hamiltonian('l')
	if(isweep==sweeps .and. nleft==(norbs+1)/2) then
		if(myid==0) then
! this is to do the chan proposed trace excited algrithom
			write(*,*) "In the last step, did not do Renormalization"
		end if
	else
		call Renormalization(nleft+1,norbs-nright,'l')
	end if

return

end subroutine
