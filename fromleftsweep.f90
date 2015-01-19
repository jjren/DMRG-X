subroutine fromleftsweep
! this subroutine is finit MPS procedure
! left space is the system
	USE mpi
	USE variables

	implicit none
	integer :: tmp

	if(myid==0) then
		write(*,*) "enter in subroutine fromleftsweep"
	end if
	
	if(nleft/=exactsite) then
		call onesitematrix(nleft+1)
		call system_bigL
		call system_constructquantaL
		call store_operatorL(nleft+1)
	!	if(myid==0) then
	!		read(*,*) tmp
	!	end if
	else
		call enviro_bigL
	end if

	call enviro_bigR
	!if(logic_spinreversal/=0) then
	!	call Spin_reversalmatL
	!end if
	!call fullmat
	call hamiltonian('l')
	call Renormalization(nleft+1,norbs-nright,'l')

return

end subroutine
