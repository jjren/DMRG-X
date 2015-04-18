subroutine fromleftsweep
! this subroutine is finit MPS procedure
! left space is the system
	USE mpi
	USE variables

	implicit none

	if(myid==0) then
		write(*,*) "enter in subroutine fromleftsweep"
	end if
	
	if(nleft/=exactsite) then
		call onesitematrix(nleft+1)
		call system_bigL
		call system_constructquantaL
		call store_operatorL(nleft+1)
	else
		call enviro_bigL
	end if
	
	if(nleft==norbs/2-1 .and. logic_C2/=0) then
		call C2_copy('l')
	else
		call enviro_bigR
	end if
	!if(logic_spinreversal/=0) then
	!	call Spin_reversalmatL
	!end if
!	call fullmat
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
