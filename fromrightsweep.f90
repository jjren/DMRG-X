subroutine fromrightsweep
! this subroutine is finit MPS procedure
! right space is the system
	USE mpi
	USE variables

	implicit none
	
	if(myid==0) then
		write(*,*) "enter in subroutine fromrightsweep"
	end if
	
	if(nright/=exactsite) then
		call onesitematrix(norbs-nright)
		if(logic_C2==0) then
			call system_bigR
			call system_constructquantaR
		else
			call system_bigRreverse
			call system_constructquantaRreverse
		end if
		call store_operatorR(norbs-nright)
	else
		call enviro_bigR
	end if

	if(nright==norbs/2-1 .and. logic_C2/=0) then
		call C2_copy('r')
	else
		call enviro_bigL
	end if
	!if(logic_spinreversal/=0) then
	!	call Spin_reversalmatR
	!end if
!	call fullmat
	call hamiltonian('r')
	call Renormalization(nleft+1,norbs-nright,'r')

return

end subroutine
