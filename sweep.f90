subroutine Sweep(direction)
! this subroutine is finit MPS procedure, sweep process

	USE mpi
	USE variables
	use communicate
	use exit_mod

	implicit none
	
	character(len=1) :: direction
	
	! local 
	character(len=1) :: domain,envirodomain
	integer :: nsysorb, &      ! nsysorb=nleft or nright
	           orbnow          ! orbnow=nleft+1 or norbs-nright

	
	call master_print_message("enter in subroutine Sweep")

	if(direction=='l') then
		domain="L"
		envirodomain="R"
		nsysorb=nleft
		orbnow=nleft+1
	else if(direction=='r') then
		domain="R"
		envirodomain="L"
		nsysorb=nright
		orbnow=norbs-nright
	else
		call exit_DMRG(sigAbort,"Sweep subroutine direction/=l .and. /=r failed!")
	end if

	if(nsysorb/=exactsite) then
		call OnesiteMatrix(orbnow)
		call System_Big(domain)
		call System_Constructquanta(domain)
		call Store_Operator(domain)
	else
		call Enviro_Big(domain)
	end if
	
	if(nleft==nright .and. logic_C2/=0) then
		call C2_Copy(direction)
	else
		call Enviro_Big(envirodomain)
	end if
	call Hamiltonian(direction)

	if(isweep==sweeps .and. direction=='l' .and. nleft==(norbs-1)/2) then
		! this is to do the chan proposed trace excited algrithom
		call master_print_message("In the last step, did not do Renormalization")
	else
		call Renormalization(nleft+1,norbs-nright,direction)
	end if

return

end subroutine Sweep
