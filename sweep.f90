subroutine Sweep(direction)
! this subroutine is finit MPS procedure, sweep process

	USE mpi
	USE variables
	use communicate
	use exit_mod
	use Renormalization_mod
	use OnesiteMatrix
	use hamiltonian_mod
	use module_sparse

	implicit none
	
	character(len=1) :: direction
	
	! local 
	character(len=1) :: domain,envirodomain
	integer :: nsysorb, &      ! nsysorb=nleft or nright
	           orbnow          ! orbnow=nleft+1 or norbs-nright
	integer :: i
	
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
		call ConstructOnesiteMatrix(orbnow)
		call System_Big(domain)
		call System_Constructquanta(domain)
		call Store_Operator(domain)
	else
		call Enviro_Big(domain)
	end if
	
	if(nleft==nright .and. logic_C2/=0) then
		call C2_Copy(direction)
		! in the C2 symmetry, we calculate the two subspace
		! first the reverse space; then the input space
		do i=1,2,1
			logic_C2=logic_C2*(-1)
			call master_print_message(logic_C2,"logic_C2=")
			call Hamiltonian(direction)
			if(i==1 .and. myid==0) then
				coeffIF(:,nstate+1:2*nstate)=coeffIF(:,1:nstate)
				coeffIFcolindex(:,nstate+1:2*nstate)=coeffIFcolindex(:,1:nstate)
				coeffIFrowindex(:,nstate+1:2*nstate)=coeffIFrowindex(:,1:nstate)
			end if
		end do
	else
		call Enviro_Big(envirodomain)
		call Hamiltonian(direction)
	end if

	if(isweep==sweeps .and. direction=='l' .and. nleft==(norbs-1)/2) then
		! this is to do the chan proposed trace excited algrithom
		call master_print_message("In the last step, did not do Renormalization")
	else
		call Renormalization(direction)
	end if

return

end subroutine Sweep
