program main
	! This is a DMRG_MPS program quantum chemistry
	! only PPP model has been written

	USE variables
	use communicate
	use exit_mod
	use mpi, only : MPI_WTIME

	implicit none
	! local 
	real(kind=r8) :: starttime,endtime

	call init_communicate
	starttime=MPI_WTIME()
	
	if(nprocs<2) then
		call exit_DMRG(sigAbort,"nprocs<2 failed!")
	end if

	! read the input files
	call ReadInput

	! allocate the operator to different process
	call LoadBalance
	
	! do infinit DMRG process
	call Infinit_MPS

	! do finit DMRG process
	call Finit_MPS
	
	! calculate transition moment between gs and ex if nstate/=1
	if(nstate/=1) then
		call TransMoment
	end if
	
	! count if the matrix operamat is sparse or not
	call countnonzero

	endtime=MPI_WTIME()
	call master_print_message(endtime-starttime,"RUMTIME:")
	call exit_DMRG(0,"Program DMRG-X end successfully")

end program main
