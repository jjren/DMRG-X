program main
! This is a DMRG_MPS program quantum chemistry
! only PPP model has been written

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Jiajun Ren                      %
!% Zhigang Shuai's Group           %
!% Tsinghua University             %
!% Email: jiajunren0522@126.com    %   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	USE variables
	use communicate
	use exit_mod
	use mpi
	use MeanField
	use analysis_mod

	implicit none
	
	real(kind=r8) :: starttime,endtime
	integer :: ierr
	
	call init_communicate
	starttime=MPI_WTIME()
	
	if(nprocs<2) then
		call exit_DMRG(sigAbort,"nprocs<2 failed!")
	end if
	
	! read the input files
	call ReadInput
	
	! SCF mean field procedure
	if(logic_meanfield==1 .and. myid==0) then
		call SCFMain
	end if
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	
	! allocate the operator to different process
	call LoadBalance
	
	! do infinit DMRG process
	call Infinit_MPS
	
	! do finit DMRG process
	call Finit_MPS
	
	! do wave function analysis
	call Analysis

	! print the sparse info
	call countnonzero

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	endtime=MPI_WTIME()
	call master_print_message(endtime-starttime,"RUMTIME:")
	call exit_DMRG(0,"Program DMRG-X end successfully")

end program main
