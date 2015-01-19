	program main
	! This is a DMRG_MPS program quantum chemistry

	USE MPI
	USE Variables

	implicit none
	real(kind=8) :: starttime,endtime,ticktime

	call MPI_INIT(ierr)
	starttime=MPI_WTIME()
	call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
	call MPI_GET_VERSION(VERSION,SUBVERSION,ierr)
	write(*,*) "myid=",myid,"version=",version,"subversion=",subversion
	call MPI_Barrier(MPI_COMM_WORLD,ierr)

	!read the input files
	call readinput
	! allocate the operator to different process
	call loadbalance
	! construct onesite parity matrix
	!if(myid==0 .and. logic_spinreversal/=0) then
	!	call parity_onesitematrix
	!end if
	! do infinit DMRG
	call infinit_MPS
	! do finit DMRG
	!call finit_MPS
  



	endtime=MPI_WTIME()
	ticktime=MPI_WTICK()
	call MPI_Barrier(MPI_COMM_WORLD,ierr)
	write(*,*) "myid=",myid,"totaltime=",endtime-starttime,"ticktime=",ticktime
	call MPI_FINALIZE(ierr)

	end program main
