module communicate
	
	use kinds_mod
	use mpi

	implicit none
	private
	save

! !PUBLIC MEMBER FUNCTIONS:
	public :: init_communicate, &
                exit_message_environment, &
                abort_message_environment, &
                master_print_message

	interface master_print_message
		module procedure master_print_str
		module procedure master_print_str_dbl
		module procedure master_print_str_int
	end interface
	
	! MPI part
	integer(kind=int_kind),public :: &
	  myid,nprocs,master_id

	contains
!=============================================================

subroutine init_communicate
	
	implicit none
	! local variables
      integer (int_kind) :: ierr ! MPI error flag
	integer(kind=int_kind) :: version,subversion

	call MPI_INIT(ierr)

	master_id = 0
	call MPI_COMM_RANK (MPI_COMM_WORLD, myid, ierr)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ierr)
	call MPI_GET_VERSION(VERSION,SUBVERSION,ierr)
	write(*,*) "myid=",myid,"version=",version,"subversion=",subversion

end subroutine init_communicate

!=============================================================
!=============================================================

subroutine exit_message_environment(ierr)
! This routine exits the message environment properly when model
! stops.
	implicit none
	integer (int_kind), intent(out) :: ierr ! MPI error flag
   
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_FINALIZE(ierr)
return
end subroutine exit_message_environment

!=============================================================
!=============================================================

subroutine abort_message_environment(ierr)

! This routine aborts the message environment when model stops.
! It will attempt to abort the entire MPI COMM WORLD.
	implicit none
	integer (int_kind), intent(out) :: ierr ! MPI error flag
	
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	call MPI_FINALIZE(ierr)
return
end subroutine abort_message_environment

!=============================================================
!=============================================================

subroutine master_print_str(message)
	implicit none
	character(*), intent(in) :: message 
      integer (int_kind) :: ierr ! MPI error flag
	
	if(myid == master_id) then 
		write(*,*) "===========  ", message, "  =============="
	end if 
return
end subroutine master_print_str

!=============================================================
!=============================================================

subroutine master_print_str_dbl(num_p, message)
   
	implicit none 
	real(r8), intent(in) :: num_p
	character(*), intent(in) :: message
	integer (int_kind) :: ierr ! MPI error flag
	
	if(myid == master_id) then 
		write(*,*) "_________________________________"
		write(*,*)  message, num_p
		write(*,*) "---------------------------------"
	end if 
return
end subroutine master_print_str_dbl

!=============================================================
!=============================================================

subroutine master_print_str_int(num_p, message)
	implicit none
	integer(int_kind), intent(in) :: num_p
	character(*), intent(in) ::  message 
	integer (int_kind) :: ierr ! MPI error flag

	if(myid == master_id) then 
		write(*,*) "====== ",message, num_p, " ======"
	end if 
return
end subroutine  master_print_str_int

!=============================================================
!=============================================================


end module communicate

