module exit_mod
! the module is to exit the total program

	use kinds_mod
	use communicate

	implicit none
	private
	save

! !PUBLIC MEMBER FUNCTIONS:
	public :: exit_DMRG

! !DEFINED PARAMETERS:

	integer (int_kind), parameter, public :: &
		sigExit = 0, &     ! signal for normal exit
		sigAbort = -1      ! signal for aborting (exit due to error)

contains
!==============================================================

subroutine exit_DMRG(exit_mode, exit_message)

! !DESCRIPTION:
! This routine prints a message, exits any message environment
! and cleans up before stopping


! !INPUT PARAMETERS:
	implicit none
	integer (int_kind), intent(in) :: exit_mode ! method for exiting (normal exit or abort)
	character (*), intent(in) :: exit_message   ! message to print before stopping

	integer (int_kind) :: ierr ! error flag

! print message
	if (myid == master_id) then
		select case(exit_mode)
		case(sigExit)
			write (*,'(a14)') 'DMRG exiting...'
		case(sigAbort)
			write (*,'(a15)') 'DMRG aborting...'
		case default
			write (*,'(a37)') 'DMRG exiting with unknown exit mode...'
		end select

            write (*,*) exit_message
	endif
!-----------------------------------------------------------------------
!
! exit or abort the message-passing environment if required
!
!-----------------------------------------------------------------------

	select case(exit_mode)
	case(sigExit)
		call exit_message_environment(ierr)
	case(sigAbort)
		call abort_message_environment(ierr)
	case default
	end select

	stop
return

end subroutine exit_DMRG

!========================================================================

end module exit_mod
