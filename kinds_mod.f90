Module kinds_mod
! define number kinds

	implicit none
	private
	save

! define parameters

	integer,parameter,public :: &
	char_len = 100 , &
	log_kind = kind (.true.) , &
	int_kind = kind(1) , &
	i1 = selected_int_kind(2) , &
	i2 = selected_int_kind(4) , &
	i4 = selected_int_kind(6) , &
	i8 = selected_int_kind(13) , &
	r4 = selected_real_kind(6) , &
	r8 = selected_real_kind(13) 


end module kinds_mod
