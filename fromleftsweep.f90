subroutine fromleftsweep
	USE mpi
	USE variables

	implicit none
	
	call onesitematrix(nleft+1)
	call system_bigL
	call enviro_bigR



return

end subroutine
