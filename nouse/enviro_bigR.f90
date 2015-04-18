Subroutine enviro_bigR
! this subroutine is to contruct the enviroment R+sigmaR space
! all the operamatbig(does not need operamatsma) read from disc
	USE mpi
	USE variables
	use communicate

	implicit none

	integer :: i,ierr
	integer :: reclength,operaindex
	logical :: alive
	character(len=50) :: filename
	integer :: thefile,status(MPI_STATUS_SIZE)
	integer(kind=MPI_OFFSET_KIND) ::offset

	! reclength is the length of direct io
	! ifort use 4 byte as 1 by default
	if(myid==0) then
		write(*,*) "enter in subroutine enviro_bigR"
	end if

	write(filename,'(i5.5,a9)') norbs-nright,'right.tmp'
	call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,thefile,ierr)
	do i=norbs,norbs-nright,-1
	if(myid==orbid(i)) then
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if
		offset=(norbs-i)*16*subM*subM*3*8
		call MPI_FILE_READ_AT(thefile,offset,operamatbig(:,:,3*(operaindex-1)+1:3*operaindex),16*subM*subM*3,mpi_real8,status,ierr)
	end if
	end do
	call MPI_FILE_CLOSE(thefile,ierr)

	if(myid==0) then
		reclength=2*16*subM*subM
!----------------open a binary file-------------
		inquire(file="0-right.tmp",exist=alive)
		if(alive) then
			open(unit=102,file="0-right.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) " 0-right.tmp doesn't exist!"
			stop
		end if
!-----------------------------------------------
! norbs is 1; norbs-1 is 2
		read(102,rec=nright+1) Hbig(:,:,2)
		close(102)
!------------------------------read the symmlink---
		if(logic_spinreversal/=0) then
		reclength=2*subM
		inquire(file="symmlink-right.tmp",exist=alive)
		if(alive) then
			open(unit=104,file="symmlink-right.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) "symmlink-right doesn't exist!"
			stop
		end if
		read(104,rec=nright+1) symmlinkbig(:,1,2)
		close(104)
		end if
!------------------------------read the quantabigR(4*subM,2)---
! every process need to read quantabigR

		reclength=4*subM*2
		! quantabigR is  integer(kind=4) so without *2
		inquire(file="quantabigR.tmp",exist=alive)
		if(alive) then
			open(unit=108,file="quantabigR.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) "quantabigR.tmp doesn't exist!"
			stop
		end if
		read(108,rec=nright+1) quantabigR
		close(108)
	end if

	call MPI_bcast(quantabigR,4*subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
return

end Subroutine
