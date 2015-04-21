Subroutine Enviro_Big(domain)
! this subroutine is to contruct the enviroment R+sigmaR/L+sigmaL space
! all the operamatbig(does not need operamatsma) read from disc

	USE mpi
	USE variables
	use communicate
	use exit_mod

	implicit none

	character(len=1) :: domain

	integer :: i
	integer :: orbnow,orbstart,orbend,orbref
	integer :: reclength,operaindex,recindex,Hindex
	character(len=50) :: filename,Hfilename,symmfilename,quantafilename
	logical :: alive

	integer :: thefile,status(MPI_STATUS_SIZE)  ! MPI flag
	integer(kind=MPI_OFFSET_KIND) ::offset
	integer :: ierr

	! reclength is the length of direct io
	! ifort use 4 byte as 1 by default
	call master_print_message("enter in subroutine Enviro_Big")
	
	if(domain=='L') then
		orbnow=nleft+1
		write(filename,'(i5.5,a8)') orbnow,'left.tmp'
		orbstart=1
		orbend=nleft+1
		orbref=1
	else if (domain=='R') then
		orbnow=norbs-nright
		write(filename,'(i5.5,a9)') orbnow,'right.tmp'
		orbstart=norbs-nright
		orbend=norbs
		orbref=norbs
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/='R' failed!")
	end if
	
	! L/R + sigmaL/R space

	call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,thefile,ierr)
	if(myid/=0) then
		do i=orbstart,orbend,1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			offset=abs(i-orbref)*16*subM*subM*3*8
			call MPI_FILE_READ_AT(thefile,offset,operamatbig(1,1,3*operaindex-2),16*subM*subM*3,mpi_real8,status,ierr)
		end if
		end do
	! 0 process
	else if(myid==0) then
		if(domain=='L') then
			Hfilename="0-left.tmp"
			symmfilename="symmlink-left.tmp"
			quantafilename="quantabigL.tmp"
			Hindex=1
		else if(domain=='R') then
			Hfilename="0-right.tmp"
			symmfilename="symmlink-right.tmp"
			quantafilename="quantabigR.tmp"
			Hindex=2
		end if
		recindex=abs(orbnow-orbref)+1

!----------------open a binary file-------------
		reclength=2*16*subM*subM
		inquire(file=trim(Hfilename),exist=alive)
		if(alive) then
			open(unit=102,file=trim(Hfilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) trim(Hfilename),"doesn't exist!"
			stop
		end if
!-----------------------------------------------
! norbs is 1; norbs-1 is 2
		read(102,rec=recindex) Hbig(:,:,Hindex)
		close(102)
!------------------------------read the symmlink---
		if(logic_spinreversal/=0) then
			reclength=2*subM
			inquire(file=trim(symmfilename),exist=alive)
			if(alive) then
				open(unit=104,file=trim(symmfilename),access="Direct",form="unformatted",recl=reclength,status="old")
			else
				write(*,*) trim(symmfilename),"doesn't exist!"
				stop
			end if
			read(104,rec=recindex) symmlinkbig(:,1,Hindex)
			close(104)
		end if
!------------------------------read the quantabigR/L(4*subM,2)---
! every process need to read quantabigR/L

		reclength=4*subM*2
		! quantabigR/L is  integer(kind=4) so without *2
		inquire(file=trim(quantafilename),exist=alive)
		if(alive) then
			open(unit=108,file=trim(quantafilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) trim(quantafilename),"doesn't exist!"
			stop
		end if
		if(domain=='L') then
			read(108,rec=recindex) quantabigL
		else if(domain=='R') then
			read(108,rec=recindex) quantabigR
		end if
		close(108)
	end if

	! broadcast the quantabigR/L to other process
	if(domain=='L') then
		call MPI_bcast(quantabigL,4*subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
	else if(domain=='R') then
		call MPI_bcast(quantabigR,4*subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
	end if

	call MPI_FILE_CLOSE(thefile,ierr)
return

end Subroutine Enviro_Big
