subroutine Store_Operator(domain)
! this subroutine is to store the operator of every site L/R space
! in fact only the 4M basis operator matrix need to be store
	use variables
	use mpi
	use communicate
	use exit_mod

	implicit none

	character(len=1) :: domain
	! local
	integer :: i
	integer :: reclength,operaindex,Hindex,recindex
	integer :: orbstart,orbend,orbref,orbnow
	! orbref is the reference orbindex 1 or norbs
	! orbnow is nleft+1 or norbs-nright
	character(len=50) :: filename,Hfilename,symmfilename,quantafilename
	logical :: alive

	integer :: thefile , status(MPI_STATUS_SIZE) ! MPI_flag
	integer(kind=MPI_OFFSET_KIND) :: offset
	integer :: ierr

	! in fortran io, reclength is the length of direct io
	! ifort use 4byte as 1 by default
	
	call master_print_message("enter in store_operator subroutine")
	
	! L+sigmaL/R+sigmaR space operator
	! MPI parallel io

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

	call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,thefile,ierr)
	do i=orbstart,orbend,1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			offset=abs(i-orbref)*16*subM*subM*3*8
			call MPI_FILE_WRITE_AT(thefile,offset,operamatbig(1,1,3*operaindex-2),16*subM*subM*3,mpi_real8,status,ierr)
		end if
	end do
	call MPI_FILE_CLOSE(thefile,ierr)

!=====================================================

	if(myid==0) then
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
			open(unit=101,file=trim(Hfilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=101,file=trim(Hfilename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
!-----------------------------------------------
		! sigmaR index =norbs is 1; norbs-1 is 2
		write(101,rec=recindex) Hbig(:,:,Hindex)
		close(101)
!------------------------------write the symmetry matrix---
		if(logic_spinreversal/=0) then
			reclength=2*subM
			! we use integer kind=2 in this symmlink matrix
			inquire(file=trim(symmfilename),exist=alive)
			if(alive) then
				open(unit=103,file=trim(symmfilename),access="Direct",form="unformatted",recl=reclength,status="old")
			else
				open(unit=103,file=trim(symmfilename),access="Direct",form="unformatted",recl=reclength,status="replace")
			end if
			write(103,rec=recindex) symmlinkbig(:,1,Hindex)
			close(103)
		end if
!------------------------------write the quantabigL/R(4*subM,2)---
		reclength=4*subM*2
		! quantabigL is  integer(kind=4) so without *2
		inquire(file=trim(quantafilename),exist=alive)
		if(alive) then
			open(unit=107,file=trim(quantafilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=107,file=trim(quantafilename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		if(domain=='L') then
			write(107,rec=recindex) quantabigL
		else if(domain=='R') then
			write(107,rec=recindex) quantabigR
		end if
		close(107)
	end if

return

end subroutine Store_Operator


