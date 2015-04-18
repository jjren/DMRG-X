subroutine store_operatorL(index1)
! this subroutine is to store the operator of every site(L space)
! in fact only the 4M basis operator matrix need to be store
! index1 is the sigmaL index
	use variables
	use mpi

	implicit none

	integer :: index1
	integer :: i,reclength,operaindex
	character(len=50) :: filename
	logical :: alive
	integer :: thefile,status(MPI_STATUS_SIZE)
	integer(kind=MPI_OFFSET_KIND) :: offset

	! reclength is the length of direct io
	! ifort use 4byte as 1 by default
	
	if(myid==0) then
		write(*,*) "enter in store_operatorL subroutine"
	end if
	
	! L+sigmaL space
! parallel io
	write(filename,'(i5.5,a8)') index1,'left.tmp'
	call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,thefile,ierr)
	do i=1,index1,1
	if(myid==orbid(i)) then
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if
		offset=(i-1)*16*subM*subM*3*8
		call MPI_FILE_WRITE_AT(thefile,offset,operamatbig(:,:,3*(operaindex-1)+1:3*operaindex),16*subM*subM*3,mpi_real8,status,ierr)
	end if
	end do
	call MPI_FILE_CLOSE(thefile,ierr)
		
	if(myid==0) then
		reclength=2*16*subM*subM
!----------------open a binary file-------------
		inquire(file="0-left.tmp",exist=alive)
		if(alive) then
			open(unit=101,file="0-left.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=101,file="0-left.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		
!-----------------------------------------------
		write(101,rec=index1) Hbig(:,:,1)
		close(101)
!------------------------------write the parity matrix---
		if(logic_spinreversal/=0) then
		reclength=2*subM
		! we use integer kind=2 in this symmlink matrix
		inquire(file="symmlink-left.tmp",exist=alive)
		if(alive) then
			open(unit=103,file="symmlink-left.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=103,file="symmlink-left.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		write(103,rec=index1) symmlinkbig(:,1,1)
		close(103)
		end if
!------------------------------write the quantabigL(4*subM,2)---
		reclength=4*subM*2
		! quantabigL is  integer(kind=4) so without *2
		inquire(file="quantabigL.tmp",exist=alive)
		if(alive) then
			open(unit=107,file="quantabigL.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=107,file="quantabigL.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		write(107,rec=index1) quantabigL
		close(107)
	end if


return

end subroutine store_operatorL



subroutine store_operatorR(index2)
! this subroutine is to store the operator of every site(R space)
! in fact only the 4M basis operator matrix need to be store
! index2 is the sigma R index
	use variables
	use mpi

	implicit none

	integer :: index2
	integer :: i,reclength,operaindex
	character(len=50) :: filename
	logical :: alive
	integer :: thefile,status(MPI_STATUS_SIZE)
	integer(kind=MPI_OFFSET_KIND) ::offset

	if(myid==0) then
		write(*,*) "enter in store_operatorR subroutine"
	end if

	! reclength is the length of direct io
	! ifort use 4byte as 1 by default

	write(filename,'(i5.5,a9)') index2,'right.tmp'
	call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,thefile,ierr)
	do i=norbs,index2,-1
	if(myid==orbid(i)) then
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if
		offset=(norbs-i)*16*subM*subM*3*8
		call MPI_FILE_WRITE_AT(thefile,offset,operamatbig(:,:,3*(operaindex-1)+1:3*operaindex),16*subM*subM*3,mpi_real8,status,ierr)
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
			open(unit=102,file="0-right.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
!-----------------------------------------------
! sigmaR index =norbs is 1; norbs-1 is 2
		write(102,rec=norbs+1-index2) Hbig(:,:,2)
		close(102)
!------------------------------write the parity matrix---
		if(logic_spinreversal/=0) then
		reclength=2*subM
		inquire(file="symmlink-right.tmp",exist=alive)
		if(alive) then
			open(unit=104,file="symmlink-right.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=104,file="symmlink-right.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		write(104,rec=norbs+1-index2) symmlinkbig(:,1,2)
		close(104)
		end if
!------------------------------write the quantabigR(4*subM,2)---
		reclength=4*subM*2
		! quantabigR is  integer(kind=4) so without *2
		inquire(file="quantabigR.tmp",exist=alive)
		if(alive) then
			open(unit=108,file="quantabigR.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=108,file="quantabigR.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		write(108,rec=norbs+1-index2) quantabigR
		close(108)
	end if


return

end subroutine store_operatorR
