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

	! reclength is the length of direct io
	! ifort use 4byte as 1 by default

	
	! L+sigmaL space
	do i=1,index1,1
	if(myid==orbid(i)) then
		reclength=2*16*subM*subM*3

!----------------open a binary file-------------
		write(filename,'(i5.5,a8)') index1,'left.tmp'
		inquire(file=trim(filename),exist=alive)
		if(alive) then
			open(unit=99,file=trim(filename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=99,file=trim(filename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
!-----------------------------------------------
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if

		write(99,rec=i) operamatbig(:,:,3*(operaindex-1)+1:3*operaindex)
		close(99)
	end if
	end do
		
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
		reclength=2*16*subM*subM
		inquire(file="paritymat-left.tmp",exist=alive)
		if(alive) then
			open(unit=103,file="paritymat-left.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=103,file="paritymat-left.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		write(103,rec=index1) adaptedbig(:,:,1)
		close(103)
		end if
!------------------------------write the quantabigL(4*subM,2)---
		reclength=4*subM*2
		! quantabigL is  integer(kind=4) so without *2
		inquire(file="quantabigL.tmp",exist=alive)
		if(alive) then
			open(unit=107,file="quantabigL.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=107,file="qunatabigL.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
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

	! reclength is the length of direct io
	! ifort use 4byte as 1 by default


	do i=norbs,index2,-1
	if(myid==orbid(i)) then
		reclength=2*16*subM*subM*3
!----------------open a binary file-------------
		write(filename,'(i5.5,a9)') index2,'right.tmp'
		inquire(file=trim(filename),exist=alive)
		if(alive) then
			open(unit=100,file=trim(filename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=100,file=trim(filename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
!-----------------------------------------------
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if
! sigmaR index = norbs is 1; norbs-1 is 2.....
		write(100,rec=norbs+1-i) operamatbig(:,:,3*(operaindex-1)+1:3*operaindex)
		close(100)
	end if
	end do
		
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
		reclength=2*16*subM*subM
		inquire(file="paritymat-right.tmp",exist=alive)
		if(alive) then
			open(unit=104,file="paritymat-right.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=104,file="paritymat-right.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		write(104,rec=norbs+1-index2) adaptedbig(:,:,2)
		close(104)
		end if
!------------------------------write the quantabigR(4*subM,2)---
		reclength=4*subM*2
		! quantabigR is  integer(kind=4) so without *2
		inquire(file="quantabigR.tmp",exist=alive)
		if(alive) then
			open(unit=108,file="quantabigR.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=108,file="qunatabigR.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		write(108,rec=norbs+1-index2) quantabigR
		close(108)
	end if


return

end subroutine store_operatorR
