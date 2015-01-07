Subroutine enviro_bigR
! this subroutine is to contruct the enviroment R+sigmaR space
! all the operamatbig(does not need operamatsma) read from disc
	USE mpi
	USE variables

	implicit none

	integer :: i
	integer :: reclength,operaindex
	logical :: alive
	character(len=50) :: filename

	! reclength is the length of direct io
	! ifort use 4 byte as 1 by default

	do i=norbs,norbs-nright,-1
	if(myid==orbid(i)) then
		reclength=2*16*subM*subM*3
!----------------open a binary file-------------
		write(filename,'(i5.5,a9)') norbs-nright,'right.tmp'
		inquire(file=trim(filename),exist=alive)
		if(alive) then
			open(unit=100,file=trim(filename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) trim(filname),"doesn't exist!"
			stop
		end if
!-----------------------------------------------
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if
! norbs is 1; norbs-1 is 2.....
		read(100,rec=norbs+1-i) operamatbig(:,:,3*(operaindex-1)+1:3*operaindex)
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
			write(*,*) " 0-right.tmp doesn't exist!"
			stop
		end if
!-----------------------------------------------
! norbs is 1; norbs-1 is 2
		read(102,rec=nright+1) Hbig(:,:,2)
		close(102)
!------------------------------read the parity matrix---
		if(logic_spinreversal/=0) then
		reclength=2*16*subM*subM
		inquire(file="paritymat-right.tmp",exist=alive)
		if(alive) then
			open(unit=104,file="paritymat-right.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) "paritymat-right doesn't exist!"
			stop
		end if
		read(104,rec=nright+1) adaptedbig(:,:,2)
		close(104)
		end if
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

return

end Subroutine
