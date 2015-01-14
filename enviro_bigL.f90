Subroutine enviro_bigL
! this subroutine is to contruct the enviroment L+sigmaL space
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
	if(myid==0) then
		write(*,*) "enter in subroutine enviro_bigL"
	end if

	! L+sigmaL space
	do i=1,nleft+1,1
	if(myid==orbid(i)) then
		reclength=2*16*subM*subM*3

!----------------open a binary file-------------
		write(filename,'(i5.5,a8)') nleft+1,'left.tmp'
		inquire(file=trim(filename),exist=alive)
		if(alive) then
			open(unit=99,file=trim(filename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) trim(filename),"doesn't exist!"
			stop
		end if
!-----------------------------------------------
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if

		read(99,rec=i) operamatbig(:,:,3*(operaindex-1)+1:3*operaindex)
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
			write(*,*) "0-left.tmp doesn't exist!"
			stop
		end if
		
!-----------------------------------------------
		read(101,rec=nleft+1) Hbig(:,:,1)
		close(101)
!------------------------------read the symmlink---
		if(logic_spinreversal/=0) then
		reclength=2*subM
		inquire(file="symmlink-left.tmp",exist=alive)
		if(alive) then
			open(unit=103,file="symmlink-left.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) "symmlink-left.tmp doesn't exist!"
			stop
		end if
		read(103,rec=nleft+1) symmlinkbig(:,1,1)
		close(103)
		end if
	end if
!------------------------------read the quantabigL(4*subM,2)---
! every process need to read quantabigL

		reclength=4*subM*2
		! quantabigL is  integer(kind=4) so without *2
		inquire(file="quantabigL.tmp",exist=alive)
		if(alive) then
			open(unit=107,file="quantabigL.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			write(*,*) "quantabigL.tmp doesn't exist!"
			stop
		end if
		read(107,rec=nleft+1) quantabigL
		close(107)

return
end Subroutine enviro_bigL
