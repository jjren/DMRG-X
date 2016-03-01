Program main
! contruct the coordinate of the polyacene for Swapan Pati's dmrg code


implicit none

real,parameter :: bonda=1.397D0,bondb=1.397D0,trans=-2.4D0
!bonda is the edge bond length, bondb is the laddar bond length
!trans is the transfer integral
integer :: block,nsites,sweep
real(kind=8),allocatable:: coord(:,:),nuclQ(:)
integer :: error,i,a,b,c,d
character(len=10) :: file1
character(len=1) :: symbol

open(unit=10,file="dst.inp",status="replace")
open(unit=12,file="eh.inp",status="replace")
open(unit=14,file="mol.inp",status="replace")
open(unit=18,file="dist.4",status="replace")
open(unit=20,file="maxsite.inp",status="replace")

write(*,*) "how many block"
read(*,*) block

allocate(coord(3,4*block+2),stat=error)
if(error/=0) stop
allocate(nuclQ(4*block+2),stat=error)
if(error/=0) stop
nuclQ=1.0D0
symbol='C'

!--------------------------------
!-----dst.inp--------------------
do i=1,4*block+2,1
	coord(3,i)=0.0D0
	if(Mod(i,4)==1) then
		coord(1,i)=2.0D0*bonda*sqrt(3.0D0)/2.0D0*dble((i-1)/4)
		coord(2,i)=0.0D0
	else if(Mod(i,4)==2) then
		coord(1,i)=2.0D0*bonda*sqrt(3.0D0)/2.0D0*dble((i-2)/4)
		coord(2,i)=-bondb
	else if(Mod(i,4)==3) then
		coord(1,i)=bonda*sqrt(3.0D0)/2.0D0*dble((i-1)/2)
		coord(2,i)=0.5D0*bonda
	else
		coord(1,i)=bonda*sqrt(3.0D0)/2.0D0*dble(i/2-1)
		coord(2,i)=-bondb-0.5D0*bonda
	end if

	if(Mod(i,4)==1 .or. Mod(i,4)==0) then
		write(12,*) i,'1'
	else
		write(12,*) i,'0'
	end if
	
	write(10,'(1A,4F16.12)') symbol,coord(1:3,i),nuclQ(i)
end do
!-------------dst.inp-----------------


!---------------dist.4------------------
write(18,*) coord(1:3,1)
write(18,*) coord(1:3,2)
write(18,*) coord(1:3,4*block+1)
write(18,*) coord(1:3,4*block+2)
!---------------dist.4-------------
!---------------------------------


!-------------------------------
!---------"mol.inp"
a=0
b=1
c=2
d=4


do i=3,2*block+1,1
	write(14,*) 2*i
	if(Mod(i,4)==3) then
		write(14,*) a,d,a
		write(14,*) i,i-2,i,1-i,-i,i-1,-i,2-i
		write(14,*) trans,trans,trans,trans
	else if(Mod(i,4)==0) then
		write(14,*) b,c,b
		write(14,*) i,i-2,-i,2-i
		write(14,*) trans,trans
		write(14,*) i-1,1-i
		write(14,*) trans
	else if(Mod(i,4)==1) then
		write(14,*) b,d,a
		write(14,*) -i,i-1,-i,2-i,i,i-2,i,1-i
		write(14,*) trans,trans,trans,trans
	else
		write(14,*) b,d,b
		write(14,*) i,i-1,i,i-2,-i,1-i,-i,2-i
		write(14,*) trans,trans,trans,trans
		write(14,*) i-1,1-i
		write(14,*) trans
	end if
end do
!---mol.inp----------------------------------
!-----------------------------------------------



!-------------------------------------------
!---------mol.nsite-------------------------
! there are two special cases because of the boundry

write(*,*) "input filename mol.nsite"
read(*,*) file1
open(unit=16,file=file1,status="replace")

nsites=4*block+2

do i=2,4*block,1
	write(16,*) i
	if(i==2) then
		write(16,'(3I2,F8.2)') a,d,a,trans
		write(16,'(10I6)') i,i-1,i,-(nsites-i-1),-(nsites-i),i-1,-(nsites-i),-(nsites-i)+2
		write(16,*) trans,trans,trans,trans
	else if(i==4*block) then
		write(16,'(3I2,F8.2)') a,d,a,trans
		write(16,'(10I6)') i,i-2,i,-(nsites-i-1),-(nsites-i),-(nsites-i-1),-(nsites-i),i-1
		write(16,*) trans,trans,trans,trans
		else if (Mod(i,4)==2 ) then
		write(16,'(3I2,F8.2)') a,d+1,a,trans
		write(16,'(10I6)') i,i-1,i,i-2,i,-(nsites-i-1),-(nsites-i),i-1,-(nsites-i),-(nsites-i)+2
		write(16,*) trans,trans,trans,trans,trans
	else if(Mod(i,4)==3) then
		write(16,'(3I2,F8.2)') a,d,a,trans
		write(16,'(10I6)') i,i-2,i,-(nsites-i-1),-(nsites-i),i-1,-(nsites-i),-(nsites-i-2)
		write(16,*) trans,trans,trans,trans
	else if(Mod(i,4)==0) then
		write(16,'(3I2,F8.2)') a,d+1,a,trans
		write(16,'(10I6)') i,i-2,i,-(nsites-i-1),-(nsites-i),-(nsites-i-1),-(nsites-i),-(nsites-i-2),-(nsites-i),i-1
		write(16,*) trans,trans,trans,trans,trans
	else
		write(16,'(3I2,F8.2)') b,d,a,trans
		write(16,'(10I6)') i,i-2,i,-(nsites-i-1),-(nsites-i),i-1,-(nsites-i),-(nsites-i-2)
		write(16,*) trans,trans,trans,trans
	end if
end do
!--------------mol.nsites--------------------------------------
!---------------------------------------------------------------

!----------------maxsite.inp-------------------------------

write(*,*) "how many sweep"
read(*,*) sweep

write(20,*) 4*block+2,sweep
!-----------------------------------------------

deallocate(coord)
close(10)
close(12)
close(14)
close(16)
close(18)
close(20)
end program
