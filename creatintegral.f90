program creatintegral

implicit none


integer :: nrings,atomindex(6),nbonds,norbs,link1,link2
real(kind=8) :: t,hubbardU
integer,allocatable :: bondlink(:,:)
integer :: i,error,j

! creat the bondlink information
	open(unit=10,file="inp",status="replace")
! for example polyacene
	nrings=2
	do i=1,6,1
		atomindex(i)=i
	end do
	do i=1,nrings,1
		write(10,*) atomindex(1),atomindex(2)
		write(10,*) atomindex(1),atomindex(3)
		write(10,*) atomindex(2),atomindex(4)
		write(10,*) atomindex(3),atomindex(5)
		write(10,*) atomindex(4),atomindex(6)
		atomindex=atomindex+4
	end do
		atomindex=atomindex-4
		write(10,*) atomindex(5),atomindex(6)

	rewind(10)


	nbonds=5*nrings+1
	norbs=4*nrings+2

	allocate(bondlink(norbs,norbs),stat=error)
	if(error/=0) stop
	

	do i=1,nbonds,1
		read(10,*) link1,link2
		bondlink(link1,link2)=1
		bondlink(link2,link1)=1
	end do

	close(10)

	open(unit=11,file="integral.inp",status="replace")
! transfer integeral
	t=-2.4D0
	do i=1,norbs,1
	do j=1,norbs,1
		if(bondlink(j,i)==1) then
			write(11,*) t
		end if
	end do
	end do

! site energy
	t=0.0D0
	do i=1,norbs,1
		write(11,*) t
	end do

! hubbardU
	hubbardU=11.26D0
	do i=1,norbs,1
		write(11,*) hubbardU
	end do
	close(11)

end 
! 

