program Tool_MO_Display
	implicit none

	real(kind=8),parameter :: decay=0.1D0,AutoA=0.52917721092D0
	real(kind=8) :: origin(3),vector(3,3)
	integer :: ngrid(3)
	real(kind=8),allocatable :: coord(:,:),coeff(:,:),grid(:,:,:)
	real(kind=8) :: x,distance
	integer :: natoms,orbbegin,orbend
	integer :: ix,iy,iz,i,j,k,l,istate
	character(len=1),allocatable :: atomc(:)
	integer,allocatable :: atomicn(:)
	character(len=20) :: filename

	open(unit=10,file="coord.xyz",status="old")
	read(10,*) natoms
	read(10,*)

	allocate(coord(3,natoms))
	allocate(atomc(natoms))
	allocate(atomicn(natoms))
	allocate(coeff(natoms,natoms))

	do i=1,natoms,1
		read(10,*) atomc(i),coord(1:3,i)
		if(atomc(i)=='C') then
			atomicn(i)=6
		else if(atomc(i)=='S') then
			atomicn(i)=16
		end if
	end do
	close(10)

	open(unit=11,file="cubegen.inp",status="old")
	read(11,*) origin
	do i=1,3,1
		read(11,*) vector(:,i)
	end do
	read(11,*) ngrid
	read(11,*) orbbegin,orbend
	close(11)

	open(unit=12,file="MO.out",status="old")
	do i=1,natoms,1
		read(12,*) 
		read(12,*) coeff(1:natoms,i)
	end do
	
	allocate(grid(0:ngrid(1)-1,0:ngrid(2)-1,0:ngrid(3)-1))

	do istate=orbbegin,orbend,1
		grid=0.0D0
		do ix=0,ngrid(1)-1,1
		do iy=0,ngrid(2)-1,1
		do iz=0,ngrid(3)-1,1
			x=0.0D0
			do i=1,natoms,1
				distance=0.0D0
				do j=1,3,1
					distance=distance+(origin(j)+ix*vector(j,1)+iy*vector(j,2)+iz*vector(j,3)-coord(j,i))**2
				end do
				x=x+exp(-1.0D0*sqrt(distance)/decay)*coeff(i,istate)
			end do
			grid(ix,iy,iz)=x
		!	write(*,* ) sqrt(distance),x
		end do
		end do
		end do

		write(filename,'(i5.5,a8)') istate,'.cube'
		open(unit=111,file=trim(filename),status="replace")
		write(111,*) "haha"
		write(111,*) "haha"
		write(111,'(I5,3f12.6)') natoms,origin(1:3)/AutoA
		do j=1,3,1
			write(111,'(I5,3f12.6)') ngrid(j),vector(:,j)/AutoA
		end do
		do j=1,natoms,1
			write(111,'(I5,4f12.6)') atomicn(j),dble(atomicn(j)),coord(:,j)/AutoA
		end do
		do j=0,ngrid(1)-1,1
		do k=0,ngrid(2)-1,1
		do l=0,ngrid(3)-1,6
			write(111,'(6f12.6)') grid(j,k,l:l+5)
		end do
		end do
		end do
		close(111)
	end do
!	do i=orbbegin,orbend,1
!		write(filename,'(i5.5,a8)') i,'.cube'
!		open(unit=111,file=trim(filename),form="unformatted",status="replace")
!		write(111) "haha"
!		write(111) "haha"
!		write(111) natoms,origin
!		do j=1,3,1
!			write(111) ngrid(j),vector(:,j)
!		end do
!		do j=1,natoms,1
!			write(111) atomicn(j),coord(:,j)
!		end do
!		do j=0,ngrid(3)-1,1
!			write(111) grid(:,:,j)
!		end do
!		close(111)
!	end do

end program Tool_MO_Display
