	Subroutine construct_exactelement(up,down,occupation)
	!    up down is the creation operator element
	USE Variables
	USE MPI

	implicit none

	integer :: exactsite,i
	integer :: i1,i2,i3,i4
	integer :: up(4,4,4,4,4),down(4,4,4,4,4),occupation(4,4,4,4,4)
	integer :: up1,down1,occupation1

	exactsite=2

	do i=1,2,1
	  if(myid==orbid(i)) then
	     do i1=0,3,1
	     do i2=0,3,1
	     do i3=0,3,1
	     do i4=0,3,1
		call exactelement(i,i1,i2,i3,i4,up1,down1,occupation1)
		up(i1+1,i2+1,i3+1,i4+1,i)=up1
		down(i1+1,i2+1,i3+1,i4+1,i)=down1
		occupation(i1+1,i2+1,i3+1,i4+1,i)=occupation1
	     end do
	     end do
	     end do
	     end do

	  else if(myid==orbid(norbs-2+i)) then
	     do i1=0,3,1
	     do i2=0,3,1
	     do i3=0,3,1
	     do i4=0,3,1
		call exactelement(i,i1,i2,i3,i4,up1,down1,occupation1)
		up(i1+1,i2+1,i3+1,i4+1,i+2)=up1
		down(i1+1,i2+1,i3+1,i4+1,i+2)=down1
		occupation(i1+1,i2+1,i3+1,i4+1,i+2)=occupation1
	     end do
	     end do
	     end do
	     end do
	 end if



	return
	end Subroutine construct_exactelement
	

!------------------------------------------------------------
	Subroutine exactelement(i,i1,i2,i3,i4,up,down,occupation)
	USE Variables
	USE MPI

	implicit none

	integer :: i,i1,i2,i3,i4
	integer :: up,down,occpuation

	if(i==1) then
		if(i2/=i4) then
			up=0
			down=0
			occupation=0
		else if(i1=1 .and. i3=0) then
			up=1
			down=0
			occupation=0
		else if(i1=2 .and. i3=0) then
			up=0
			down=1
			occupation=0
		else if(i1=1 .and. i3=1) then
			up=0
			down=0
			occupation=1
		else if(i1=3 .and. i3=1) then
			up=0
			down=-1
			occupation=0
		else if(i1=2 .and. i3=2 ) then
			up=0
			down=0
			occupation=1
		else if(i1=3 .and. i3=2) then
			up=1
			down=0
			occupation=0
		else if(i1=3 .and. i3=3) then
			up=0
			down=0
			occupation=2
		else
			up=0
			down=0
			occupation=0
		end if
	else if(i==2) then
		if(i1/=i3) then
			up=0
			down=0
			occupation=0
		else if(i2=1 .and. i4=0) then
			up=1*(-1)**i1
			down=0
			occupation=0
		else if(i2=2 .and. i4=0) then
			up=0
			down=1*(-1)**i1
			occupation=0
		else if(i2=1 .and. i4=1) then
			up=0
			down=0
			occupation=1
		else if(i2=3 .and. i4=1) then
			up=0
			down=-1*(-1)**i1
			occupation=0
		else if(i2=2 .and. i4=2 ) then
			up=0
			down=0
			occupation=1
		else if(i2=3 .and. i4=2) then
			up=1*(-1)**i1
			down=0
			occupation=0
		else if(i2=3 .and. i4=3) then
			up=0
			down=0
			occupation=2
		else
			up=0
			down=0
			occupation=0
		end if
	end if

	return

	end Subroutine exactelement

!--------------------------------------------------------------
	Subroutine Construct_exactH(up,down,occpuation)

	USE variables
	USE mpi
	implicit none














