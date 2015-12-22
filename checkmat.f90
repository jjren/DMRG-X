module checkmat_mod
	use variables
	use communicate
	implicit none

contains
subroutine checkmat(Hmat,Hcol,Hrow,mat1,mat1col,mat1row,quantaL,quantaR,csize,domain)
! this subroutine is used to output the important matrix 
	implicit none
	
	character(len=1),intent(in) :: domain
	character(len=3),intent(in) :: csize
	integer(kind=i4),intent(in) :: Hcol(:,:),Hrow(:,:),&
	mat1col(:,:),mat1row(:,:),quantaL(:,:),quantaR(:,:)
	real(kind=r8),intent(in) :: Hmat(:,:),mat1(:,:)


	! local
	character(len=50) filename
	integer :: orbstart,orbend,dim1,Hindex,operaindex,units
	logical :: alive
	integer :: i,j

	call master_print_message("enter subroutine checkmat")

	if(domain=="L") then
		orbstart=1
		Hindex=1
		if(csize=="big") then
			orbend=nleft+1
			dim1=4*Lrealdim
		else if(csize=="sma") then
			orbend=nleft
			dim1=Lrealdim
		else
			stop
		end if
	else if(domain=='R') then
		orbend=norbs
		Hindex=2
		if(csize=="big") then
			orbstart=norbs-nright
			dim1=4*Rrealdim
		else if(csize=="sma") then
			orbstart=norbs-nright+1
			dim1=Rrealdim
		else
			stop
		end if
	else 
		stop
	end if
	
	units=myid+100

	write(filename,'(i5,a12)') myid,'checkmat.tmp'
	inquire(file=trim(filename),exist=alive)
	if(alive) then
		open(unit=units,file=trim(filename),status="old",position="append")
	else
		open(unit=units,file=trim(filename),status="replace")
	end if

	write(units,*) "isweep=",isweep
	do i=orbstart,orbend,1
		if(myid==orbid1(i,1)) then
			write(units,*) i
			do j=1,3,1
				operaindex=orbid1(i,2)*3-3+j
				write(units,*) mat1row(1:dim1+1,operaindex)
			end do
			write(units,*) 
		end if
	end do

	if(myid==0) then
		write(units,*) nleft
		write(units,*) Hrow(1:dim1+1,Hindex)
		write(units,*)
	end if
	
	if(domain=="L") then
		write(units,*) quantaL(1:dim1,1:2)
	else
		write(units,*) quantaR(1:dim1,1:2)
	end if

	close(units)
return

end subroutine checkmat


end module checkmat_mod

