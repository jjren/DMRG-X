subroutine C2_copy(direction)
! this subroutine is to copy the L space operator and R space operator in
! the C2 like symmetry condition
! direction is 'i' 'l' 'r'
	use mpi
	use variables

	implicit none
	character(len=1) :: direction
	integer :: i,j
	integer :: operaindex1,operaindex2
	integer :: status(MPI_STATUS_SIZE)
	
	if(myid==0) then 
		write(*,*) "enter in C2_copy subroutine"
	end if

	if(direction=='i') then
		do i=1,nleft+1,1
			if(myid==orbid(i)) then
				if(mod(i,nprocs-1)==0) then
					operaindex1=i/(nprocs-1)
				else
					operaindex1=i/(nprocs-1)+1
				end if
				if(myid==orbid(norbs-i+1)) then
					if(mod(norbs-i+1,nprocs-1)==0) then
						operaindex2=(norbs-i+1)/(nprocs-1)
					else
						operaindex2=(norbs-i+1)/(nprocs-1)+1
					end if
					operamatbig(:,:,3*(operaindex2-1)+1:3*operaindex2)=&
						operamatbig(:,:,3*(operaindex1-1)+1:3*operaindex1)
				else
					call MPI_SEND(operamatbig(:,:,3*(operaindex1-1)+1:3*operaindex1),3*16*subM*subM,mpi_real8,orbid(norbs-i+1),i,MPI_COMM_WORLD,ierr)
				end if
			else if(myid==orbid(norbs-i+1)) then
					if(mod(norbs-i+1,nprocs-1)==0) then
						operaindex2=(norbs-i+1)/(nprocs-1)
					else
						operaindex2=(norbs-i+1)/(nprocs-1)+1
					end if
					call MPI_RECV(operamatbig(:,:,3*(operaindex2-1)+1:3*operaindex2),3*16*subM*subM,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
			end if
		end do
	end if

return
end subroutine


