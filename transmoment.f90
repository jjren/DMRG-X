subroutine transmoment
! this subroutine is to calculate the transition moment between gs and ex
! only used in the state average method
	use mpi
	use variables
	use blas95
	use f95_precision

	implicit none

	integer :: operaindex
	integer :: i,j,k,error
	integer :: status(MPI_STATUS_SIZE)
	real(kind=8),allocatable :: midmat(:,:)
	real(kind=8),allocatable :: imoment(:),zeromoment(:,:)
	! imoment(1) is a mediate number
	! imoment(2:nstate) is <faiex|ni|faigs>
	real(kind=8) :: moment(4,nstate)


	if(myid==0) then
		write(*,*) "enter transmoment subroutine"
		if(nstate==1) then
			write(*,*) "---------------------------"
			write(*,*) "nstate==1 wrong!"
			write(*,*) "---------------------------"
			stop
		end if
	end if
	
	if(myid/=0) then
		allocate(coeffIF(4*subM,4*subM,nstate),stat=error)
		if(error/=0) stop
		allocate(midmat(4*subM,4*subM),stat=error)
		if(error/=0) stop
		allocate(imoment(nstate),stat=error)
		if(error/=0) stop
	else
		allocate(zeromoment(nstate,norbs),stat=error)
		if(error/=0) stop
		moment=0.0D0
	end if

	call MPI_bcast(coeffIF,16*subM*subM*nstate,MPI_real8,0,MPI_COMM_WORLD,ierr)

	do i=1,norbs,1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
		! it does not need to remove the nuclearQ
		! because the gs and ex state is orthogonal
		!	do j=1,4*subM,1
		!		operamatbig(j,j,3*operaindex)=operamatbig(j,j,3*operaindex)+nuclQ(i)
		!	end do
			if(i<=nleft+1) then
				call gemm(operamatbig(:,:,3*operaindex),coeffIF(:,:,1),&
				midmat,'N','N',1.0D0,0.0D0)
				imoment=0.0D0
				do j=2,nstate,1
					do k=1,4*subM,1
						imoment(1)=dot(coeffIf(k,:,j),midmat(k,:))
						imoment(j)=imoment(j)+imoment(1)
					end do
				end do
			else
				call gemm(operamatbig(:,:,3*operaindex),coeffIF(:,:,1),&
				midmat,'N','T',1.0D0,0.0D0)
				imoment=0.0D0
				do j=2,nstate,1
					do k=1,4*subM,1
						imoment(1)=dot(coeffIf(k,:,j),midmat(:,k))
						imoment(j)=imoment(j)+imoment(1)
					end do
				end do
			end if

			call MPI_SEND(imoment,nstate,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(zeromoment(:,i),nstate,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
			do j=2,nstate,1
				do k=1,3,1
					moment(k,j)=moment(k,j)+zeromoment(j,i)*coord(k,i)
				end do
			end do
		end if
	end do
	
	if(myid==0) then
		write(*,*) "transition moment"
		do i=2,nstate,1
			moment(4,i)=sqrt(moment(1,i)**2+moment(2,i)**2+moment(3,i)**2)
			write(*,'(4F10.5)') moment(:,i)
		end do
	end if


end subroutine
