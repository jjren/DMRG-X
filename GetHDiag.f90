Subroutine GetHDiag(HDIAGnosymm)
! This subroutine is to get the diagonal element of the Hamiltonian(no symmetry)
! from all the process which can be used in davidson diagonalization
! the hopping term did not contribute anyting to the diagnal term 
! because the number of electrons is not equal
! so the diagnol term only need the PPP term operator

	use mpi
	use variables
	use communicate

	implicit none
	
	real(kind=8) :: HDIAGnosymm(ngoodstates)
	real(kind=8),allocatable :: buffermat(:),buffermat0(:,:),Hdiagdummy(:)
	integer :: operaindex
	integer :: status(MPI_STATUS_SIZE),ierr
	integer :: i,error,j,k,m

	call master_print_message("enter in GetHDiag subroutine")

	if(myid/=0) then
		allocate(buffermat(4*subM),stat=error)
		if(error/=0) stop
		buffermat=0.0D0
	else
		allocate(buffermat0(4*subM,norbs),stat=error)
		if(error/=0) stop
		buffermat0=0.0D0
		allocate(Hdiagdummy(16*Lrealdim*Rrealdim),stat=error)
		if(error/=0) stop
	end if

	! L space
	do i=1,nleft+1,1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
! copy the diagonal element to the buffermat
			do j=1,4*Lrealdim,1
				buffermat(j)=operamatbig(j,j,operaindex*3)
			end do
! send the diagonal element to the 0 process
			call MPI_SEND(buffermat,4*subM,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(1,i),4*subM,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
		end if
	end do

	! R space
	do j=norbs,norbs-nright,-1
		if(myid==orbid(j)) then
			if(mod(j,nprocs-1)==0) then
				operaindex=j/(nprocs-1)
			else
				operaindex=j/(nprocs-1)+1
			end if
! copy the diagonal element to the buffermat
			do k=1,4*Rrealdim,1
				buffermat(k)=operamatbig(k,k,operaindex*3)
			end do
! send the diagonal element to the 0 process
			call MPI_SEND(buffermat,4*subM,mpi_real8,0,j,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(1,j),4*subM,mpi_real8,orbid(j),j,MPI_COMM_WORLD,status,ierr)
		end if
	end do

	if(myid==0) then
		Hdiagdummy=0.0D0
! HL contribution
		do i=1,4*Rrealdim,1
			do j=1,4*Lrealdim,1
				Hdiagdummy((i-1)*4*Lrealdim+j)=Hbig(j,j,1)
			end do
		end do
! HR contribution
		do i=1,4*Rrealdim,1
			Hdiagdummy((i-1)*4*Lrealdim+1:i*4*Lrealdim)=Hbig(i,i,2)+Hdiagdummy((i-1)*4*Lrealdim+1:i*4*Lrealdim)
		end do
! transfer integral contribution the contribute is zero
! PPP term contribution
		if(logic_PPP==1) then
			do i=1,nleft+1,1
			do j=norbs,norbs-nright,-1
				do k=1,4*Rrealdim,1
					Hdiagdummy((k-1)*4*Lrealdim+1:k*4*Lrealdim)=buffermat0(1:4*Lrealdim,i)*buffermat0(k,j)*pppV(i,j)&
											+Hdiagdummy((k-1)*4*Lrealdim+1:k*4*Lrealdim)
				end do
			end do
			end do
		end if
! copy Hdiagdummy to Hdiag
! and ignore those diag element corresponding states without good quantum number
		m=1
		do i=1,4*Rrealdim,1
		do j=1,4*Lrealdim,1
			if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
				quantabigL(j,2)+quantabigR(i,2)==totalSz) then
				HDIAGnosymm(m)=Hdiagdummy((i-1)*4*Lrealdim+j)
				m=m+1
			end if
		end do
		end do
		
		m=m-1
		if(m/=ngoodstates) then
			call master_print_message(m,"HDIAGnosymm number wrong! failed!,m=")
		end if
	end if

	if(myid==0) then
		deallocate(buffermat0)
		deallocate(Hdiagdummy)
	else
		deallocate(buffermat)
	end if

return

end subroutine GetHDiag
