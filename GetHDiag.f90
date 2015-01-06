Subroutine GetHDiag(HDIAG)
! This subroutine is to get the diagonal element of the Hamiltonian
! from all the process which can be used in davidson diagonalization

	use mpi
	use variables

	implicit none

	real(kind=8) :: HDIAG(ngoodstates)
	real(kind=8),allocatable :: buffermat(:,:),buffermat0(:,:,:),Hdiagdummy(:)
	integer :: operaindex
	integer :: status(MPI_STATUS_SIZE)
	integer :: i,error,j,k,m


	if(myid==0) then
		write(*,*) "enter in GetHDiag subroutine"
	end if

	if(myid/=0) then
		allocate(buffermat(4*subM,3),stat=error)
		if(error/=0) stop
		buffermat=0.0D0
	else
		allocate(buffermat0(4*subM,3,norbs),stat=error)
		if(error/=0) stop
		buffermat0=0.0D0
		allocate(Hdiagdummy(16*Lrealdim*Rrealdim),stat=error)
		if(error/=0) stop
		Hdiagdummy=0.0D0
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
				buffermat(j,:)=operamatbig(j,j,(operaindex-1)*3+1:operaindex*3)
			end do
! send the diagonal element to the 0 process
			call MPI_SEND(buffermat,12*subM,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(:,:,i),12*subM,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
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
				buffermat(k,:)=operamatbig(k,k,(operaindex-1)*3+1:operaindex*3)
			end do
! send the diagonal element to the 0 process
			call MPI_SEND(buffermat,12*subM,mpi_real8,0,j,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(:,:,j),12*subM,mpi_real8,orbid(j),j,MPI_COMM_WORLD,status,ierr)
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
! transfer integral contribution
		do i=1,nleft+1,1
		do j=norbs,norbs-nright,-1
			if(bondlink(i,j)==1) then
! add the phase al^+*ar
				do m=1,4*Lrealdim,1
				buffermat0(m,1:2,i)=buffermat0(m,1:2,i)*((-1.0D0)**(mod(quantabigL(m,1),2)))
				end do
! m represents up and down
				do m=1,2,1
				do k=1,4*Rrealdim,1
					Hdiagdummy((k-1)*4*Lrealdim+1:k*4*Lrealdim)=buffermat0(1:4*Lrealdim,m,i)*buffermat0(k,m,j)*t(i,j)*2.0D0&
											+Hdiagdummy((k-1)*4*Lrealdim+1:k*4*Lrealdim)
				end do
				end do
			end if
		end do
		end do
! PPP term contribution
		if(logic_PPP==1) then
			do i=1,nleft+1,1
			do j=norbs,norbs-nright,-1
				do k=1,4*Rrealdim,1
					Hdiagdummy((k-1)*4*Lrealdim+1:k*4*Lrealdim)=buffermat0(1:4*Lrealdim,3,i)*buffermat0(k,3,j)*pppV(i,j)&
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
			if((quantabigL(j,1)+quantabigR(i,1)==nelecs+ncharges) .and. &
				quantabigL(j,2)+quantabigR(i,2)==totalSz) then
				HDIAG(m)=Hdiagdummy((i-1)*4*Lrealdim+j)
				m=m+1
			end if
		end do
		end do
		
		m=m-1
		if(m/=ngoodstates) then
			write(*,*) "----------------------------------------------"
			write(*,*) "HDIAG number wrong! failed!,m=",m
			write(*,*) "----------------------------------------------"
			stop
		end if

	
	end if



	if(myid==0) then
		deallocate(buffermat0)
		deallocate(Hdiagdummy)
	else
		deallocate(buffermat)
	end if


return

end subroutine
