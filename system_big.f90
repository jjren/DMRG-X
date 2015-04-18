Subroutine System_Big(domain)

! construct the L+sigmaL/R+sigmaR subspace operator matrix in 4M basis
! two different conditions
! R space and logic_C2==0  ;  L space or (R space and logic_C2/=0)
	
	use variables
	use mpi
	use mathlib
	use symmetry
	use communicate
	use exit_mod

	implicit none
	
	character(len=1) :: domain   !L/R
	
	! local
	integer :: error,ierr
	integer :: i,j,k,j1,j2
	integer :: operaindex,orbstart,orbend,orbadd,Hindex,realdim
	! orbstart is from 1 or norbs-nright+1
	! orbend is nleft or norbs
	! orbadd is nleft+1 or norbs-nright
	! Hindex : HL/1 HR/2
	real(kind=r8) :: II(4,4)
	integer(kind=i4),allocatable :: phase(:)
	real(kind=r8),allocatable :: Hbuffer(:,:),operabuffer(:,:,:),recvoperabuffer(:,:,:)
	
	integer :: status(MPI_STATUS_SIZE),recvtag,sendrequest(norbs),recvrequest ! MPI flag

	call master_print_message("enter in subroutine system_big")

	if(domain=='L') then
		orbstart=1
		orbend=nleft
		orbadd=nleft+1
		realdim=Lrealdim
		Hindex=1
	else if(domain=='R') then
		orbstart=norbs-nright+1
		orbend=norbs
		orbadd=norbs-nright
		realdim=Rrealdim
		Hindex=2
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
	end if

!==================================================================
! set phase value
	allocate(phase(4*subM),stat=error)
	if(error/=0) stop
	if(domain=='R' .and. logic_C2==0) then
		phase(1:4*realdim:4)=1
		phase(2:4*realdim:4)=-1
		phase(3:4*realdim:4)=-1
		phase(4:4*realdim:4)=1
	else
		do j=1,4*realdim,1
			if(mod(j,realdim)==0) then
				k=realdim
			else
				k=mod(j,realdim)
			end if
			if(domain=='L') then
				phase(j)=(-1)**(mod(quantasmaL(k,1),2))
			else
				phase(j)=(-1)**(mod(quantasmaR(k,1),2))
			end if
		end do
	end if

! construct the unit matrix
	II=0.0D0
	do i=1,4,1
		II(i,i)=1.0D0
	end do
!=========================================================================

! construct the L/R(without sigmaL/R) subspace operator matrix in 4M basis
	do i=orbstart,orbend,1
	if(myid==orbid(i)) then
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if
		
		if(bondlink(i,orbadd)==1) then
		! block
		!	call MPI_SEND(operamatsma(1,1,3*operaindex-2),3*subM*subM,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		! non-block
			call MPI_ISEND(operamatsma(1,1,3*operaindex-2),3*subM*subM,mpi_real8,0,i,MPI_COMM_WORLD,sendrequest(i),ierr)
		else
		! block
		!	call MPI_SEND(operamatsma(1,1,3*operaindex),subM*subM,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		! non-block
			call MPI_ISEND(operamatsma(1,1,3*operaindex),subM*subM,mpi_real8,0,i,MPI_COMM_WORLD,sendrequest(i),ierr)
		end if
		
		if(domain=='R' .and. logic_C2==0 ) then
			do j=1,3,1
				if(j<=2) then
					call DirectProduct(II,4, &
						operamatsma(1:realdim,1:realdim,3*(operaindex-1)+j),realdim, &
						operamatbig(1:realdim*4,1:realdim*4,3*(operaindex-1)+j))
					do k=1,4*realdim
						operamatbig(:,k,3*(operaindex-1)+j)=operamatbig(:,k,3*(operaindex-1)+j)*DBLE(phase(k))
					end do
				else
					call DirectProduct(II,4, &
						operamatsma(1:realdim,1:realdim,3*(operaindex-1)+j),realdim, &
						operamatbig(1:realdim*4,1:realdim*4,3*(operaindex-1)+j))
				end if
			end do
		else
			do j=1,3,1
				call DirectProduct(operamatsma(1:realdim,1:realdim,3*(operaindex-1)+j),realdim, &
					II,4,operamatbig(1:realdim*4,1:realdim*4,3*(operaindex-1)+j))
			end do
		end if
	end if
	end do


! construct the sigmaL/R subspace operator matrix in 4M basis
	if(myid==orbid(orbadd)) then
		if(mod(orbadd,nprocs-1)==0) then
			operaindex=orbadd/(nprocs-1)
		else
			operaindex=orbadd/(nprocs-1)+1
		end if
		
		operamatbig(1:4*realdim,1:4*realdim,3*operaindex-2:3*operaindex)=0.0D0
		if(domain=='R' .and. logic_C2==0) then
			do i=1,3,1
				do j=1,4*realdim,4
					operamatbig(j:j+3,j:j+3,3*(operaindex-1)+i)=onesitemat(:,:,i)
				end do
			end do
		else 
			do i=1,3,1
				do j1=1,4,1 ! column
				do j2=1,4,1 ! row
					do k=1,realdim,1
						if(i<=2) then
							operamatbig(k+(j2-1)*realdim,k+(j1-1)*realdim,3*operaindex-3+i)=onesitemat(j2,j1,i)*DBLE(phase(k))
						else
							operamatbig(k+(j2-1)*realdim,k+(j1-1)*realdim,3*operaindex-3+i)=onesitemat(j2,j1,i)
						end if
					end do
				end do
				end do
			end do
		end if
	end if

! cosntruct the L+sigmaL/R+sigmaR Hamiltonian operator in 4M basis
	if(myid==0) then
		allocate(Hbuffer(4*subM,4*subM),stat=error)
		if(error/=0) stop
		
		Hbig(:,:,Hindex)=0.0D0
!===========================================================

		allocate(recvoperabuffer(subM,subM,3),stat=error)
		if(error/=0) stop
		allocate(operabuffer(subM,subM,3),stat=error)
		if(error/=0) stop
		
		call MPI_IRECV(recvoperabuffer(1,1,1),3*subM*subM,mpi_real8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,recvrequest,ierr)

!===========================================================

!     L/R Hamiltonian contribute
		if(domain=='R' .and. logic_C2==0) then
			call directproduct(II,4, &
				Hsma(1:realdim,1:realdim,Hindex),realdim, &
				Hbuffer(1:4*realdim,1:4*realdim))
		else
			call directproduct(Hsma(1:realdim,1:realdim,Hindex),realdim, &
				II,4,Hbuffer(1:4*realdim,1:4*realdim))
		end if
		Hbig(1:4*realdim,1:4*realdim,Hindex)=Hbuffer(1:4*realdim,1:4*realdim)

!===========================================================

!     sigmaL Hamiltonian contribute. site energy+HubbardU
		Hbuffer=0.0D0
		if(domain=='R' .and. logic_C2==0) then
			do i=1,4*realdim,4
				Hbuffer(i+1,i+1)=t(orbadd,orbadd)
				Hbuffer(i+2,i+2)=t(orbadd,orbadd)
				Hbuffer(i+3,i+3)=2.0D0*t(orbadd,orbadd)+HubbardU(orbadd)
			end do
		else
			do i=1,realdim,1
				Hbuffer(1*realdim+i,1*realdim+i)=t(orbadd,orbadd)
				Hbuffer(2*realdim+i,2*realdim+i)=t(orbadd,orbadd)
				Hbuffer(3*realdim+i,3*realdim+i)=2.0D0*t(orbadd,orbadd)+HubbardU(orbadd)
			end do
		end if
		Hbig(1:4*realdim,1:4*realdim,Hindex)=Hbuffer(1:4*realdim,1:4*realdim)+Hbig(1:4*realdim,1:4*realdim,Hindex)

!===========================================================

!     L sigmaL/R sigmaR interaction operator contribute

		do i=orbstart,orbend,1
		
		! block
		!	if(bondlink(orbadd,i)==1) then
		!		call MPI_RECV(operabuffer(1,1,1),3*subM*subM,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
		!	else
		!		call MPI_RECV(operabuffer(1,1,3),subM*subM,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
		!	end if

		! nonblock
			call MPI_WAIT(recvrequest,status,ierr)
			recvtag=status(MPI_TAG)
			if(bondlink(recvtag,orbadd)==1) then
				operabuffer=recvoperabuffer
			else
				operabuffer(:,:,3)=recvoperabuffer(:,:,1)
			end if

			if(i<orbend) then
				call MPI_IRECV(recvoperabuffer(1,1,1),3*subM*subM,mpi_real8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,recvrequest,ierr)
			end if
			
		!transfer integral term
			if(bondlink(recvtag,orbadd)==1) then
			do j=1,2,1
				if(domain=='R' .and. logic_C2==0) then
					call directproduct(onesitemat(:,:,j+3),4, &
						operabuffer(1:realdim,1:realdim,j),realdim, &
						Hbuffer(1:4*realdim,1:4*realdim))
					do k=1,4*realdim,1
						Hbuffer(:,k)=Hbuffer(:,k)*DBLE(phase(k))*(-1.0D0)
					end do
				else
					call directproduct(operabuffer(1:realdim,1:realdim,j),realdim, &
						onesitemat(:,:,3+j),4, &
						Hbuffer(1:4*realdim,1:4*realdim))
					do k=1,4*realdim,1
						Hbuffer(:,k)=Hbuffer(:,k)*DBLE(phase(k))
					end do
				end if
				Hbig(1:4*realdim,1:4*realdim,Hindex)=Hbig(1:4*realdim,1:4*realdim,Hindex)+&
						(Hbuffer(1:4*realdim,1:4*realdim)+transpose(Hbuffer(1:4*realdim,1:4*realdim)))*t(recvtag,orbadd)
			end do
			end if
		! ppp term
			if(domain=='R' .and. logic_C2==0) then
				call directproduct(onesitemat(:,:,3),4,operabuffer(1:realdim,1:realdim,3),realdim,Hbuffer(1:4*realdim,1:4*realdim))
			else
				call directproduct(operabuffer(1:realdim,1:realdim,3),realdim,onesitemat(:,:,3),4,Hbuffer(1:4*realdim,1:4*realdim))
			end if
			Hbig(1:4*realdim,1:4*realdim,Hindex)=Hbig(1:4*realdim,1:4*realdim,Hindex)+&
				Hbuffer(1:4*realdim,1:4*realdim)*pppV(recvtag,orbadd)
		end do
!=========================================================================
	
	! construct the symmmlinkbig
		if(logic_spinreversal/=0) then
			call Creatsymmlinkbig(realdim,domain,Hindex)
		end if

		deallocate(Hbuffer)
		deallocate(operabuffer)
		deallocate(recvoperabuffer)
	end if

! free the sendrequest
	do i=orbstart,orbend,1
		if(myid==orbid(i)) then
			call MPI_REQUEST_FREE(sendrequest(i),ierr)
		end if
	end do

	deallocate(phase)
	return
end subroutine System_Big

