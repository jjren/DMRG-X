Subroutine system_bigR
! construct the R+sigmaR subspace operator matrix in 4M basis

	use variables
	use mpi
	use mathlib
	use communicate

	implicit none
	
	integer :: operaindex,error,i,j,k,l,ierr
	real(kind=r8),allocatable :: Hbuffer(:,:),operabuffer(:,:,:)
	integer :: status(MPI_STATUS_SIZE)
	real(kind=r8) :: II(4,4),Im(subM,subM)
	integer(kind=i4),allocatable :: phase(:,:)
	

	if(myid==0) then
		write(*,*) "enter in subroutine system_bigR"
	end if
	
	allocate(phase(4*subM,4*subM),stat=error)
	if(error/=0) stop

	! construct the unit matrix
	II=0.0D0
	do i=1,4,1
		II(i,i)=1.0D0
	end do
	Im=0.0D0
	do i=1,subM,1
		Im(i,i)=1.0D0
	end do

	
! construct the R subspace operator matrix in 4M basis
	do i=norbs,norbs-nright+1,-1
	if(myid==orbid(i)) then
		if(mod(i,nprocs-1)==0) then
			operaindex=i/(nprocs-1)
		else
			operaindex=i/(nprocs-1)+1
		end if
		
		do j=1,3,1
		if(j<=2) then
			phase(:,1:4*Rrealdim:4)=1
			phase(:,2:4*Rrealdim:4)=-1
			phase(:,3:4*Rrealdim:4)=-1
			phase(:,4:4*Rrealdim:4)=1
			call directproduct(II,4,operamatsma(1:Rrealdim,1:Rrealdim,3*(operaindex-1)+j),Rrealdim,operamatbig(1:Rrealdim*4,1:Rrealdim*4,3*(operaindex-1)+j),phase(1:4*Rrealdim,1:4*Rrealdim))
		else
			call directproduct(II,4,operamatsma(1:Rrealdim,1:Rrealdim,3*(operaindex-1)+j),Rrealdim,operamatbig(1:Rrealdim*4,1:Rrealdim*4,3*(operaindex-1)+j))
		end if
		end do
		call MPI_SEND(operamatsma(:,:,3*(operaindex-1)+1:3*(operaindex-1)+3),3*subM*subM,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
	end if
	end do


! construct the sigmaR subspace operator matrix in 4M basis
	if(myid==orbid(norbs-nright)) then
		if(mod(norbs-nright,nprocs-1)==0) then
			operaindex=(norbs-nright)/(nprocs-1)
		else
			operaindex=(norbs-nright)/(nprocs-1)+1
		end if
		
		do i=1,3,1
		call directproduct(onesitemat(:,:,i),4,Im(1:Rrealdim,1:Rrealdim),Rrealdim,operamatbig(1:Rrealdim*4,1:Rrealdim*4,3*(operaindex-1)+i))
		end do
	end if

! cosntruct the R+sigmaR Hamiltonian operator in 4M basis
	if(myid==0) then
		allocate(Hbuffer(4*subM,4*subM),stat=error)
		if(error/=0) stop
		
	!	R space Hamiltonian contribute
		Hbig(:,:,2)=0.0D0
		call directproduct(II,4,Hsma(1:Rrealdim,1:Rrealdim,2),Rrealdim,Hbuffer(1:4*Rrealdim,1:4*Rrealdim))
		Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbuffer(1:4*Rrealdim,1:4*Rrealdim)
!-------------------------------------------------------
!     R sigmaL interaction operator contribute
		allocate(operabuffer(subM,subM,3),stat=error)
		if(error/=0) stop

		do i=norbs,norbs-nright+1,-1
		!	write(*,*) "myid",myid,"i will recv operamatsmaR",orbid(i)
			call MPI_RECV(operabuffer,3*subM*subM,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
		 !     write(*,*) "myid",myid,"i have recv operamatsmaR",orbid(i)
			
			!     transfer integral term
			if(bondlink(i,norbs-nright)==1) then
			do j=1,2,1
				Hbuffer=0.0D0
				phase(:,1:4*Rrealdim:4)=1
				phase(:,2:4*Rrealdim:4)=-1
				phase(:,3:4*Rrealdim:4)=-1
				phase(:,4:4*Rrealdim:4)=1
				call directproduct(onesitemat(:,:,j+3),4,operabuffer(1:Rrealdim,1:Rrealdim,j),Rrealdim,Hbuffer(1:4*Rrealdim,1:4*Rrealdim),phase(1:4*Rrealdim,1:4*Rrealdim))
				Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbig(1:4*Rrealdim,1:4*Rrealdim,2)+(Hbuffer(1:4*Rrealdim,1:4*Rrealdim)+transpose(Hbuffer(1:4*Rrealdim,1:4*Rrealdim)))*(-1.0D0)*t(i,norbs-nright)
				! here the -1.0D0 transfer from ai*aj^+ to aj^+*ai
			end do
			end if
				Hbuffer=0.0D0
			!     ppp term
				call directproduct(onesitemat(:,:,3),4,operabuffer(1:Rrealdim,1:Rrealdim,3),Rrealdim,Hbuffer(1:4*Rrealdim,1:4*Rrealdim))
				Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbig(1:4*Rrealdim,1:4*Rrealdim,2)+Hbuffer(1:4*Rrealdim,1:4*Rrealdim)*pppV(i,norbs-nright)
		end do
!--------------------------------------------------------------
!     sigmaL Hamiltonian contribute. site energy+HubbardU
		Hbuffer=0.0D0
		do i=1,4*Rrealdim,4
			Hbuffer(i+1,i+1)=t(norbs-nright,norbs-nright)
			Hbuffer(i+2,i+2)=t(norbs-nright,norbs-nright)
			Hbuffer(i+3,i+3)=2.0D0*t(norbs-nright,norbs-nright)+HubbardU(norbs-nright)
		end do
		Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbuffer(1:4*Rrealdim,1:4*Rrealdim)+Hbig(1:4*Rrealdim,1:4*Rrealdim,2)
!-------------------------------------------------------------------
!	do i=1,16,1
!		write(*,'(16F5.1)') Hbig(i,1:Rrealdim*4,2)
!	end do
	
		if(logic_spinreversal/=0) then
			do i=1,Rrealdim,1
				symmlinkbig((i-1)*4+1,1,2)=((abs(symmlinksma(i,1,2))-1)*4+1)*&
				sign(1,symmlinksma(i,1,2))
				symmlinkbig((i-1)*4+2,1,2)=((abs(symmlinksma(i,1,2))-1)*4+3)*&
				sign(1,symmlinksma(i,1,2))
				symmlinkbig((i-1)*4+3,1,2)=((abs(symmlinksma(i,1,2))-1)*4+2)*&
				sign(1,symmlinksma(i,1,2))
				symmlinkbig((i-1)*4+4,1,2)=((abs(symmlinksma(i,1,2))-1)*4+4)*&
				sign(1,symmlinksma(i,1,2))*(-1)
			end do
		end if
	deallocate(Hbuffer)
	deallocate(operabuffer)
	end if

	deallocate(phase)
	
	return
end subroutine system_bigR

