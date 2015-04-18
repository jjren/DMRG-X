Subroutine system_bigRreverse
! construct the R+sigmaR subspace operator matrix in 4M basis
! the difference between bigRreverse and bigR is the order 
! in bigR the order is sigmaR,sigmaR+1....norbs
! in bigRreverse the order is norbs,norbs-1,...sigmaR
! bigRreverse is easy to include the C2 symmetry

	
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
		write(*,*) "enter in subroutine system_bigRreverse"
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
		call directproduct(operamatsma(1:Rrealdim,1:Rrealdim,3*(operaindex-1)+j),Rrealdim,II,4,operamatbig(1:Rrealdim*4,1:Rrealdim*4,3*(operaindex-1)+j))
		end do
		!write(*,*) "myid=",myid,"i will send operamatsma"
		call MPI_SEND(operamatsma(:,:,3*(operaindex-1)+1:3*(operaindex-1)+3),3*subM*subM,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		!write(*,*) "myid=",myid,"i have send operamatsma"
	end if
	end do


! construct the sigmaR subspace operator matrix in 4M basis
	if(myid==orbid(norbs-nright)) then
		!write(*,*) "i am",myid
		if(mod(norbs-nright,nprocs-1)==0) then
			operaindex=(norbs-nright)/(nprocs-1)
		else
			operaindex=(norbs-nright)/(nprocs-1)+1
		end if
		
		do i=1,3,1
		if(i<=2) then
			do j=1,4*Rrealdim,1
				if(mod(j,Rrealdim)==0) then
					k=Rrealdim
				else
					k=mod(j,Rrealdim)
				end if
				phase(:,j)=(-1)**(mod(quantasmaR(k,1),2))
			end do
		call directproduct(Im(1:Rrealdim,1:Rrealdim),Rrealdim,onesitemat(:,:,i),4,operamatbig(1:Rrealdim*4,1:Rrealdim*4,3*(operaindex-1)+i),phase(1:Rrealdim*4,1:Rrealdim*4))
		else
		call directproduct(Im(1:Rrealdim,1:Rrealdim),Rrealdim,onesitemat(:,:,i),4,operamatbig(1:Rrealdim*4,1:Rrealdim*4,3*(operaindex-1)+i))
		end if
		end do
	end if

! cosntruct the R+sigmaR Hamiltonian operator in 4M basis
	if(myid==0) then
		allocate(Hbuffer(4*subM,4*subM),stat=error)
		if(error/=0) stop
		
!     R Hamiltonian contribute
		Hbig(:,:,2)=0.0D0
		call directproduct(Hsma(1:Rrealdim,1:Rrealdim,2),Rrealdim,II,4,Hbuffer(1:4*Rrealdim,1:4*Rrealdim))
		Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbuffer(1:4*Rrealdim,1:4*Rrealdim)
!-------------------------------------------------------
!     R sigmaR interaction operator contribute
		allocate(operabuffer(subM,subM,3),stat=error)
		if(error/=0) stop

		do i=norbs,norbs-nright+1,-1
			!write(*,*) "myid",myid,"i will recv operamatsma",orbid(i)
			call MPI_RECV(operabuffer,3*subM*subM,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
		      !write(*,*) "myid",myid,"i have recv operamatsma",orbid(i)
			
			!     transfer integral term
			if(bondlink(i,norbs-nright)==1) then
			do j=1,2,1
				Hbuffer=0.0D0
				do k=1,4*Rrealdim,1
					if(mod(k,Rrealdim)==0) then
						l=Rrealdim
					else
						l=mod(k,Rrealdim)
					end if
					phase(:,k)=(-1)**(mod(quantasmaR(l,1),2))
				end do
				call directproduct(operabuffer(1:Rrealdim,1:Rrealdim,j),Rrealdim,onesitemat(:,:,3+j),&
				4,Hbuffer(1:4*Rrealdim,1:4*Rrealdim),phase(1:4*Rrealdim,1:4*Rrealdim))
				Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbig(1:4*Rrealdim,1:4*Rrealdim,2)+(Hbuffer(1:4*Rrealdim,1:4*Rrealdim)+transpose(Hbuffer(1:4*Rrealdim,1:4*Rrealdim)))*t(i,norbs-nright)
			end do
			end if
			!     ppp term
				Hbuffer=0.0D0
				call directproduct(operabuffer(1:Rrealdim,1:Rrealdim,3),Rrealdim,onesitemat(:,:,3),&
				4,Hbuffer(1:4*Rrealdim,1:4*Rrealdim))
				Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbig(1:4*Rrealdim,1:4*Rrealdim,2)+Hbuffer(1:4*Rrealdim,1:4*Rrealdim)*pppV(i,norbs-nright)
		end do
!--------------------------------------------------------------
!     sigmaR Hamiltonian contribute. site energy+HubbardU
		Hbuffer=0.0D0
		do i=1,Rrealdim,1
			Hbuffer(1*Rrealdim+i,1*Rrealdim+i)=t(norbs-nright,norbs-nright)
			Hbuffer(2*Rrealdim+i,2*Rrealdim+i)=t(norbs-nright,norbs-nright)
			Hbuffer(3*Rrealdim+i,3*Rrealdim+i)=t(norbs-nright,norbs-nright)*2.0D0+hubbardU(norbs-nright)
		end do
		Hbig(1:4*Rrealdim,1:4*Rrealdim,2)=Hbuffer(1:4*Rrealdim,1:4*Rrealdim)+Hbig(1:4*Rrealdim,1:4*Rrealdim,2)
!-------------------------------------------------------------------
!	do i=1,16,1
!		write(*,'(16F5.1)') Hbig(i,1:Rrealdim*4,1)
!	end do
	
		if(logic_spinreversal/=0) then
			symmlinkbig(1:Rrealdim,1,2)=symmlinksma(1:Rrealdim,1,2)
			symmlinkbig(Rrealdim+1:2*Rrealdim,1,2)=(abs(symmlinksma(1:Rrealdim,1,1))+2*&
				Rrealdim)*sign(1,symmlinksma(1:Rrealdim,1,2))
			symmlinkbig(2*Rrealdim+1:3*Rrealdim,1,2)=(abs(symmlinksma(1:Rrealdim,1,1))+&
				Rrealdim)*sign(1,symmlinksma(1:Rrealdim,1,2))
			symmlinkbig(3*Rrealdim+1:4*Rrealdim,1,2)=(abs(symmlinksma(1:Rrealdim,1,1))+3*&
				Rrealdim)*sign(1,symmlinksma(1:Rrealdim,1,2))*(-1)
		end if

	deallocate(Hbuffer)
	deallocate(operabuffer)
	end if
	
	deallocate(phase)
	return
end subroutine system_bigRreverse

