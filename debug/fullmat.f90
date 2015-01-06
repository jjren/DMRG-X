subroutine fullmat
! this subroutine only used in test mode
! construct the fullmat in the 16M^2 basis

	use mpi
	use variables
	use mathlib
	USE LAPACK95
	USE F95_PRECISION

	implicit none
	integer :: i,j,k,l,m,operaindex,tmp
	integer :: error
	integer :: status(MPI_STATUS_SIZE)
	real(kind=8),allocatable :: buffermat(:,:,:),buffermat0(:,:,:,:),fullH(:,:),fullHdummy(:,:),I4M(:,:),goodH(:,:)
	real(kind=8),allocatable :: eigenvalue(:)!,z(:,:)
	!real(kind=8) :: vl,vu,abstol
	!integer,allocatable :: isuppz(:)
	integer :: info
	
	

	if(myid==0) then 
		write(*,*) "enter in fullmat subroutine"
	end if
	
	allocate(buffermat(4*subM,4*subM,3),stat=error)
	if(error/=0) stop
	buffermat=0.0D0
	if(myid==0) then
		allocate(buffermat0(4*subM,4*subM,3,norbs),stat=error)
		if(error/=0) stop
		allocate(fullH(16*Rrealdim*Lrealdim,16*Rrealdim*Lrealdim),stat=error)
		if(error/=0) stop
		allocate(fullHdummy(16*Rrealdim*Lrealdim,16*Rrealdim*Lrealdim),stat=error)
		if(error/=0) stop
		buffermat0=0.0D0
		fullH=0.0D0
		fullHdummy=0.0D0
! identity matrix
		allocate(I4M(4*subM,4*subM),stat=error)
		if(error/=0) stop
		I4M=0.0D0
		do i=1,4*subM,1
			I4M(i,i)=1.0D0
		end do
	end if

	do i=1,nleft+1,1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			buffermat=operamatbig(:,:,3*(operaindex-1)+1:3*operaindex)
			call MPI_SEND(buffermat,16*subM*subM*3,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(:,:,:,i),16*subM*subM*3,mpi_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
		end if
	end do

	do j=norbs,norbs-nright,-1
		if(myid==orbid(j)) then
			if(mod(j,nprocs-1)==0) then
				operaindex=j/(nprocs-1)
			else
				operaindex=j/(nprocs-1)+1
			end if
			buffermat=operamatbig(:,:,3*(operaindex-1)+1:3*operaindex)
			call MPI_SEND(buffermat,16*subM*subM*3,mpi_real8,0,j,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(:,:,:,j),16*subM*subM*3,mpi_real8,orbid(j),j,MPI_COMM_WORLD,status,ierr)
		end if
	end do
	


	if(myid==0) then
		do i=1,nleft+1,1
		do j=norbs,norbs-nright,-1
			if(bondlink(i,j)==1) then
				!do k=1,2,1
				!	buffermat(:,:,k)=transpose(buffermat0(:,:,k,i))
				!	call directproduct(buffermat(1:Lrealdim*4,1:Lrealdim*4,k),4*Lrealdim,&
				!	buffermat0(1:Rrealdim*4,1:Rrealdim*4,k,j),4*Rrealdim,fullHdummy)
				!	do l=1,4*Lrealdim,1
				!		fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Rrealdim)=&
				!		fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Rrealdim)*((-1.0D0)**(mod(quantabigL(l,1),2)))
				!	end do
				! transfer from al*ar^+ to ar^+*al
				!	fullH=(-1.0D0*fullHdummy+transpose(-1.0D0*fullHdummy))*t(i,j)+fullH
				
				do k=1,2,1
					buffermat(:,:,k)=transpose(buffermat0(:,:,k,j))
					call directproduct(buffermat0(1:Lrealdim*4,1:Lrealdim*4,k,i),4*Lrealdim,&
					buffermat(1:Rrealdim*4,1:Rrealdim*4,k),4*Rrealdim,fullHdummy)
					do l=1,4*Lrealdim,1
						fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Rrealdim)=&
						fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Rrealdim)*((-1.0D0)**(mod(quantabigL(l,1),2)))
					end do
					fullH=(fullHdummy+transpose(fullHdummy))*t(i,j)+fullH
				end do
				
			end if
			! PPP term
				call directproduct(buffermat0(1:Lrealdim*4,1:Lrealdim*4,3,i),4*Lrealdim,&
				buffermat0(1:Rrealdim*4,1:Rrealdim*4,3,j),4*Rrealdim,fullHdummy)
				fullH=fullH+fullHdummy*pppV(i,j)
		end do
		end do
! Hl contribute
		call directproduct(Hbig(1:Lrealdim*4,1:Lrealdim*4,1),4*Lrealdim,&
		I4M(1:Rrealdim*4,1:Rrealdim*4),4*Rrealdim,fullHdummy)
		fullH=fullH+fullHdummy
! HR contribute
		call directproduct(I4M(1:Lrealdim*4,1:Lrealdim*4),4*Lrealdim,&
		Hbig(1:Rrealdim*4,1:Rrealdim*4,2),4*Rrealdim,fullHdummy)
		fullH=fullH+fullHdummy

! write the full Hamiltonian out-----------------------------
		!write(*,*) "the full hamiltonian"
		!do i=1,16*Rrealdim*Lrealdim,1
		!	do j=1,16*Lrealdim*Rrealdim,1
		!		if(abs(fullH(j,i))>1.0D-2) then
		!			write(*,*) fullH(j,i),j,i
		!		end if
		!	end do
		!end do
! ---------------------------------------------------------
! direct diagonalizaiton the full Hamiltonian matrix
	
	m=0
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs+ncharges) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
		m=m+1
		fullH(:,m)=fullH(:,(i-1)*4*Lrealdim+j)
		fullH(m,:)=fullH((i-1)*4*Lrealdim+j,:)
		end if
	end do
	end do
	fullHdummy(1:m,1:m)=fullH(1:m,1:m)
	write(*,*) "direct ngoodstates=",m
	
	allocate(eigenvalue(m),stat=error)
	if(error/=0) stop
	!allocate(goodH(m,m),stat=error)
	!if(error/=0) stop
	!goodH=fullHdummy(1:m,1:m)
	!write(*,*) "fullH",fullH(1:m,1:m)
	!call syevr(fullH(1:m,1:m),eigenvalue)
	!write(*,*) "syevr,direct diagonalizaiton result,E=",eigenvalue
	write(*,*) "fullH"
	write(*,*) fullH(1:m,1:m)
	call syevd(fullHdummy(1:m,1:m),eigenvalue,'V','U',info)
	!'U',z,vl,vu,1,1,m,isuppz,abstol,info)
	write(*,*) "info",info
	write(*,*) "syevd,direct diagonalizaiton result,E=",eigenvalue
	write(*,*) "eigenstate"
	write(*,*) fullHdummy(1:m,1)
	deallocate(eigenvalue)
	!deallocate(goodH)
		write(*,*) "fullmat ends"
	end if
	
	if(myid==0) then
	deallocate(buffermat0)
	deallocate(fullH)
	deallocate(fullHdummy)
	deallocate(I4M)
	end if
	deallocate(buffermat)


return
end subroutine
