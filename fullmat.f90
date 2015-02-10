subroutine fullmat
! this subroutine only used in test mode
! construct the fullmat in the 16M^2 basis

	use mpi
	use variables
	use mathlib
	USE BLAS95
	USE LAPACK95
	USE F95_PRECISION

	implicit none
	integer :: i,j,k,l,m,operaindex,tmp
	integer :: error
	integer :: status(MPI_STATUS_SIZE)
	real(kind=8),allocatable :: buffermat(:,:,:),buffermat0(:,:,:,:),fullH(:,:),fullH2(:,:),fullHdummy(:,:),I4M(:,:),goodH(:,:)
	real(kind=8),allocatable :: eigenvalue(:)!,z(:,:)
	!real(kind=8) :: vl,vu,abstol
	!integer,allocatable :: isuppz(:)
	integer :: info
	integer :: i1,j1
	
	

	if(myid==0) then 
		write(*,*) "enter in fullmat subroutine"
	end if
	!if(myid==0) then
	!read(*,*) a
	!end if
	
	allocate(buffermat(4*subM,4*subM,3),stat=error)
	if(error/=0) stop
	buffermat=0.0D0
	if(myid==0) then
		allocate(buffermat0(4*subM,4*subM,3,norbs),stat=error)
		if(error/=0) stop
		allocate(fullH(16*Rrealdim*Lrealdim,16*Rrealdim*Lrealdim),stat=error)
		if(error/=0) stop
		allocate(fullH2(16*Rrealdim*Lrealdim,16*Rrealdim*Lrealdim),stat=error)
		if(error/=0) stop
		allocate(fullHdummy(16*Rrealdim*Lrealdim,16*Rrealdim*Lrealdim),stat=error)
		if(error/=0) stop
		buffermat0=0.0D0
		fullH=0.0D0
		fullH2=0.0D0
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
			!if(nright==1 .and. j==norbs-1) then
			!	write(*,*) buffermat
			!end if
			call MPI_SEND(buffermat,16*subM*subM*3,mpi_real8,0,j,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(:,:,:,j),16*subM*subM*3,mpi_real8,orbid(j),j,MPI_COMM_WORLD,status,ierr)
		end if
	end do
	


	if(myid==0) then
!debug			
	!	read(*,*) tmp
		open(unit=11,file="imme.tmp",status="old",position="append")
		do i=1,nleft+1,1
		do j=norbs,norbs-nright,-1
			if(bondlink(i,j)==1) then
				open(unit=120,file="imme4.tmp",status="replace")
				do k=1,2,1
					!write(120,*) buffermat0(:,:,k,j)
					!read(*,*) tmp
					buffermat(1:4*Lrealdim,1:4*Lrealdim,k)=transpose(buffermat0(1:4*Lrealdim,1:4*Lrealdim,k,i))
					call directproduct(buffermat(1:Lrealdim*4,1:Lrealdim*4,k),4*Lrealdim,&
					buffermat0(1:Rrealdim*4,1:Rrealdim*4,k,j),4*Rrealdim,fullHdummy)
					do l=1,4*Lrealdim,1
						fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Lrealdim)=&
						fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Lrealdim)*((-1.0D0)**(mod(quantabigL(l,1),2)))
					end do

					!write(11,*) fullHdummy*t(i,j)*(-1.0D0)
				! transfer from al*ar^+ to ar^+*al
					fullHdummy=fullHdummy*(-1.0D0)
				!	write(120,*) fullHdummy
					!write(11,*) "ar+al",k
					!do i1=1,16*Lrealdim*Rrealdim,1
					!do j1=1,16*Lrealdim*Rrealdim,1
					!	if(abs(fullHdummy(j1,i1))>0.1D-3) then
					!	write(11,*) fullHdummy(j1,i1),j1,i1
					!end if
					!end do
					!end do
					fullH=(fullHdummy+transpose(fullHdummy))*t(i,j)+fullH
				end do
				close(120)
				
				do k=1,2,1
					buffermat(1:4*Rrealdim,1:4*Rrealdim,k)=transpose(buffermat0(1:4*Rrealdim,1:4*Rrealdim,k,j))
					call directproduct(buffermat0(1:Lrealdim*4,1:Lrealdim*4,k,i),4*Lrealdim,&
					buffermat(1:Rrealdim*4,1:Rrealdim*4,k),4*Rrealdim,fullHdummy)
					do l=1,4*Lrealdim,1
						fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Lrealdim)=&
						fullHdummy(:,l:16*Lrealdim*Rrealdim:4*Lrealdim)*((-1.0D0)**(mod(quantabigL(l,1),2)))
					end do
					!write(11,*) "al+ar",k
					!do i1=1,16*Lrealdim*Rrealdim,1
					!do j1=1,16*Lrealdim*Rrealdim,1
					!	if(abs(fullHdummy(j1,i1))>0.1D-3) then
					!	write(11,*) fullHdummy(j1,i1),j1,i1
					!end if
					!end do
					!end do
					!write(11,*) fullHdummy*t(i,j) 
					fullH2=(fullHdummy+transpose(fullHdummy))*t(i,j)+fullH2
				end do

				fullHdummy=fullH-fullH2
			!	do i1=1,16*Lrealdim*Rrealdim,1
			!	do j1=1,16*Lrealdim*Rrealdim,1
			!		if(abs(fullHdummy(j1,i1))>0.1D-3) then
			!			write(11,*) fullHdummy(j1,i1),j1,i1
			!		end if
			!	end do
			!	end do
				
			end if
			! PPP term
				call directproduct(buffermat0(1:Lrealdim*4,1:Lrealdim*4,3,i),4*Lrealdim,&
				buffermat0(1:Rrealdim*4,1:Rrealdim*4,3,j),4*Rrealdim,fullHdummy)
				fullH=fullH+fullHdummy*pppV(i,j)
				fullH2=fullH2+fullHdummy*pppV(i,j)
				if(bondlink(i,j)==1) then
				!	write(11,*) fullHdummy*pppV(i,j)
				end if

		end do
		end do
		close(11)
! Hl contribute
		call directproduct(Hbig(1:Lrealdim*4,1:Lrealdim*4,1),4*Lrealdim,&
		I4M(1:Rrealdim*4,1:Rrealdim*4),4*Rrealdim,fullHdummy)
		fullH=fullH+fullHdummy
		fullH2=fullH2+fullHdummy
! HR contribute
		call directproduct(I4M(1:Lrealdim*4,1:Lrealdim*4),4*Lrealdim,&
		Hbig(1:Rrealdim*4,1:Rrealdim*4,2),4*Rrealdim,fullHdummy)
		fullH=fullH+fullHdummy
		fullH2=fullH2+fullHdummy

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
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
		m=m+1
		if((i-1)*4*Lrealdim+j/=m) then
		call swap(fullH(:,(i-1)*4*Lrealdim+j),fullH(:,m))
		call swap(fullH((i-1)*4*Lrealdim+j,:),fullH(m,:))
		end if
		fullH2(:,m)=fullH2(:,(i-1)*4*Lrealdim+j)
		fullH2(m,:)=fullH2((i-1)*4*Lrealdim+j,:)
		end if
	end do
	end do
	fullHdummy(1:m,1:m)=fullH(1:m,1:m)
	write(*,*) "direct ngoodstates=",m
	!open(unit=999,file="H.tmp",status="replace")
	!write(999,*) fullH(1:m,1:m)
	!close(999)


	allocate(eigenvalue(m),stat=error)
	if(error/=0) stop
	!allocate(goodH(m,m),stat=error)
	!if(error/=0) stop
	!goodH=fullHdummy(1:m,1:m)
	!write(*,*) "fullH",fullH(1:m,1:m)
	!call syevr(fullH(1:m,1:m),eigenvalue)
	!write(*,*) "syevr,direct diagonalizaiton result,E=",eigenvalue
write(*,*) "fullH"
!write(*,*) fullH(1:m,1:m)
!do i=1,m,1
!do j=1,m,1
!	if(abs(fullH(j,i))>1.0D-7) then
!		write(*,*) fullH(j,i),j,i
!	end if
!end do
!end do
call syevd(fullH(1:m,1:m),eigenvalue,'V','U',info)
	!'U',z,vl,vu,1,1,m,isuppz,abstol,info)
write(*,*) "info",info
write(*,*) "syevd,direct diagonalizaiton result,energy=",eigenvalue(1:10)
write(*,*) "eigenstate"
write(*,*) fullH(1:m,1)
	!call syevd(fullH2(1:m,1:m),eigenvalue,'V','U',info)
	!'U',z,vl,vu,1,1,m,isuppz,abstol,info)
	!write(*,*) "info",info
	!write(*,*) "syevd2,direct diagonalizaiton result,E=",eigenvalue
	!write(*,*) "eigenstate2"
	!write(*,*) fullH2(1:m,1)
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
