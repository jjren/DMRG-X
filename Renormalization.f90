Subroutine Renormalization(indexLp1,indexRm1,direction)
! after diaganolizaiton we need to renormalizaiton the many body states
! in fact only renormalization all the operator matrix
! direction=l means l block is the system
! direction=i means is the infinit MPS
! direction=r means r block is the system

! indexLp1 is the nleft+1=sigmaL index
! indexRm1 is the indexR-1=sigmaR index

	USE mpi
	USE variables
	use communicate
	USE BLAS95
	USE LAPACK95
	USE F95_PRECISION
    use excitedbase

	implicit none
	
	integer :: i,error,info,j,k,ierr
	integer :: operaindex
	! integer :: mindim
	! mindim is the SVD minimun dimension
	real(kind=8),allocatable ::  leftu(:,:),rightv(:,:),singularvalue(:),&
	dummymat(:,:),leftubuffer(:,:),rightvbuffer(:,:)
	! leftu is the left transfer unitary matrix
	! rightv is the right transfer unitary matrix
	! be careful that leftu is column like,U(+)U=1
	! rightv is row like,vv(+)=1
	! singularvalue
	character(len=1) :: direction 
	integer :: reclength,indexLp1,indexRm1
	logical :: alive
	real(kind=8) :: norm
! exscheme2 space
	integer,allocatable :: quantasmaL2(:,:),quantasmaR2(:,:),symmlinksma2(:,:,:)
	real(kind=8),allocatable ::leftu2(:,:),rightv2(:,:)

! debug
!	integer :: mindim
!	real(kind=8),allocatable :: ww(:)
!	real(kind=8),allocatable :: coeffbuffer(:,:)
	


	if(myid==0)  then
		write(*,*) "enter Renormalization subroutine"
	end if


	if(4*Lrealdim>subM .or. 4*Rrealdim>subM .or. (mode=='d' .and. &
	modeindex==1) ) then
	if(myid==0) then
! when nstate=1, we use SVD method to renormalize the 4M basis
		!if(nstate==1) then
!			mindim=min(4*Lrealdim,4*Rrealdim)
!			allocate(leftu(4*Lrealdim,mindim),stat=error)
!			if(error/=0) stop
!			allocate(singularvalue(mindim),stat=error)
!			if(error/=0) stop
!			allocate(rightv(mindim,4*Rrealdim),stat=error)
!			if(error/=0) stop
!			allocate(ww(mindim-1),stat=error)
!			if(error/=0) stop
!			allocate(coeffbuffer(4*Lrealdim,4*Rrealdim),stat=error)
!			if(error/=0) stop
!			coeffbuffer=coeffIF(1:4*Lrealdim,1:4*Rrealdim,1)
			! gesvd sigular value is the descending order
!			call gesvd(coeffbuffer(1:4*Lrealdim,1:4*Rrealdim),singularvalue,leftu,rightv,ww,'N',info)
!			if(info/=0) then
!				write(*,*) "SVD failed!",indexLp1,indexRm1,info
!			stop
!			end if
!			write(*,*) "singularvalue"
!			write(*,*) singularvalue
		!	write(*,*) leftu
		!	write(*,*) rightv
!			deallocate(singularvalue)
!			deallocate(leftu)
!			deallocate(rightv)
!			deallocate(ww)
!			deallocate(coeffbuffer)

! the reduced density matrix diagnal element is singular value^2
		!	singularvalue=singularvalue*singularvalue
! when nstate/=1, two different scheme to target the excited states
! average method
		if(nstate==1 .or. exscheme==1) then
			allocate(singularvalue(subM),stat=error)
			if(error/=0) stop
!---------------left transfer unitary matrix-------
			allocate(leftu(4*Lrealdim,subM),stat=error)
			if(error/=0) stop
			!allocate(valueL(4*Lrealdim),stat=error)
			!if(error/=0) stop
			!allocate(coeffbufferL(4*Lrealdim,4*Lrealdim,nstate+1),stat=error)
			!if(error/=0) stop
			!coeffbufferL(:,:,nstate+1)=0.0D0
			!do i=1,nstate,1
			!call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffbufferL(:,:,i),'N','T',1.0D0,0.0D0)
			!coeffbufferL(:,:,nstate+1)=coeffbufferL(:,:,nstate+1)+coeffbufferL(:,:,i)*nweight(i)
			!end do
			! syevd eigenvalue is the ascending order
			!call syevd(coeffbufferL(:,:,nstate+1),valueL,'V','U',info)
			!if(info/=0) then
			!	write(*,*) "left diagnolization failed!",indexLp1,info
			!	stop
			!end if
			!leftu=coeffbufferL(:,4*Lrealdim-subM+1:4*Lrealdim,nstate+1)
			!singularvalue=valueL(4*Lrealdim-subM+1:4*Lrealdim)
!---------------right transfer unitary matrix-------
			allocate(rightv(subM,4*Rrealdim),stat=error)
			if(error/=0) stop
			!allocate(valueR(4*Rrealdim),stat=error)
			!if(error/=0) stop
			!allocate(coeffbufferR(4*Rrealdim,4*Rrealdim,nstate+1),stat=error)
			!if(error/=0) stop
			!coeffbufferR(:,:,nstate+1)=0.0D0
			!do i=1,nstate,1
			!call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffbufferR(:,:,i),'T','N',1.0D0,0.0D0)
			!coeffbufferR(:,:,nstate+1)=coeffbufferR(:,:,nstate+1)+coeffbufferR(:,:,i)*nweight(i)
			!end do
			! syevd eigenvalue is the ascending order
			!call syevd(coeffbufferR(:,:,nstate+1),valueR,'V','U',info)
			!if(info/=0) then
			!	write(*,*) "right diagnolization failed!",indexRm1,info
			!	stop
			!end if
			call splitsvdL(singularvalue,leftu,1,nstate,indexlp1)
			!write(*,*) "Lsingular"
			!write(*,*) singularvalue(1:subM)
			call splitsvdR(singularvalue,rightv,1,nstate,indexRm1)
			!write(*,*) "Rsingular"
			!write(*,*) singularvalue(1:subM)
			!rightv=transpose(coeffbufferR(:,4*Rrealdim-subM+1:4*Rrealdim,nstate+1))
			!singularvalue=valueR(4*Rrealdim-subM+1:4*Rrealdim)
		!	write(*,*) "quanta"
		!	write(*,*) quantasmaL
		!	write(*,*) quantasmaR
		!	write(*,*) symmlinksma(:,:,1)
		!	write(*,*) symmlinksma(:,:,2)
			

!--------------------------------------------------
		else if(exscheme==2) then
			allocate(singularvalue(subM*nstate),stat=error)
			if(error/=0) stop
			allocate(leftu2(4*Lrealdim,subM*nstate),stat=error)
			if(error/=0) stop
			allocate(rightv2(subM*nstate,4*Rrealdim),stat=error)
			if(error/=0) stop
			allocate(quantasmaL2(subM*nstate,2),stat=error)
			if(error/=0) stop
			allocate(quantasmaR2(subM*nstate,2),stat=error)
			if(error/=0) stop
			if(logic_spinreversal/=0) then
				allocate(symmlinksma2(subM*nstate,1,2),stat=error)
				if(error/=0) stop
			end if

			write(*,*) "my new method to calulate the excited states!"
			do i=1,nstate,1
				call splitsvdL(singularvalue((i-1)*subM+1:i*subM),leftu2(:,(i-1)*subM+1:i*subM),i,i,indexlp1)
				call splitsvdR(singularvalue((i-1)*subM+1:i*subM),rightv2((i-1)*subM+1:i*subM,:),i,i,indexRm1)
				quantasmaL2((i-1)*subM+1:i*subM,:)=quantasmaL
				quantasmaR2((i-1)*subM+1:i*subM,:)=quantasmaR
				if(logic_spinreversal/=0) then
					symmlinksma2((i-1)*subM+1:i*subM,:,:)=symmlinksma
				end if
			end do
			
			allocate(leftu(4*Lrealdim,subM),stat=error)
			if(error/=0) stop
			allocate(rightv(subM,4*Rrealdim),stat=error)
			if(error/=0) stop

			if(logic_spinreversal==0) then
				call excitedbasis(leftu,rightv,singularvalue,leftu2,rightv2,quantasmaL2,quantasmaR2)
			else
				call excitedbasis(leftu,rightv,singularvalue,leftu2,rightv2,quantasmaL2,quantasmaR2,symmlinksma2)
			end if

			deallocate(leftu2)
			deallocate(rightv2)
			deallocate(quantasmaL2)
			deallocate(quantasmaR2)
			if(logic_spinreversal/=0) then
				deallocate(symmlinksma2)
			end if
		end if
!--------------------------------------------------------------------------
! there are some numerical inaccuracy when we do svd or diagonalizaiton and this
! will make different good quantum number stat to mix. This is wrong. So we need
! to Reconstruct the transfer matrix leftu and rightv
	!	do i=1,subM,1
	!		do j=1,4*Lrealdim,1
	!			if(leftu(j,i)>quantaconst) then
	!				do k=1,4*Lrealdim,1
	!					if(quantabigL(k,1)/=quantabigL(j,1) .or. &
	!					quantabigL(k,2)/=quantabigL(j,2)) then
	!					leftu(k,i)=0.0D0
	!					end if
	!				end do
	!				exit
	!			end if
	!		end do
	!		norm=dot(leftu(:,i),leftu(:,i))
	!		if(norm<1.0D-10) then
	!			write(*,*) "-------------------------------------"
	!			write(*,*) "Renormalizaiton norm<1.0D-10 caution!"
	!			write(*,*) "-------------------------------------"
	!		end if
	!		leftu(:,i)=leftu(:,i)/sqrt(norm)
	!	end do
		
	!	do i=1,subM,1
	!		do j=1,4*Rrealdim,1
	!			if(rightv(i,j)>quantaconst) then
	!				do k=1,4*Rrealdim,1
	!					if(quantabigL(k,1)/=quantabigL(j,1) .or. &
	!					quantabigL(k,2)/=quantabigL(j,2)) then
	!					rightv(i,k)=0.0D0
	!					end if
	!				end do
	!				exit
	!			end if
	!		end do
	!		norm=dot(rightv(i,:),rightv(i,:))
	!		rightv(i,:)=rightv(i,:)/sqrt(norm)
	!	end do

!--------------------------------------------------
! write the wavefunction matrix---------------------------------------
		if(modeindex/=1) then
		reclength=2*subM*subM
		inquire(file="wavefunction.tmp",exist=alive)
		if(alive) then
			open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		open(unit=106,file="singularvalue.tmp",status="replace")
		if(direction=='i' .or. direction=='l') then
			! divid the leftu and rightv to small M*M matrix
			do i=1,4,1
			write(105,rec=4*(IndexLp1-1)+i) leftu((i-1)*Lrealdim+1:i*Lrealdim,1:subM)
			end do
			if(direction=='l') then
				do i=1,4,1
				write(105,rec=4*indexLp1+i) rightv(1:subM,i:4*Rrealdim:4)
				end do
			end if
		end if
		if(direction=='i' .or. direction=='r') then
			do i=1,4,1
			write(105,rec=4*(indexRm1-1)+i) rightv(1:subM,i:4*Rrealdim:4)
			end do
			if(direction=='r') then
				do i=1,4,1
				write(105,rec=4*(indexRm1-2)+i) leftu((i-1)*Lrealdim+1:i*Lrealdim,1:subM)
				end do
			end if
		end if
		! write singularvalue though only used in finit MPS
		! the singularvalue we use here is the exactly singlarvalue^2
		write(106,*) singularvalue
		close(105)
		close(106)
		end if
		
	end if

	call MPI_bcast(quantasmaL(:,:),subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(quantasmaR(:,:),subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
!------------------------------------------------------------
	
!   the two buffer are used to transfer the unitary matrix to other process
! remember test good quantum number
!------------------------------------------------------------------
		allocate(dummymat(subM,4*subM),stat=error)
		if(error/=0) stop
		
	if(direction=='i' .or. direction=='l') then
		
		allocate(leftubuffer(4*Lrealdim,subM),stat=error)
		if(error/=0) stop
		if(myid==0) then
		leftubuffer=leftu
		end if
		
		call MPI_BCAST(leftubuffer,4*Lrealdim*subM,mpi_real8,0,MPI_COMM_WORLD,ierr)
! left space operator renormalization
		do i=1,indexLp1,1
			if(myid==orbid(i)) then
				if(mod(i,nprocs-1)==0) then
					operaindex=i/(nprocs-1)
				else
					operaindex=i/(nprocs-1)+1
				end if
					do j=1,3,1
					call gemm(leftubuffer,operamatbig(1:4*Lrealdim,1:4*Lrealdim,3*(operaindex-1)+j),&
					dummymat(:,1:4*Lrealdim),'T','N',1.0D0,0.0D0)
					call gemm(dummymat(:,1:4*Lrealdim),leftubuffer,operamatsma(:,:,3*(operaindex-1)+j),'N','N',1.0D0,0.0D0)
					end do
			end if
		end do
! left HL renormalizaiton,adaptedmat renormalization
		if(myid==0) then
			call gemm(leftubuffer,Hbig(1:4*Lrealdim,1:4*Lrealdim,1),dummymat(:,1:4*Lrealdim),'T','N',1.0D0,0.0D0)
			call gemm(dummymat(:,1:4*Lrealdim),leftubuffer,Hsma(:,:,1),'N','N',1.0D0,0.0D0)
			!if(logic_spinreversal/=0) then
			!	call gemm(leftubuffer,adaptedbig(1:4*Lrealdim,1:4*Lrealdim,1),dummymat(:,1:4*Lrealdim),'T','N',1.0D0,0.0D0)
			!	call gemm(dummymat(:,1:4*Lrealdim),leftubuffer,adaptedsma(:,:,1),'N','N',1.0D0,0.0D0)
			!end if
		end if
! good quantum number renormalization
!		do i=1,subM,1
!			do j=1,4*Lrealdim,1
!				if (leftubuffer(j,i)>=quantaconst) then
!					quantasmaL(i,:)=quantabigL(j,:)
!					exit
!				end if
!			end do
!		end do
	end if
! right space
	if(direction=='i' .or. direction=='r') then
		
		allocate(rightvbuffer(subM,4*Rrealdim),stat=error)
		if(error/=0) stop
		if(myid==0) then
		rightvbuffer=rightv
		end if
		
		call MPI_BCAST(rightvbuffer,4*Rrealdim*subM,mpi_real8,0,MPI_COMM_WORLD,ierr)
! right space operator
		do i=norbs,indexRm1,-1
			if(myid==orbid(i)) then
				if(mod(i,nprocs-1)==0) then
					operaindex=i/(nprocs-1)
				else
					operaindex=i/(nprocs-1)+1
				end if
					do j=1,3,1
					call gemm(rightvbuffer,operamatbig(1:4*Rrealdim,1:4*Rrealdim,3*(operaindex-1)+j),dummymat(:,1:4*Rrealdim),'N','N',1.0D0,0.0D0)
					call gemm(dummymat(:,1:4*Rrealdim),rightvbuffer,operamatsma(:,:,3*(operaindex-1)+j),'N','T',1.0D0,0.0D0)
					end do
			end if
		end do
! right space HR,adaptedmat
		if(myid==0) then
			call gemm(rightvbuffer,Hbig(1:4*Rrealdim,1:4*Rrealdim,2),dummymat(:,1:4*Rrealdim),'N','N',1.0D0,0.0D0)
			call gemm(dummymat(:,1:4*Rrealdim),rightvbuffer,Hsma(:,:,2),'N','T',1.0D0,0.0D0)
			!if(logic_spinreversal/=0) then
			!	call gemm(rightvbuffer,adaptedbig(1:4*Rrealdim,1:4*Rrealdim,2),dummymat(:,1:4*Rrealdim),'N','N',1.0D0,0.0D0)
			!	call gemm(dummymat(:,1:4*Rrealdim),rightvbuffer,adaptedsma(:,:,2),'N','T',1.0D0,0.0D0)
			!end if
		end if
	!	do i=1,subM,1
	!		do j=1,4*Rrealdim,1
	!			if (rightvbuffer(i,j)>=quantaconst) then
	!				quantasmaR(i,:)=quantabigR(j,:)
	!				exit
	!			end if
	!		end do
	!	end do
		end if
	end if

!-----------------------------------------------------------------------------------------------------------------------
! when we need not renormalization
! operamatbig -> operamatsma 
! Hbig -> Hsma
! adaptedbig -> adaptedsma
! symmlinkbig -> symmlinksma
! only in the infinit MPS process
if(4*Lrealdim<=subM .and. 4*Rrealdim<=subM .and. mode/='d') then
	do i=1,indexLp1,1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			operamatsma(1:4*Lrealdim,1:4*Lrealdim,3*(operaindex-1)+1:3*operaindex)=&
			operamatbig(1:4*Lrealdim,1:4*Lrealdim,3*(operaindex-1)+1:3*operaindex)
		end if
	end do
	
	do i=norbs,indexRm1,-1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			
			operamatsma(1:4*Rrealdim,1:4*Rrealdim,3*(operaindex-1)+1:3*operaindex)=&
			operamatbig(1:4*Rrealdim,1:4*Rrealdim,3*(operaindex-1)+1:3*operaindex)
		end if
	end do
	
		quantasmaL(1:4*Lrealdim,:)=quantabigL(1:4*Lrealdim,:)
		quantasmaR(1:4*Rrealdim,:)=quantabigR(1:4*Rrealdim,:)

	if(myid==0) then
		Hsma(1:4*Lrealdim,1:4*Lrealdim,1)=Hbig(1:4*Lrealdim,1:4*Lrealdim,1)
		Hsma(1:4*Rrealdim,1:4*Rrealdim,2)=Hbig(1:4*Rrealdim,1:4*Rrealdim,2)
		if(logic_spinreversal/=0) then
		!	adaptedsma(1:4*Lrealdim,1:4*Lrealdim,1)=adaptedbig(1:4*Lrealdim,1:4*Lrealdim,1)
		!	adaptedsma(1:4*Rrealdim,1:4*Rrealdim,2)=adaptedbig(1:4*Rrealdim,1:4*Rrealdim,2)
			symmlinksma(1:4*Lrealdim,1,1)=symmlinkbig(1:4*Lrealdim,1,1)
			symmlinksma(1:4*Rrealdim,1,2)=symmlinkbig(1:4*Rrealdim,1,2)

		end if
	end if

end if

if(4*Lrealdim>subM .or. 4*Rrealdim>subM .or. (mode=='d' .and. &
modeindex==1)) then
if(myid==0) then
	deallocate(singularvalue)
	deallocate(leftu)
	deallocate(rightv)
end if
	if(direction=='i' .or. direction=='l') then
	deallocate(leftubuffer)
	end if
	if(direction=='i' .or. direction=='r') then
	deallocate(rightvbuffer)
	end if
	deallocate(dummymat)
end if



!if(myid==0) then
!	deallocate(coeffIF)
!end if

return

end Subroutine Renormalization

