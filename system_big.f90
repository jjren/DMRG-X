module Construct_System_Big

contains
!=================================================================
!=================================================================

Subroutine System_Big(domain,&
cap_sma,cap_smacol,cap_smarow,&
cap_big,cap_bigcol,cap_bigrow,&
cap_Hsma,cap_Hsmacol,cap_Hsmarow,&
cap_Hbig,cap_Hbigcol,cap_Hbigrow,&
cap_quantasmaL,cap_quantasmaR,&
LRdim,isubM,ifperturbation)

! construct the L+sigmaL/R+sigmaR subspace operator matrix in 4M basis
! two different conditions
! R space and logic_C2==0  ;  L space or (R space and logic_C2/=0)
! all matrix in CSR sparse matrix format
	
	use variables
	use mpi
	use mathlib
	use symmetry
	use communicate
	use exit_mod
	use module_sparse
	use OnesiteMatrix
	use BLAS95
	use f95_precision

	implicit none

	character(len=1),intent(in) :: domain   !L/R
	real(kind=r8),intent(in) :: cap_sma(:,:),cap_Hsma(:,:)
	real(kind=r8),intent(out) :: cap_big(:,:),cap_Hbig(:,:)
	integer(kind=i4),intent(in) :: cap_smacol(:,:),cap_smarow(:,:),&
		cap_Hsmacol(:,:),cap_Hsmarow(:,:),cap_quantasmaL(:,:),cap_quantasmaR(:,:)
	integer(kind=i4),intent(out) :: cap_bigcol(:,:),cap_bigrow(:,:),&
		cap_Hbigcol(:,:),cap_Hbigrow(:,:)
	integer,intent(in) :: LRdim,isubM
	logical,intent(in) :: ifperturbation
	
	! local
	integer :: orbstart,orbend,orbadd,Hindex
	! orbstart is from 1 or norbs-nright+1
	! orbend is nleft or norbs
	! orbadd is nleft+1 or norbs-nright
	! Hindex : HL/1 HR/2
	integer :: iHbigdim ,& ! the biggest H dimension
			ibigdim1
	integer :: operaindex,operaindex2,operaindex3
	real(kind=r8) :: II(4),IM(isubM)
	integer(kind=i4) :: IIrowindex(5),IIcolindex(4),IMrowindex(isubM+1),IMcolindex(isubM)

	integer(kind=i4),allocatable :: phase(:)

	real(kind=r8),allocatable :: Hbufmat(:),H0mat(:)  ! Hbuffer
	integer(kind=i4),allocatable :: Hbufrowindex(:),Hbufcolindex(:),&
						H0rowindex(:),H0colindex(:)
	
	integer :: nnonzero
	logical :: havecount
	integer :: count1,shouldrecv(nprocs-1)
	integer :: i,j,k
	integer :: error

	logical :: ifbondord,iflocalspin
	integer :: onesitematindex,systemmatindex

	integer :: status(MPI_STATUS_SIZE),sendrequest  ! MPI flag
	character(len=1),allocatable :: packbuf(:)
	integer(kind=i4) :: packsize,position1
	integer :: ierr

	call master_print_message("enter in subroutine system_big")
	
	! define the biggest H dimension
	if(ifperturbation==.true.) then
		iHbigdim=Hbigdimp
		ibigdim1=bigdim1p
	else
		iHbigdim=Hbigdim
		ibigdim1=bigdim1
	end if

! set the parameters
	if(domain=='L') then
		orbstart=1
		orbend=nleft
		orbadd=nleft+1
		Hindex=1
		if(nleft<=(norbs-1)/2) then
			ifbondord=.true.
			iflocalspin=.true.
		else
			ifbondord=.false.
			iflocalspin=.false.
		end if
	else if(domain=='R') then
		orbstart=norbs-nright+1
		orbend=norbs
		orbadd=norbs-nright
		Hindex=2
		if(nright<=(norbs-2)/2) then
			ifbondord=.true.
			iflocalspin=.true.
		else
			ifbondord=.false.
			iflocalspin=.false.
		end if
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
	end if

! allocate workarray
	
	do i=orbstart,orbend,1
		if(myid==orbid1(i,1)) then
		! the intermediate H0mat matrix
			allocate(Hbufmat(iHbigdim),stat=error)
			if(error/=0) stop
			allocate(Hbufrowindex(4*isubM+1),stat=error)
			if(error/=0) stop
			allocate(Hbufcolindex(iHbigdim),stat=error)
			if(error/=0) stop
			Hbufrowindex=1

			! store the H0mat matrix
			allocate(H0mat(iHbigdim),stat=error)
			if(error/=0) stop
			allocate(H0colindex(iHbigdim),stat=error)
			if(error/=0) stop
			allocate(H0rowindex(4*isubM+1),stat=error)
			if(error/=0) stop
			H0rowindex=1
			
			! pack the H0mat matrix
			packsize=iHbigdim*12+16*isubM+1000
			allocate(packbuf(packsize),stat=error) 
			if(error/=0) stop
			
			exit
		end if
	end do

	if(myid==0) then
		allocate(H0mat(iHbigdim),stat=error)
		if(error/=0) stop
		allocate(H0colindex(iHbigdim),stat=error)
		if(error/=0) stop
		allocate(H0rowindex(4*isubM+1),stat=error)
		if(error/=0) stop
		H0rowindex=1
		! pack the H0mat matrix
		packsize=iHbigdim*12+16*isubM+1000
		allocate(packbuf(packsize),stat=error) 
		if(error/=0) stop
	end if

!==================================================================
! set phase value
	allocate(phase(4*isubM),stat=error)
	if(error/=0) stop
	if(domain=='R' .and. logic_C2==0) then
		phase(1:4*LRdim:4)=1
		phase(2:4*LRdim:4)=-1
		phase(3:4*LRdim:4)=-1
		phase(4:4*LRdim:4)=1
	else
		do j=1,4*LRdim,1
			if(mod(j,LRdim)==0) then
				k=LRdim
			else
				k=mod(j,LRdim)
			end if
			if(domain=='L') then
				phase(j)=(-1)**(mod(cap_quantasmaL(k,1),2))
			else
				phase(j)=(-1)**(mod(cap_quantasmaR(k,1),2))
			end if
		end do
	end if

! construct the unit matrix
	II=1.0D0
	do i=1,4,1
		IIrowindex(i)=i
		IIcolindex(i)=i
	end do
	IIrowindex(5)=5
	IM=1.0D0
	do i=1,isubM,1
		IMrowindex(i)=i
		IMcolindex(i)=i
	end do
	IMrowindex(isubM+1)=isubM+1

!=========================================================================
	! calculate the hopping term and PPP term
	
	H0rowindex=1
	
	do i=orbstart,orbend,1
	if(myid==orbid1(i,1)) then
		!transfer integral term
		if(bondlink(i,orbadd)==1) then
			do j=1,2,1
				Hbufrowindex=0
				operaindex=orbid1(i,2)*3-3+j
				if(domain=='R' .and. logic_C2==0) then
					call SparseDirectProduct(4,4,onesitemat(:,j+3),osmcolindex(:,j+3),osmrowindex(:,j+3),&
									LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
									Hbufmat,Hbufcolindex,Hbufrowindex,iHbigdim)
					do k=1,Hbufrowindex(4*LRdim+1)-1,1
						Hbufmat(k)=Hbufmat(k)*DBLE(phase(Hbufcolindex(k)))*(-1.0D0)
					end do
				else
					call SparseDirectProduct(LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
									4,4,onesitemat(:,j+3),osmcolindex(:,j+3),osmrowindex(:,j+3),&
									Hbufmat,Hbufcolindex,Hbufrowindex,iHbigdim)
					do k=1,Hbufrowindex(4*LRdim+1)-1,1
						Hbufmat(k)=Hbufmat(k)*DBLE(phase(Hbufcolindex(k)))
					end do
				end if
				!======================================================================
				! store the bondorder mat  ai^+aj ( i is near the boundary)
				! hopping operator is the intermediate of PPP model
				if(logic_bondorder/=0 .and. ifbondord==.true.) then
					if(myid/=orbid2(i,orbadd,1)) then
						write(*,*) "=============================================="
						write(*,*) "myid/=orbid2(i,orbadd,1)",myid,orbid2(i,orbadd,1)
						write(*,*) "=============================================="
						stop
					end if
					operaindex2=orbid2(i,orbadd,2)*2-2+j
					bigrowindex2(1:4*LRdim+1,operaindex2)=Hbufrowindex(1:4*LRdim+1)
					nnonzero=Hbufrowindex(4*LRdim+1)-1
					bigcolindex2(1:nnonzero,operaindex2)=Hbufcolindex(1:nnonzero)
					operamatbig2(1:nnonzero,operaindex2)=Hbufmat(1:nnonzero)
				end if
				!======================================================================

				! ai^+*aj
				call SpmatAdd(4*LRdim,4*LRdim,H0mat,H0colindex,H0rowindex,&
				'N',t(i,orbadd),4*LRdim,4*LRdim,Hbufmat,Hbufcolindex,Hbufrowindex,iHbigdim)
				! aj^+*ai
				call SpmatAdd(4*LRdim,4*LRdim,H0mat,H0colindex,H0rowindex,&
				'T',t(i,orbadd),4*LRdim,4*LRdim,Hbufmat,Hbufcolindex,Hbufrowindex,iHbigdim)
			end do
		end if
		
		if(logic_PPP==1) then
			! ppp term
			Hbufrowindex=0
			operaindex=orbid1(i,2)*3
			if(domain=='R' .and. logic_C2==0) then
				call SparseDirectProduct(4,4,onesitemat(:,3),osmcolindex(:,3),osmrowindex(:,3),&
								LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
								Hbufmat,Hbufcolindex,Hbufrowindex,iHbigdim)
			else
				call SparseDirectProduct( LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
								4,4,onesitemat(:,3),osmcolindex(:,3),osmrowindex(:,3),&
								Hbufmat,Hbufcolindex,Hbufrowindex,iHbigdim)
			end if
			call SpmatAdd(4*LRdim,4*LRdim,H0mat,H0colindex,H0rowindex,&
				'N',pppV(i,orbadd),4*LRdim,4*LRdim,Hbufmat,Hbufcolindex,Hbufrowindex,iHbigdim)
		end if
	end if
	end do
	
	! every process pack the H0mat mat and send to 0 process
	do j=orbstart,orbend,1
		if(myid==orbid1(j,1)) then
			position1=0
			call MPI_PACK(H0rowindex,4*isubM+1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(H0colindex,H0rowindex(4*LRdim+1)-1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(H0mat,H0rowindex(4*LRdim+1)-1,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			call MPI_ISEND(packbuf,position1,MPI_PACKED,0,myid,MPI_COMM_WORLD,sendrequest,ierr)
			exit  ! only send once
		end if
	end do

!==============================================================================================================

	! construct the L/R(without sigmaL/R) subspace operator matrix in 4M basis
	do i=orbstart,orbend,1
	if(myid==orbid1(i,1)) then
		
		cap_bigrow(:,orbid1(i,2)*3-2:orbid1(i,2)*3)=0
		
		if(domain=='R' .and. logic_C2==0 ) then
			do j=1,3,1
				operaindex=(orbid1(i,2)-1)*3+j
				call SparseDirectProduct(4,4,II,IIcolindex,IIrowindex,&
								LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
								cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),ibigdim1)
				if(j<=2) then
					do k=1,cap_bigrow(4*LRdim+1,operaindex)-1,1
						cap_big(k,operaindex)=cap_big(k,operaindex)*DBLE(phase(cap_bigcol(k,operaindex)))
					end do
				end if
			end do
		else
			do j=1,3,1
				operaindex=(orbid1(i,2)-1)*3+j
				call SparseDirectProduct(LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
								4,4,II,IIcolindex,IIrowindex,&
								cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),ibigdim1)
			end do
		end if
	end if
	end do

!==============================================================================================================

! construct the sigmaL/R subspace operator matrix in 4M basis
	
	if(myid==orbid1(orbadd,1)) then
		
		cap_bigrow(:,orbid1(orbadd,2)*3-2:orbid1(orbadd,2)*3)=0
		
		if(domain=='R' .and. logic_C2==0) then
			do j=1,3,1
				operaindex=(orbid1(orbadd,2)-1)*3+j
				call SparseDirectProduct(4,4,onesitemat(:,j),osmcolindex(:,j),osmrowindex(:,j),&
								LRdim,LRdim,IM,IMcolindex,IMrowindex,&
								cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),ibigdim1)
			end do
		else 
			do j=1,3,1
				operaindex=(orbid1(orbadd,2)-1)*3+j
				call SparseDirectProduct(LRdim,LRdim,IM,IMcolindex,IMrowindex,&
								4,4,onesitemat(:,j),osmcolindex(:,j),osmrowindex(:,j),&
								cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),ibigdim1)
				if(j<=2) then
					do k=1,cap_bigrow(4*LRdim+1,operaindex)-1,1
						cap_big(k,operaindex)=cap_big(k,operaindex)*DBLE(phase(cap_bigcol(k,operaindex)))
					end do
				end if
			end do
		end if
	end if

	! cosntruct the L+sigmaL/R+sigmaR Hamiltonian operator in 4M basis
	if(myid==0) then
		count1=0    ! 0 process should recv how many times
		shouldrecv=0
		do i=orbstart,orbend,1
			havecount=.false.
			do j=1,count1,1
				if(orbid1(i,1)==shouldrecv(j)) then
					havecount=.true.
					exit
				end if
			end do
			if(havecount==.false.) then
				count1=count1+1
				shouldrecv(count1)=orbid1(i,1)
			end if
		end do

		! initiate the cap_Hbig matrix
		cap_Hbigrow(:,Hindex)=1
		
		! get the ppp term and transfer integral term
		do i=1,count1,1
			call MPI_RECV(packbuf,packsize,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
			position1=0
			call MPI_UNPACK(packbuf,packsize,position1,H0rowindex,4*isubM+1,MPI_integer4,MPI_COMM_WORLD,ierr)
			call MPI_UNPACK(packbuf,packsize,position1,H0colindex,H0rowindex(4*LRdim+1)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
			call MPI_UNPACK(packbuf,packsize,position1,H0mat,H0rowindex(4*LRdim+1)-1,MPI_real8,MPI_COMM_WORLD,ierr)
			
			call SpmatAdd(4*LRdim,4*LRdim,cap_Hbig(:,Hindex),cap_Hbigcol(:,Hindex),cap_Hbigrow(:,Hindex),&
			'N',1.0D0,4*LRdim,4*LRdim,H0mat,H0colindex,H0rowindex,iHbigdim)
		end do

!===========================================================

	!     L/R Hamiltonian contribute
		H0rowindex=1

		if(domain=='R' .and. logic_C2==0) then
			call SparseDirectProduct(4,4,II,IIcolindex,IIrowindex,&
							LRdim,LRdim,cap_Hsma(:,Hindex),cap_Hsmacol(:,Hindex),cap_Hsmarow(:,Hindex),&
							H0mat,H0colindex,H0rowindex,iHbigdim)
		else
			call SparseDirectProduct(LRdim,LRdim,cap_Hsma(:,Hindex),cap_Hsmacol(:,Hindex),cap_Hsmarow(:,Hindex),&
							4,4,II,IIcolindex,IIrowindex,&
							H0mat,H0colindex,H0rowindex,iHbigdim)
		end if
			
		call SpmatAdd(4*LRdim,4*LRdim,cap_Hbig(:,Hindex),cap_Hbigcol(:,Hindex),cap_Hbigrow(:,Hindex),&
			'N',1.0D0,4*LRdim,4*LRdim,H0mat,H0colindex,H0rowindex,iHbigdim)

!===========================================================

	!     sigmaL Hamiltonian contribute. site energy+HubbardU

		H0rowindex=1

		if(domain=='R' .and. logic_C2==0) then
			call SparseDirectProduct(4,4,onesitemat(:,6),osmcolindex(:,6),osmrowindex(:,6),&
							LRdim,LRdim,IM,IMcolindex,IMrowindex,&
							H0mat,H0colindex,H0rowindex,iHbigdim)
		else
			call SparseDirectProduct(LRdim,LRdim,IM,IMcolindex,IMrowindex,&
							4,4,onesitemat(:,6),osmcolindex(:,6),osmrowindex(:,6),&
							H0mat,H0colindex,H0rowindex,iHbigdim)
		end if

		call SpmatAdd(4*LRdim,4*LRdim,cap_Hbig(:,Hindex),cap_Hbigcol(:,Hindex),cap_Hbigrow(:,Hindex),&
			'N',1.0D0,4*LRdim,4*LRdim,H0mat,H0colindex,H0rowindex,iHbigdim)
		
!=========================================================================
	
	! construct the symmmlinkbig
		if(logic_spinreversal/=0) then
			call Creatsymmlinkbig(LRdim,domain,Hindex)
		end if
		
	end if

!===================================================================
! bond order matrix construct
! no phase in two operator matrix
	if(logic_bondorder/=0 .and. ifbondord==.true.) then
		! transfer the L/R space(without sigmaL/R) from M basis to 4M basis
		do i=orbstart,orbend,1
		do j=i,orbend,1
			! two conditions bond order term and one partical density matrix
			if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
				if(myid==orbid2(i,j,1)) then
					do k=1,2,1
						operaindex2=orbid2(i,j,2)*2-2+k
						bigrowindex2(:,operaindex2)=0
						if(domain=='R' .and. logic_C2==0) then
							call SparseDirectProduct(4,4,II,IIcolindex,IIrowindex,&
								LRdim,LRdim,operamatsma2(:,operaindex2),smacolindex2(:,operaindex2),smarowindex2(:,operaindex2),&
								operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2),bigdim2)
						else
							call SparseDirectProduct(LRdim,LRdim,operamatsma2(:,operaindex2),smacolindex2(:,operaindex2),smarowindex2(:,operaindex2),&
								4,4,II,IIcolindex,IIrowindex,&
								operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2),bigdim2)
						end if
					end do
				end if
			end if
		end do
		end do

		! construct the sigmaL/R subspace operator matrix in 4M basis
		if(myid==orbid2(orbadd,orbadd,1)) then
			do j=1,2,1
				operaindex2=(orbid2(orbadd,orbadd,2)-1)*2+j
				bigrowindex2(:,operaindex2)=0
				if(domain=='R' .and. logic_C2==0) then
					! in onesite matrix the niup/down index is 7/8
					call SparseDirectProduct(4,4,onesitemat(:,j+6),osmcolindex(:,j+6),osmrowindex(:,j+6),&
							LRdim,LRdim,IM,IMcolindex,IMrowindex,&
							operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2),bigdim2)
				else
					call SparseDirectProduct(LRdim,LRdim,IM,IMcolindex,IMrowindex,&
							4,4,onesitemat(:,j+6),osmcolindex(:,j+6),osmrowindex(:,j+6),&
							operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2),bigdim2)
				end if
			end do
		end if

		! contruct the L/R+sigmaL/R one partical operator matrix
		if(logic_bondorder==2) then
			do i=orbstart,orbend,1
				if(myid==orbid2(i,orbadd,1)) then
					if(myid/=orbid1(i,1)) then
						write(*,*) "=============================================="
						write(*,*) "myid/=orbid1(i,1)",myid,orbid1(i,1)
						write(*,*) "=============================================="
						stop
					end if
					if(bondlink(i,orbadd)==0) then ! bondorder term has calculated 
						do j=1,2,1
							operaindex2=orbid2(i,orbadd,2)*2-2+j
							operaindex=orbid1(i,2)*3-3+j
							bigrowindex2(:,operaindex2)=0
							if(domain=='R' .and. logic_C2==0) then
								call SparseDirectProduct(4,4,onesitemat(:,j+3),osmcolindex(:,j+3),osmrowindex(:,j+3),&
												LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
												operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2),bigdim2)
								do k=1,bigrowindex2(4*LRdim+1,operaindex2)-1,1
									operamatbig2(k,operaindex2)=operamatbig2(k,operaindex2)*DBLE(phase(bigcolindex2(k,operaindex2)))*(-1.0D0)
								end do
							else
								call SparseDirectProduct(LRdim,LRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),&
												4,4,onesitemat(:,j+3),osmcolindex(:,j+3),osmrowindex(:,j+3),&
												operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2),bigdim2)
								do k=1,bigrowindex2(4*LRdim+1,operaindex2)-1,1
									operamatbig2(k,operaindex2)=operamatbig2(k,operaindex2)*DBLE(phase(bigcolindex2(k,operaindex2)))
								end do
							end if
						end do
					end if
				end if
			end do
		end if
	end if
!===================================================================
! local spin operator matrix contruct
	if(logic_localspin==1 .and. iflocalspin==.true.) then
		
		! transfer the L/R space(without sigmaL/R) from M basis to 4M basis
		do i=orbstart,orbend,1
		do j=i,orbend,1
			if(myid==orbid3(i,j,1)) then
				operaindex3=orbid3(i,j,2)
				do k=operaindex3-1,operaindex3,1
					if(domain=='R' .and. logic_C2==0) then
						call SparseDirectProduct(4,4,II,IIcolindex,IIrowindex,&
							LRdim,LRdim,operamatsma3(:,k),smacolindex3(:,k),smarowindex3(:,k),&
							operamatbig3(:,k),bigcolindex3(:,k),bigrowindex3(:,k),bigdim3)
					else
						call SparseDirectProduct(LRdim,LRdim,operamatsma3(:,k),smacolindex3(:,k),smarowindex3(:,k),&
							4,4,II,IIcolindex,IIrowindex,&
							operamatbig3(:,k),bigcolindex3(:,k),bigrowindex3(:,k),bigdim3)
					end if
				end do
			end if
		end do
		end do

		! construct the sigmaL/R subspace operator matrix in 4M basis
		if(myid==orbid3(orbadd,orbadd,1)) then
			do j=1,2,1
				operaindex3=orbid3(orbadd,orbadd,2)-2+j
				if(j==1) then
					onesitematindex=10
				else
					onesitematindex=11
				end if
				if(domain=='R' .and. logic_C2==0) then
					call SparseDirectProduct(4,4,onesitemat(:,onesitematindex),osmcolindex(:,onesitematindex),osmrowindex(:,onesitematindex),&
							LRdim,LRdim,IM,IMcolindex,IMrowindex,&
							operamatbig3(:,operaindex3),bigcolindex3(:,operaindex3),bigrowindex3(:,operaindex3),bigdim3)
				else
					call SparseDirectProduct(LRdim,LRdim,IM,IMcolindex,IMrowindex,&
							4,4,onesitemat(:,onesitematindex),osmcolindex(:,onesitematindex),osmrowindex(:,onesitematindex),&
							operamatbig3(:,operaindex3),bigcolindex3(:,operaindex3),bigrowindex3(:,operaindex3),bigdim3)
				end if
			end do
		end if

		! construct the L/R+sigmaL/sigmaR operator matrix
		! like PPP term no phase
		do i=orbstart,orbend,1
			if(myid==orbid3(i,orbadd,1)) then
				if(orbid3(i,orbadd,1)/=orbid2(i,i,1) .or. orbid3(i,orbadd,1)/=orbid3(i,i,1)) then
					write(*,*) "========================================================"
					write(*,*) "orbid3(i,orbadd,1)/=orbid2(i,i,1) orbid3(i,i,1) failed!",orbid3(i,orbadd,1),orbid2(i,i,1),orbid3(i,i,1)
					write(*,*) "========================================================"
					stop
				end if

				operaindex3=orbid3(i,orbadd,2)
				operaindex2=orbid2(i,i,2)*2
				systemmatindex=orbid3(i,i,2)-1
				
				if(domain=='R' .and. logic_C2==0) then
					! ai^+down*aiup*aj^+up*aj^down
					call SparseDirectProduct(4,4,onesitemat(:,9),osmcolindex(:,9),osmrowindex(:,9),&
							LRdim,LRdim,operamatsma3(:,systemmatindex),smacolindex3(:,systemmatindex),smarowindex3(:,systemmatindex),&
							operamatbig3(:,operaindex3-1),bigcolindex3(:,operaindex3-1),bigrowindex3(:,operaindex3-1),bigdim3)
					! (niup-nidown)*(njup-njdown)
					call SparseDirectProduct(4,4,onesitemat(:,8),osmcolindex(:,8),osmrowindex(:,8),&
							LRdim,LRdim,operamatsma2(:,operaindex2),smacolindex2(:,operaindex2),smarowindex2(:,operaindex2),&
							operamatbig3(:,operaindex3),bigcolindex3(:,operaindex3),bigrowindex3(:,operaindex3),bigdim3)
				else
					! ai^+down*aiup*aj^+up*aj^down
					call SparseDirectProduct( LRdim,LRdim,operamatsma3(:,systemmatindex),smacolindex3(:,systemmatindex),smarowindex3(:,systemmatindex),&
							4,4,onesitemat(:,9),osmcolindex(:,9),osmrowindex(:,9),&
							operamatbig3(:,operaindex3-1),bigcolindex3(:,operaindex3-1),bigrowindex3(:,operaindex3-1),bigdim3)
					! (niup-nidown)*(njup-njdown)
					call SparseDirectProduct( LRdim,LRdim,operamatsma2(:,operaindex2),smacolindex2(:,operaindex2),smarowindex2(:,operaindex2),&
							4,4,onesitemat(:,8),osmcolindex(:,8),osmrowindex(:,8),&
							operamatbig3(:,operaindex3),bigcolindex3(:,operaindex3),bigrowindex3(:,operaindex3),bigdim3)
				end if
			end if
		end do
!
	end if
!===================================================================

	! wait for the sendrequest
	do i=orbstart,orbend,1
		if(myid==orbid1(i,1)) then
			call MPI_WAIT(sendrequest,status,ierr)
			exit
		end if
	end do
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	if(allocated(packbuf)) deallocate(packbuf)
	if(allocated(H0mat)) deallocate(H0mat,H0colindex,H0rowindex)
	if(allocated(Hbufmat)) deallocate(Hbufmat,Hbufcolindex,Hbufrowindex)
	deallocate(phase)

	return
end subroutine System_Big

!=================================================================
!=================================================================
end module construct_System_Big
