subroutine op(bigdim,smadim,coeff,newcoeff)
! this is the core subroutine to calculate the S*H*S*C or H*C
! the parallel schema follow JCP 12 3174(2004) garnet chan
! if want to save memory, then can write a wrapper, to send one coeff every time

!--------------------------------------------------------
! input bigdim,smadim,coeff
! bigdim is the totaldim 16M*M
! bigdim may be < 16M*M because we use good quantum number and spin symmetry,
! and the dimension  may be half
! if groud state smadim=1
! if gs+ex smadim may be >1
!--------------------------------------------------------
! output newcoeff
! coeff is the input coefficient and in the 1-d arrary format
! new coeff is H cross C result
!---------------------------------------------------------

	use mpi
	use variables
	use symmetry
	use mathlib
	use module_sparse
	use BLAS95
	use F95_PRECISION

	implicit none
	include "mkl_spblas.fi"

	integer :: bigdim,smadim
	real(kind=r8) :: coeff(bigdim*smadim),newcoeff(bigdim*smadim)
	
	! local
	integer :: operaindex
	
	real(kind=r8),allocatable :: LRcoeffin(:,:),LRcoeffout(:,:),coeffnosymm(:),coeffnosymmreduce(:)
	integer(kind=i4),allocatable :: LRcoeffincol(:,:),LRcoeffinrow(:,:),&
	LRcoeffoutcol(:,:),LRcoeffoutrow(:,:),&
	LRcoeffinrowdummy(:,:)
	
	real(kind=r8),allocatable :: hopmat(:,:,:),pppVmat(:,:),buffmat(:)
	integer(kind=i4),allocatable :: &
	hopmatcol(:,:,:),pppVmatcol(:,:),buffmatcol(:),&
	hopmatrow(:,:,:),pppVmatrow(:,:),buffmatrow(:)

	real(kind=r8),allocatable :: operamatbuf(:)
	integer(kind=i4),allocatable :: operamatbufcol(:),operamatbufrow(:)

	integer(kind=i4) :: pppnelement,hopnelement,LRoutnelement

	character(len=1),allocatable :: hoppackbuf(:),pppVpackbuf(:)
	integer :: position1,pppVpacksize,hoppacksize

	real(kind=r8),allocatable :: phase(:)
	integer :: error,i,j,k,l,m,n
	integer :: info
	
	logical :: ifhop,ifhopsend,ifpppVsend
	integer :: hoptouched(nprocs-1),hopntouched,pppVtouched(nprocs-1),pppVntouched
	
	! MPI flag
	integer :: status(MPI_STATUS_SIZE),hopsendrequest(nprocs-1),hoprecvrequest
	integer :: ierr
	
!============================================================
! allocate workspace
	! store nosymmetry coeff
	allocate(coeffnosymm(ngoodstates*smadim),stat=error)
	if(error/=0) stop

	! set the sparse matrix dim
	pppnelement=NINT(DBLE(16*Lrealdim*Rrealdim)/pppmatratio,i4)
	hopnelement=NINT(DBLE(16*Lrealdim*Rrealdim)/hopmatratio,i4)
	LRoutnelement=NINT(DBLE(16*Lrealdim*Rrealdim)/LRoutratio,i4)
	
	if(myid/=0) then
		do i=norbs,norbs-nright,-1
			if(myid==orbid1(i,1)) then
				if( .not. allocated(LRcoeffin)) then
					! transform the 1-array to 4M*4M form 
					allocate(LRcoeffin(ngoodstates,smadim),stat=error)   ! coeff to LR format
					if(error/=0) stop
					allocate(LRcoeffincol(ngoodstates,smadim),stat=error)   
					if(error/=0) stop
					allocate(LRcoeffinrow(4*Lrealdim+1,smadim),stat=error)   
					if(error/=0) stop
					exit
				end if
			end if
		end do
	
		do i=1,nleft+1,1
			if(myid==orbid1(i,1)) then
				if( .not. allocated(LRcoeffout)) then
					allocate(LRcoeffout(LRoutnelement,smadim),stat=error)  ! newcoeff to LR format
					if(error/=0) stop
					allocate(LRcoeffoutcol(LRoutnelement,smadim),stat=error)  
					if(error/=0) stop
					allocate(LRcoeffoutrow(4*Lrealdim+1,smadim),stat=error) 
					if(error/=0) stop
					! intermediate array store LRcoeffout
					allocate(buffmat(LRoutnelement),stat=error)
					if(error/=0) stop
					allocate(buffmatcol(LRoutnelement),stat=error)
					if(error/=0) stop
					allocate(buffmatrow(4*Lrealdim+1),stat=error)
					if(error/=0) stop
					exit
				end if
			end if
		end do

		do i=1,nleft+1,1
		do j=norbs,norbs-nright,-1
			if(bondlink(i,j)==1) then
				if(myid==orbid1(i,1) .or. myid==orbid1(j,1)) then
					if(.not. allocated(hopmat)) then
						allocate(hopmat(hopnelement,4,smadim),stat=error) ! store the hopping matrix
						if(error/=0) stop
						allocate(hopmatcol(hopnelement,4,smadim),stat=error) 
						if(error/=0) stop
						allocate(hopmatrow(4*Lrealdim+1,4,smadim),stat=error) 
						if(error/=0) stop

						hoppacksize=(hopnelement*12+4*(4*Lrealdim+1))*smadim*4+1000  ! 1000 is redundant
						allocate(hoppackbuf(hoppacksize),stat=error) ! packbuf to send hopping matrix
						if(error/=0) stop
					end if
				end if
			end if
		end do
		end do

		allocate(pppVmat(pppnelement,smadim),stat=error) ! store the pppV matrix
		if(error/=0) stop
		allocate(pppVmatcol(pppnelement,smadim),stat=error) 
		if(error/=0) stop
		allocate(pppVmatrow(4*Lrealdim+1,smadim),stat=error) 
		if(error/=0) stop

		pppVpacksize=(pppnelement*12+4*(4*Lrealdim+1))*smadim+1000
		allocate(pppVpackbuf(pppVpacksize),stat=error) ! packbuf to send the pppV matrix
		if(error/=0) stop
	else  
		! 0 process
		allocate(LRcoeffin(ngoodstates,smadim),stat=error)   ! coeff to LR format
		if(error/=0) stop
		allocate(LRcoeffincol(ngoodstates,smadim),stat=error)   
		if(error/=0) stop
		allocate(LRcoeffinrow(4*Lrealdim+1,smadim),stat=error)   
		if(error/=0) stop
		allocate(LRcoeffout(LRoutnelement,smadim),stat=error)  ! newcoeff to LR format
		if(error/=0) stop
		allocate(LRcoeffoutcol(LRoutnelement,smadim),stat=error)  
		if(error/=0) stop
		allocate(LRcoeffoutrow(4*Lrealdim+1,smadim),stat=error) 
		if(error/=0) stop
		! intermediate array store LRcoeffout
		allocate(buffmat(LRoutnelement),stat=error)
		if(error/=0) stop
		allocate(buffmatcol(LRoutnelement),stat=error)
		if(error/=0) stop
		allocate(buffmatrow(4*Lrealdim+1),stat=error)
		if(error/=0) stop
	end if

!=================================================================================================
	
	! unsymmetrize the coeff if needed
	if( myid==0 ) then
		! if symmetry==.true. then transform the symmetry coeff to the unsymmetry coeff
		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
			if(bigdim/=nsymmstate) then
				call master_print_message(bigdim,"In symmetry, op bigdim/=nsymmstate wrong!")
				stop
			end if
			do i=1,smadim,1
				call symmetrizestate(ngoodstates,coeffnosymm(ngoodstates*(i-1)+1:i*ngoodstates),&
					coeff(bigdim*(i-1)+1:i*bigdim),'u')
			end do
		else
			if(bigdim/=ngoodstates) then
				call master_print_message(bigdim,"op bigdim/=ngoodstates wrong!")
				stop
			end if
			coeffnosymm=coeff
		end if
	end if

	! send coeffnosymm to R space process
	do i=1,nprocs-1,1
		if(myid==0 .or. myid==i) then
			do j=norbs,norbs-nright,-1
				if(i==orbid1(j,1)) then
					if(myid==0) then
						call MPI_SEND(coeffnosymm,ngoodstates*smadim,MPI_real8,i,1,MPI_COMM_WORLD,ierr)
					else if(myid==i) then
						call MPI_RECV(coeffnosymm,ngoodstates*smadim,MPI_real8,0,1,MPI_COMM_WORLD,status,ierr)
					end if
					exit
				end if
			end do
		end if
	end do

!-----------------------------------------------------------------------------

! R space process do it
! to transform the 16M*M coeff to 4M*4M(L*R) format ; coeff(16M^2,n) to coeff(4M,4M,n) 
! since the input coeff is ngoodstates and other nongoodstates sets to 0
	if(allocated(LRcoeffin)) then
		! CSR format and fortran's column major format not corresponds
		allocate(LRcoeffinrowdummy(4*Rrealdim+1,smadim),stat=error)
		if(error/=0) stop
		LRcoeffinrowdummy(1,:)=1

		! in the CSC form
		n=0
		do k=1,smadim,1
			m=0
			do i=1,4*Rrealdim,1
			do j=1,4*Lrealdim,1
				if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
					quantabigL(j,2)+quantabigR(i,2)==totalSz) then
					m=m+1
					LRcoeffin(m,k)=coeffnosymm(m+n)
					LRcoeffincol(m,k)=j
				end if
			end do
			LRcoeffinrowdummy(i+1,k)=m+1
			end do
			n=m+n
		end do

		! CSC transfer to CSR form
		do k=1,smadim,1
			call CSCtoCSR('CR',4*Lrealdim,4*Rrealdim,LRcoeffin(:,k),LRcoeffincol(:,k),LRcoeffinrowdummy(:,k),LRcoeffinrow(:,k))
		end do

		deallocate(LRcoeffinrowdummy)
	end if

! L space process and 0 process initializaiton
! L space process have LRcoeffout and send to 0 process at last

	if(allocated(LRcoeffout)) then
		LRcoeffoutrow=1  ! define the LRcoeffout matrix is 0
	end if

!  calculate HL*1 and 1*HR 
	if(myid==0) then 
		do i=1,smadim,1
			do k=1,2,1
				if(k==1) then ! HL*1
					call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Lrealdim,4*Rrealdim, &
						Hbig(:,1),Hbigcolindex(:,1),Hbigrowindex(:,1), &
						LRcoeffin(:,i),LRcoeffincol(:,i),LRcoeffinrow(:,i), &
						buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
					call checkinfo(info)
				else ! 1*HR
					call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Rrealdim,4*Rrealdim, &
						LRcoeffin(:,i),LRcoeffincol(:,i),LRcoeffinrow(:,i), &
						Hbig(:,2),Hbigcolindex(:,2),Hbigrowindex(:,2), &
						buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
					call checkinfo(info)
				end if
				! add LRcoeffout and bufmat
				call SpMatAdd(4*Rrealdim,4*Lrealdim,LRcoeffout(:,i),LRcoeffoutcol(:,i),LRcoeffoutrow(:,i),&
				'N',1.0D0,4*Rrealdim,4*Lrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
			end do
		end do
	end if

!------------------------------------------------
! vlr=Hlrl'r'*Cl'r'=sum(opt,l',r')=sum(Lopt,l') parity*Oll'*sum(Ropt,r') Orr'cl'r'
! the parallel schema is that 0 process bcast the coeff matrix to other process
! and 0 process gather the result
	if(myid/=0) then
	do i=norbs,norbs-nright,-1
		if(myid==orbid1(i,1)) then
			operaindex=orbid1(i,2)
!=====================================================================================
			
			! construct the pppVmat
			do j=1,smadim,1
				call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Rrealdim,4*Rrealdim, &
					LRcoeffin(:,j),LRcoeffincol(:,j),LRcoeffinrow(:,j), &
					operamatbig1(:,operaindex*3),bigcolindex1(:,operaindex*3),bigrowindex1(:,operaindex*3), &
					pppVmat(:,j),pppVmatcol(:,j),pppVmatrow(:,j),pppnelement,info)
				call checkinfo(info)
			end do
			! pack the pppVmat
			position1=0
			do k=1,smadim,1
				call MPI_PACK(pppVmatrow(1,k),(4*Lrealdim+1),MPI_integer4,pppVpackbuf,pppVpacksize,position1,MPI_COMM_WORLD,ierr)
				call MPI_PACK(pppVmat(1,k),pppVmatrow(4*Lrealdim+1,k)-1,MPI_real8,pppVpackbuf,pppVpacksize,position1,MPI_COMM_WORLD,ierr)
				call MPI_PACK(pppVmatcol(1,k),pppVmatrow(4*Lrealdim+1,k)-1,MPI_integer4,pppVpackbuf,pppVpacksize,position1,MPI_COMM_WORLD,ierr)
			end do
			
			! send the pppVmat
			pppVtouched=0
			pppVntouched=0
			do l=1,nleft+1,1
				if(orbid1(l,1)/=myid) then ! if l==myid need not send hopmat
					ifpppVsend=.false.
					do m=1,pppVntouched,1
						if(orbid1(l,1)==pppVtouched(m)) then
							ifpppVsend=.true.   ! have send pppVmat to this process
							exit
						end if
					end do
					if(ifpppVsend==.false.) then
						pppVntouched=pppVntouched+1
						pppVtouched(pppVntouched)=orbid1(l,1)
						call MPI_SEND(pppVpackbuf,position1,MPI_PACKED,orbid1(l,1),i,MPI_COMM_WORLD,ierr)
					!	call MPI_ISEND(pppVmat,16*Lrealdim*Rrealdim*smadim,MPI_real8,orbid(l),i,MPI_COMM_WORLD,pppVsendrequest(pppVntouched),ierr)
					! some problem in the isend, maybe the system buffer size limitation
					! use bsend is possible
					end if
				end if
			end do
!====================================================================
			! check if need hopping matrix
			ifhop=.false.
			do j=1,nleft+1,1
				if(bondlink(i,j)==1) then
					ifhop=.true.
					exit
				end if
			end do
			
			if(ifhop==.true.) then
				allocate(operamatbuf(bigdim1),stat=error)
				if(error/=0) stop
				allocate(operamatbufcol(bigdim1),stat=error)
				if(error/=0) stop
				allocate(operamatbufrow(4*Rrealdim+1),stat=error)
				if(error/=0) stop

				! construct hopmat
				do j=1,smadim,1
					do k=1,4,1
						if(k<=2) then
						! firstly needed to transfer from a(+) to a, so 'T' and 'T' is 'N', (ni-1)^+=(ni-1)
						! k=1 a up,k=2 a down,k=3 n,k=4 a+ up,k=5 a+ down;
							call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Rrealdim,4*Rrealdim, &
									LRcoeffin(:,j),LRcoeffincol(:,j),LRcoeffinrow(:,j), &
									operamatbig1(:,(operaindex-1)*3+k),bigcolindex1(:,operaindex*3-3+k),bigrowindex1(:,operaindex*3-3+k), &
									hopmat(:,k,j),hopmatcol(:,k,j),hopmatrow(:,k,j),hopnelement,info)
							call checkinfo(info)
						else
							!call gemm(LRcoeffin(:,:,j),operamatbig(1:4*Rrealdim,1:4*Rrealdim,(operaindex-1)*3+k-2)&
							!,hopmat(:,:,k,j),'N','T',1.0D0,0.0D0)
							operamatbuf=operamatbig1(:,operaindex*3-5+k)
							operamatbufcol=bigcolindex1(:,operaindex*3-5+k)
							! CSR to CSC need transfer N to T ; the operamatbig here needs T
							call CSCtoCSR('RC',4*Rrealdim,4*Rrealdim, &
							operamatbuf,operamatbufcol,bigrowindex1(:,operaindex*3-5+k),operamatbufrow)

							call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Rrealdim,4*Rrealdim, &
									LRcoeffin(:,j),LRcoeffincol(:,j),LRcoeffinrow(:,j), &
									operamatbuf,operamatbufcol,operamatbufrow, &
									hopmat(:,k,j),hopmatcol(:,k,j),hopmatrow(:,k,j),hopnelement,info)
							call checkinfo(info)
						end if
					end do
				end do

				deallocate(operamatbuf,operamatbufcol,operamatbufrow)
			!---------------------------------------------------------
				! the +1 -1 phase added to l' of hopmat
				allocate(phase(4*Lrealdim),stat=error)
				if(error/=0) stop

				do j=1,4*Lrealdim,1
					phase(j)=(-1.0D0)**(mod(quantabigL(j,1),2))
				end do

				do m=1,smadim,1
				do l=1,2,1
				do j=1,4*Lrealdim,1
				do k=hopmatrow(j,l,m),hopmatrow(j+1,l,m)-1,1
					hopmat(k,l,m)=hopmat(k,l,m)*phase(j)
				end do
				end do
				end do
				end do

				phase=phase*(-1.0D0)
				do m=1,smadim,1
				do l=3,4,1
				do j=1,4*Lrealdim,1
				do k=hopmatrow(j,l,m),hopmatrow(j+1,l,m)-1,1
					!transfer from al*ar^(+) to ar^(+)*al
					hopmat(k,l,m)=hopmat(k,l,m)*phase(j)
				end do
				end do
				end do
				end do

				deallocate(phase)
			!-----------------------------------------------------------
				! pack hopmatrix
				position1=0
				do k=1,smadim,1
					do l=1,4,1
						call MPI_PACK(hopmatrow(1,l,k),(4*Lrealdim+1),MPI_integer4,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
						call MPI_PACK(hopmat(1,l,k),(hopmatrow(4*Lrealdim+1,l,k)-1),MPI_real8,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
						call MPI_PACK(hopmatcol(1,l,k),(hopmatrow(4*Lrealdim+1,l,k)-1),MPI_integer4,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
					end do
				end do

				! send the hopmat
				hoptouched=0
				hopntouched=0
				do l=1,nleft+1,1
					if(bondlink(i,l)==1 .and. orbid1(l,1)/=myid) then ! if orbid(l)==myid need not send hopmat
						ifhopsend=.false.
						do m=1,hopntouched,1
							if(orbid1(l,1)==hoptouched(m)) then
								ifhopsend=.true.   ! have send hopmat to this process
								exit
							end if
						end do
						if(ifhopsend==.false.) then
							hopntouched=hopntouched+1
							hoptouched(hopntouched)=orbid1(l,1)
							call MPI_ISEND(hoppackbuf,position1,MPI_PACKED,orbid1(l,1),i,MPI_COMM_WORLD,hopsendrequest(hopntouched),ierr)
						end if
					end if
				end do
			end if
		end if
!===============================================================================
	
! L space recv pppVmat hopmat---------------------------------------
		if(myid/=orbid1(i,1)) then
			! pppVmat recv
			do l=1,nleft+1,1
				if(myid==orbid1(l,1)) then
					call MPI_RECV(pppVpackbuf,pppVpacksize,MPI_PACKED,orbid1(i,1),i,MPI_COMM_WORLD,status,ierr)
					position1=0
					do k=1,smadim,1
						call MPI_UNPACK(pppVpackbuf,pppVpacksize,position1,pppVmatrow(1,k),(4*Lrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
						call MPI_UNPACK(pppVpackbuf,pppVpacksize,position1,pppVmat(1,k),pppVmatrow(4*Lrealdim+1,k)-1,MPI_real8,MPI_COMM_WORLD,ierr)
						call MPI_UNPACK(pppVpackbuf,pppVpacksize,position1,pppVmatcol(1,k),pppVmatrow(4*Lrealdim+1,k)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
					end do
					exit  ! only recv once
				end if
			end do
			! hopmat recv
			do l=1,nleft+1,1
				if(bondlink(l,i)==1 .and. myid==orbid1(l,1)) then  
					call MPI_IRECV(hoppackbuf,hoppacksize,MPI_PACKED,orbid1(i,1),i,MPI_COMM_WORLD,hoprecvrequest,ierr)
					exit  ! only recv once
				end if
			end do
		end if
!---------------------------------------------------------------------
! pppV calculation
		do l=1,nleft+1,1
			if(myid==orbid1(l,1)) then
				operaindex=orbid1(l,2)

				do j=1,smadim,1
					! buffmat is to save the intermediate matrix
					call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Lrealdim,4*Rrealdim, &
						operamatbig1(:,operaindex*3),bigcolindex1(:,operaindex*3),bigrowindex1(:,operaindex*3), &
						pppVmat(:,j),pppVmatcol(:,j),pppVmatrow(:,j), &
						buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
					call checkinfo(info)
					! add LRcoeffout and buffmat
					call SpMatAdd(4*Rrealdim,4*Lrealdim,LRcoeffout(:,j),LRcoeffoutcol(:,j),LRcoeffoutrow(:,j),&
					'N',pppV(i,l),4*Rrealdim,4*Lrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
				end do
			end if
		end do
!---------------------------------------------------------------------
! hopping term calculation
		! wait to recv hopmat
		if(myid/=orbid1(i,1)) then
			do l=1,nleft+1,1
				if(bondlink(l,i)==1 .and. myid==orbid1(l,1)) then  
					call MPI_WAIT(hoprecvrequest,status,ierr)
					position1=0
					do k=1,smadim,1
						do j=1,4,1
							call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmatrow(1,j,k),(4*Lrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
							call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmat(1,j,k),(hopmatrow(4*Lrealdim+1,j,k)-1),MPI_real8,MPI_COMM_WORLD,ierr)
							call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmatcol(1,j,k),(hopmatrow(4*Lrealdim+1,j,k)-1),MPI_integer4,MPI_COMM_WORLD,ierr)
						end do
					end do
					exit  ! every process only wait once
				end if
			end do
		end if

		do l=1,nleft+1,1
			if(bondlink(i,l)==1 .and. myid==orbid1(l,1)) then
				operaindex=orbid1(l,2)

				do j=1,smadim,1
					do k=1,4,1
						!k<=2 al^+*ar,k>2 al*ar^(+) 
						if(k<=2) then
							call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Lrealdim,4*Rrealdim, &
								operamatbig1(:,operaindex*3-3+k),bigcolindex1(:,operaindex*3-3+k),bigrowindex1(:,operaindex*3-3+k), &
								hopmat(:,k,j),hopmatcol(:,k,j),hopmatrow(:,k,j), &
								buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
							call checkinfo(info)
						else
							call mkl_dcsrmultcsr('T',0,8,4*Lrealdim,4*Lrealdim,4*Rrealdim, &
								operamatbig1(:,operaindex*3-5+k),bigcolindex1(:,operaindex*3-5+k),bigrowindex1(:,operaindex*3-5+k), &
								hopmat(:,k,j),hopmatcol(:,k,j),hopmatrow(:,k,j), &
								buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
							call checkinfo(info)
						end if
						call SpMatAdd(4*Rrealdim,4*Lrealdim,LRcoeffout(:,j),LRcoeffoutcol(:,j),LRcoeffoutrow(:,j),&
						'N',t(i,l),4*Rrealdim,4*Lrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
					end do
				end do
			end if
		end do
		
		! confirm that the pppVmat and hopmat can be used again without problem
		if(myid==orbid1(i,1) .and. ifhop==.true.) then
		!	do j=1,pppVntouched,1
		!		call MPI_WAIT(pppVsendrequest(j),status,ierr)
		!	end do
			do j=1,hopntouched,1
				call MPI_WAIT(hopsendrequest(j),status,ierr)
			end do
		end if
	end do
	end if   ! myid/=0 ends

	! every process transfer LRcoeffout to coeffnosymm
	if(allocated(LRcoeffout)) then
		m=0
		do k=1,smadim,1
		do i=1,4*Rrealdim,1
		do j=1,4*Lrealdim,1
			if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
				quantabigL(j,2)+quantabigR(i,2)==totalSz) then
				m=m+1
				call SpMatIJ(4*Lrealdim,j,i,LRcoeffout(:,k),LRcoeffoutcol(:,k),LRcoeffoutrow(:,k),coeffnosymm(m))
			end if
		end do
		end do
		end do
		if(m/=ngoodstates*smadim) then
			write(*,*) "========================"
			write(*,*) "m/=ngoodstates*k failed!",m,smadim
			write(*,*) "========================"
			stop
		end if
	else
		coeffnosymm=0.0D0   ! other process coeffnosymm does not sum up
	end if
	
	if(allocated(LRcoeffin)) deallocate(LRcoeffin,LRcoeffinrow,LRcoeffincol)
	if(allocated(LRcoeffout)) deallocate(LRcoeffout,LRcoeffoutrow,LRcoeffoutcol)
	if(allocated(buffmat)) deallocate(buffmat,buffmatcol,buffmatrow)
	if(allocated(pppVmat)) deallocate(pppVmat,pppVmatrow,pppVmatcol)
	if(allocated(hopmat)) deallocate(hopmat,hopmatrow,hopmatcol)
	if(allocated(hoppackbuf)) deallocate(hoppackbuf)
	if(allocated(pppVpackbuf)) deallocate(pppVpackbuf)
	
	if(myid==0) then
		allocate(coeffnosymmreduce(ngoodstates*smadim),stat=error)
		if(error/=0) stop
	end if

	call MPI_REDUCE(coeffnosymm,coeffnosymmreduce,ngoodstates*smadim,mpi_real8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	
	if(myid==0) then
		newcoeff=0.0D0
		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
			do i=1,smadim,1
				call symmetrizestate(ngoodstates,coeffnosymmreduce(ngoodstates*(i-1)+1:i*ngoodstates),&
					newcoeff(bigdim*(i-1)+1:i*bigdim),'s')
			end do
		else
			newcoeff=coeffnosymmreduce
		end if
	end if

	if(allocated(coeffnosymm)) deallocate(coeffnosymm)
	if(allocated(coeffnosymmreduce)) deallocate(coeffnosymmreduce)
	
return

end subroutine op




