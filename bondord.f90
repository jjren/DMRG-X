module BondOrder_mod
! this module calculate the bond order matrix;
! is the same as the one body density matrix a(i,sigma)^+ a(j,sigma)

	use module_sparse
	use variables
	use communicate
	
	implicit none

	real(kind=r8),allocatable ::  bondordmat(:,:,:,:),transDM0(:,:,:,:),&
		transDMMO(:,:,:,:)
	
	contains
!===========================================================================
!===========================================================================

subroutine init_BOmat(orbindex)
! initiate the on site niup,nidown matrix
	use onesitematrix

	implicit none

	integer :: orbindex
	! local
	integer :: operaindex2
	
	if(myid==orbid2(orbindex,orbindex,1)) then
		call ConstructOnesiteMatrix(orbindex)
		operaindex2=orbid2(orbindex,orbindex,2)
		smarowindex2(:,operaindex2*2-1:operaindex2*2)=0

		operamatsma2(1:4,2*operaindex2-1:2*operaindex2)=onesitemat(1:4,7:8)
		smacolindex2(1:4,2*operaindex2-1:2*operaindex2)=osmcolindex(1:4,7:8)
		smarowindex2(1:5,2*operaindex2-1:2*operaindex2)=osmrowindex(1:5,7:8)
	end if
return

end subroutine init_BOmat

!===========================================================================
!===========================================================================

subroutine BondOrder
! this subroutine calculate the BondOrder Matrix in the last step
	implicit none

	integer :: error
	integer :: i,j,k
	real(kind=8) :: tmp

	call master_print_message("enter BondOrder subroutine")
	
	if(myid==0) then
		allocate(bondordmat(norbs,norbs,2,nstate),stat=error)
		if(error/=0) stop
		bondordmat=0.0D0
		
		allocate(transDM0(norbs,norbs,2,nstate),stat=error)
		if(error/=0) stop
		transDM0=0.0D0
	end if
	
	call Calc_BOmat_link
	call Calc_BOmat_subspace('L')
	call Calc_BOmat_subspace('R')
	
	if(myid==0) then
		! bondorder/one partical density matrix
		call master_print_message("bondorder matrix")
		open(unit=399,file="bondord.out",status="replace")
		do k=1,nstate,1
		do i=1,norbs,1
		do j=i,norbs,1
			if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
				if(i==j) then    ! recover niup+nidown,niup-nidown to niup,nidown
					tmp=bondordmat(i,i,1,k)-bondordmat(i,i,2,k)
					bondordmat(i,i,1,k)=(bondordmat(i,i,1,k)+bondordmat(i,i,2,k))/2.0D0
					bondordmat(i,i,2,k)=tmp/2.0D0
				end if
				write(*,*) i,j,bondordmat(i,j,:,k)
				if(bondlink(i,j)==1) then
					write(399,*) i,j,bondordmat(i,j,:,k)
				end if
			end if
		end do
		end do
		end do
		close(399)
		!  transition density matrix
		call master_print_message("transition density matrix")
		do k=2,nstate,1
			do i=1,norbs,1
			do j=1,norbs,1
				if(i==j) then    ! recover niup+nidown,niup-nidown to niup,nidown
					tmp=transDM0(i,i,1,k)-transDM0(i,i,2,k)
					transDM0(i,i,1,k)=(transDM0(i,i,1,k)+transDM0(i,i,2,k))/2.0D0
					transDM0(i,i,2,k)=tmp/2.0D0
				end if
				write(*,*) i,j,transDM0(i,j,:,k)
			end do
			end do
		end do
		call transDMAO2MO
		if(nstate>1) then
			call NatTraOrb
		end if
	end if
return
end subroutine BondOrder

!===========================================================================
!===========================================================================

subroutine transDMAO2MO
! the transition density matrix or one partical density matrix from PPP-AO to MO
! MOij=CDC^+   C is MO*AO column format
	use MeanField
	use MKL95_BLAS
	use MKL95_PRECISION
	implicit none
	integer :: i,j,k,l
	real(kind=r8),allocatable :: midmat(:,:)
	real(kind=r8) tmp

	allocate(transDMMO(norbs,norbs,2,nstate))
	allocate(midmat(norbs,norbs))

	! transition density matrix
	open(unit=151,file="MO-Opdm.out",status="replace")
	do i=2,nstate,1
		do j=1,2,1   ! spin up down
			call gemm(coeffC,transDM0(:,:,j,i),midmat,'T','N',1.0D0,0.0D0)
			call gemm(midmat,coeffC,transDMMO(:,:,j,i),'N','N',1.0D0,0.0D0)
			write(151,*) j,i
			write(151,*) transDMMO(:,:,j,i)
		end do
		do k=1,norbs,1
		do l=1,norbs,1
			tmp=abs(transDMMO(k,l,1,i)+transDMMO(k,l,2,i))
			if(tmp>0.44) then
				write(151,*) k,"<<--",l,sqrt(2.0D0)/2.0D0*tmp
			end if
		end do
		end do
	end do
	close(151)

	deallocate(midmat)
return
end subroutine transDMAO2MO

!===========================================================================
!===========================================================================

subroutine NatTraOrb
! this subroutine do natural transition orbital(NTO) analysis if the ex is
! single excitation
! <ex|a^+a*ai|gs>
	
	use MKL95_LAPACK
	USE MKL95_PRECISION
	use meanfield
	implicit none
	real(kind=r8),allocatable :: Tai(:,:,:),svdvalue(:),&
	rightv(:,:),leftu(:,:),ww(:)
	integer :: info,tmp
	integer :: i

	allocate(Tai(norbs-nocc,nocc,nstate))
	
	! the operator is the singlet excitation operator 
	! sqrt(2)/2*(apup^+*aqup+apdown^+*aqdown)
	Tai(:,:,:)=sqrt(2.0D0)/2.0D0*(TransDMMO(nocc+1:norbs,1:nocc,1,:)+TransDMMO(nocc+1:norbs,1:nocc,2,:))
	
	tmp=min(nocc,norbs-nocc)
	allocate(svdvalue(tmp))
	allocate(leftu(norbs-nocc,tmp))
	allocate(rightv(tmp,nocc))
	allocate(ww(tmp-1))
	
	call master_print_message("NTO analysis result")
	open(unit=150,file="NTO.out",status="replace")
	do i=2,nstate,1
		call gesvd(Tai(:,:,i),svdvalue,leftu,rightv,ww,'N',info)
		if(info/=0) then
			write(*,*) "NatTraOrb info/=0",info
			stop
		end if
		write(*,*) "stateindex",i
		write(*,*) "svdvalue",svdvalue
		write(150,*) i
		write(150,*) svdvalue
		write(150,*) leftu
		write(150,*) rightv
	end do
	
	deallocate(Tai,leftu,rightv,ww,svdvalue)

return
end subroutine NatTraOrb

!===========================================================================
!===========================================================================

subroutine Calc_BOmat_subspace(domain)
! calculate operator in the L/R subspace i,j<=nleft+1,or i,j>=norbs-nright
! <R|<L|CLR ai^+aj CL'R'|L'>|R'>
	use exit_mod
	use mpi
	use mathlib
	implicit none
	
	character(len=1) :: domain
	
	! local
	integer :: i,j,istate,k,l
	integer :: operaindex2,orbstart,orbend
	real(kind=r8) :: ibondord(2,nstate),itransDM(2,2,nstate)
	! ibondord :: 2 means up and down ; nstate means the specific bondorder
	integer :: error,ierr
	integer :: nmid
	real(kind=r8),allocatable :: midmat(:)
	integer(kind=i4),allocatable :: midcolindex(:),midrowindex(:),coeffIFrowindexdummy(:,:)
	
	integer :: status(MPI_STATUS_SIZE) ! MPI flag

	if(myid/=0) then
		nmid=CEILING(DBLE(16*subM*subM)/pppmatratio)
		allocate(midmat(nmid),stat=error)
		if(error/=0) stop
		allocate(midcolindex(nmid),stat=error)
		if(error/=0) stop
		allocate(midrowindex(4*subM+1),stat=error)
		if(error/=0) stop
		allocate(coeffIFrowindexdummy(4*subM+1,nstate),stat=error)
		if(error/=0) stop
	end if

	! set the parameters
	if(domain=='L') then
		orbstart=1
		orbend=nleft+1
	else if(domain=='R') then
		orbstart=norbs-nright
		orbend=norbs
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
	end if
	
	! two operator matrix => no phase
	
	! for example :: CLR*CL'R'OLL'IRR'=CLR*(OLL'CL'R') the same as transition moment
	! in the R domain need to tranpose the coeffIF
	if(myid/=0) then
		if(domain=='R') then
			do istate=1,nstate,1
				call CSCtoCSR('RC',4*Rrealdim,4*Lrealdim,&
				coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindex(:,istate),&
				coeffIFrowindexdummy(:,istate))
			end do
		else
			coeffIFrowindexdummy=coeffIFrowindex
		end if
	end if

	do i=orbstart,orbend,1
	do j=i,orbend,1
		if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
			if(myid==orbid2(i,j,1)) then
				do istate=1,nstate,1
				do k=1,2,1
					operaindex2=orbid2(i,j,2)*2-2+k
					call SpMMtoSp('N','N',4*subM,4*subM,4*subM,4*subM,4*subM,&
						operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2), &
						coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
						midmat,midcolindex,midrowindex,nmid)
					! trace(CLR*QLR) or trace(CRL*QRL)
					call SpMMtrace('T',4*subM,&
						coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
						midmat,midcolindex,midrowindex,ibondord(k,istate))
					!=====================================================================================
					! calculate the transition density matrix
					! <psai1|ai^+*aj|psai2>/=<psai1|aj^+*ai|psai2>
					! <ex|ai^+*aj|gs>=<gs|aj^+*ai|ex>
					if(nstate>1) then
						if(istate==1) then
							do l=2,nstate,1
								! trace(CLR*QLR) or trace(CRL*QRL)
								call SpMMtrace('T',4*subM,&
									coeffIF(:,l),coeffIFcolindex(:,l),coeffIFrowindexdummy(:,l),&
									midmat,midcolindex,midrowindex,itransDM(1,k,l))
							end do
						else
							call SpMMtrace('T',4*subM,&
								coeffIF(:,1),coeffIFcolindex(:,1),coeffIFrowindexdummy(:,1),&
								midmat,midcolindex,midrowindex,itransDM(2,k,istate))
						end if
					end if
					!=====================================================================================
				end do
				end do
				call MPI_SEND(ibondord,nstate*2,mpi_real8,0,orbid2(i,j,2),MPI_COMM_WORLD,ierr)
				!=====================================================================================
				! calculate the transition density matrix
				! <psai1|ai^+*aj|psai2>/=<psai1|aj^+*ai|psai2>
				! <ex|ai^+*aj|gs>=<gs|aj^+*ai|ex>
				if(nstate>1) then
			!		do l=1,2,1  ! ai^+*aj   aj^+*ai
			!			if(l==1) then
			!				trans='N'
			!			else
			!				trans='T'
			!			end if
			!			do k=1,2,1
			!				operaindex2=orbid2(i,j,2)*2-2+k
			!				call SpMMtoSp(trans,'N',4*subM,4*subM,4*subM,4*subM,4*subM,&
			!					operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2), &
			!					coeffIF(:,1),coeffIFcolindex(:,1),coeffIFrowindexdummy(:,1),&
			!					midmat,midcolindex,midrowindex,nmid)
			!				do istate=2,nstate,1
			!					! trace(CLR*QLR) or trace(CRL*QRL)
			!					call SpMMtrace('T',4*subM,&
			!						coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
			!						midmat,midcolindex,midrowindex,itransDM(l,k,istate))
			!				end do
			!			end do
			!		end do
					call MPI_SEND(itransDM,nstate*4,mpi_real8,0,orbid2(i,j,2),MPI_COMM_WORLD,ierr)
				end if
				!=====================================================================================
			
			else if(myid==0) then

				call MPI_RECV(ibondord,nstate*2,mpi_real8,orbid2(i,j,1),orbid2(i,j,2),MPI_COMM_WORLD,status,ierr)

				do istate=1,nstate,1
				do k=1,2,1
					bondordmat(i,j,k,istate)=ibondord(k,istate)
					bondordmat(j,i,k,istate)=ibondord(k,istate)
				end do
				end do
				
				! transition density matrix 
				if(nstate>1) then
					call MPI_RECV(itransDM,nstate*4,mpi_real8,orbid2(i,j,1),orbid2(i,j,2),MPI_COMM_WORLD,status,ierr)
					do istate=2,nstate,1
					do k=1,2,1
						if(domain=='L') then   ! in the L space l=1 means (i,j) pair
							transDM0(i,j,k,istate)=itransDM(1,k,istate)
							transDM0(j,i,k,istate)=itransDM(2,k,istate)
						else if(domain=='R') then   ! in the R space l=1 means (j,i) pair
							transDM0(i,j,k,istate)=itransDM(2,k,istate)
							transDM0(j,i,k,istate)=itransDM(1,k,istate)
						end if
					end do
					end do
				end if
			end if
		end if
	end do
	end do

	! recovery the coeffIF
	if(myid/=0) then
		if(domain=='R') then
			do istate=1,nstate,1
				call CSCtoCSR('CR',4*Lrealdim,4*Rrealdim,&
				coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
				coeffIFrowindex(:,istate))
			end do
		end if
	end if

	if(myid/=0) deallocate(midmat,midcolindex,midrowindex,coeffIFrowindexdummy)

return

end subroutine Calc_BOmat_subspace

!===========================================================================
!===========================================================================

subroutine Calc_BOmat_link
! this subroutine calculate bond order matrix belong to different subspace
! i<=nleft+1,j>=norbs-nright 
! the same as op subrouitne transfer integral algorithm
	use mathlib
	use mpi
	use blas95
	use F95_PRECISION
	
	implicit none
	include "mkl_spblas.fi"

	integer :: hopnelement,midnelement
	real(kind=r8),allocatable :: hopmat(:,:),midmat(:),midmat2(:)
	integer(kind=i4),allocatable :: &
	hopmatcol(:,:),hopmatrow(:,:),&
	midmatcol(:),midmatrow(:),&
	midmatcol2(:),midmatrow2(:),&
	phase(:)

	character(len=1),allocatable :: hoppackbuf(:)
	integer :: position1,hoppacksize
	integer :: operaindex,nnonzero
	real(kind=r8) :: ibondord(2,nstate),itransDM(2,2,nstate)
	integer :: i,j,k,l,m,istate,p
	logical :: ifhop,ifhopsend
	integer :: hoptouched(nprocs-1),hopntouched
	integer :: error,info
	
	integer :: status(MPI_STATUS_SIZE),hopsendrequest(nprocs-1)
	integer :: ierr

	hopnelement=CEILING(DBLE(16*subM*subM)/bigratio1)
	midnelement=CEILING(DBLE(16*subM*subM)/hopmatratio)
	
	do i=1,nleft+1,1
	do j=norbs,norbs-nright,-1
		if(bondlink(i,j)==1 .or. logic_bondorder==2) then
			if(myid==orbid1(i,1) .or. myid==orbid1(j,1)) then
				if(.not. allocated(hopmat)) then
					allocate(hopmat(hopnelement,2),stat=error) ! store the hopping matrix
					if(error/=0) stop
					allocate(hopmatcol(hopnelement,2),stat=error) 
					if(error/=0) stop
					allocate(hopmatrow(4*subM+1,2),stat=error) 
					if(error/=0) stop

					hoppacksize=(hopnelement*12+4*(4*subM+1))*2+1000  ! 1000 is redundant
					allocate(hoppackbuf(hoppacksize),stat=error) ! packbuf to send hopping matrix
					if(error/=0) stop
				end if
			end if

			if(myid==orbid1(i,1)) then
				if(.not. allocated(midmat)) then
					allocate(midmat(midnelement),stat=error) ! store the intermediate matrix
					if(error/=0) stop
					allocate(midmatcol(midnelement),stat=error) 
					if(error/=0) stop
					allocate(midmatrow(4*subM+1),stat=error) 
					if(error/=0) stop
					allocate(midmat2(midnelement),stat=error) ! store the intermediate matrix
					if(error/=0) stop
					allocate(midmatcol2(midnelement),stat=error) 
					if(error/=0) stop
					allocate(midmatrow2(4*subM+1),stat=error) 
					if(error/=0) stop
				end if
			end if
		end if
	end do
	end do
	
	do i=norbs,norbs-nright,-1
		if(myid==orbid1(i,1)) then
			! check if need hopping matrix
			ifhop=.false.
			do j=1,nleft+1,1
				if(bondlink(i,j)==1 .or. logic_bondorder==2) then
					ifhop=.true.
					exit
				end if
			end do
			
			if(ifhop==.true.) then
				operaindex=orbid1(i,2)
				
				! copy operamatbig to hopmat
				! send operamatbig, not R*C
				do l=1,2,1
					! integer can not call copy
					hopmatrow(:,l)=bigrowindex1(:,operaindex*3-3+l)
					nnonzero=bigrowindex1(4*subM+1,operaindex*3-3+l)-1
					hopmatcol(1:nnonzero,l)=bigcolindex1(1:nnonzero,operaindex*3-3+l)
					call copy(operamatbig1(1:nnonzero,operaindex*3-3+l),hopmat(1:nnonzero,l))
				end do

				! pack hopmatrix
				position1=0
				do l=1,2,1
					call MPI_PACK(hopmatrow(1,l),(4*subM+1),MPI_integer4,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
					call MPI_PACK(hopmat(1,l),(hopmatrow(4*subM+1,l)-1),MPI_real8,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
					call MPI_PACK(hopmatcol(1,l),(hopmatrow(4*subM+1,l)-1),MPI_integer4,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
				end do

				! send the hopmat
				hoptouched=0
				hopntouched=0
				do l=1,nleft+1,1
					if((bondlink(i,l)==1 .or. logic_bondorder==2) .and. orbid1(l,1)/=myid) then ! if orbid(l)==myid need not send hopmat
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

		! recv hopmat
		! hopping term calculation
		if(myid/=orbid1(i,1) .and. myid/=0) then
			do l=1,nleft+1,1
				if((bondlink(l,i)==1 .or. logic_bondorder==2) .and. myid==orbid1(l,1)) then  
					call MPI_RECV(hoppackbuf,hoppacksize,MPI_PACKED,orbid1(i,1),i,MPI_COMM_WORLD,status,ierr)
					position1=0
					do j=1,2,1
						call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmatrow(1,j),(4*subM+1),MPI_integer4,MPI_COMM_WORLD,ierr)
						call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmat(1,j),(hopmatrow(4*subM+1,j)-1),MPI_real8,MPI_COMM_WORLD,ierr)
						call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmatcol(1,j),(hopmatrow(4*subM+1,j)-1),MPI_integer4,MPI_COMM_WORLD,ierr)
					end do
					exit  ! only recv once
				end if
			end do
		end if
		
		do l=1,nleft+1,1
			if(bondlink(i,l)==1 .or. logic_bondorder==2) then
			if(myid==orbid1(l,1)) then
				operaindex=orbid1(l,2)
				
				!-----------------------------------------------------
				! the +1 -1 phase added to l' 
				allocate(phase(4*subM),stat=error)
				if(error/=0) stop

				do j=1,4*subM,1
					phase(j)=(-1.0D0)**(mod(quantabigL(j,1),2))
				end do
				!----------------------------------------------------

				! construct hopmat
				do j=1,nstate,1
					do k=1,2,1
						! firstly needed to transfer from a(+) to a, so 'T' and 'T' is 'N'
						! k=1 a up,k=2 a down
						call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
								coeffIF(:,j),coeffIFcolindex(:,j),coeffIFrowindex(:,j), &
								hopmat(:,k),hopmatcol(:,k),hopmatrow(:,k), &
								midmat,midmatcol,midmatrow,midnelement,info)
						call checkinfo(info)

						! add phase
						do p=1,4*subM,1
						do m=midmatrow(p),midmatrow(p+1)-1,1
							midmat(m)=midmat(m)*phase(p)
						end do
						end do

						!k<=2 al^+*ar
						call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
								operamatbig1(:,operaindex*3-3+k),bigcolindex1(:,operaindex*3-3+k),bigrowindex1(:,operaindex*3-3+k), &
								midmat,midmatcol,midmatrow,&
								midmat2,midmatcol2,midmatrow2,midnelement,info)
						call checkinfo(info)

						! trace(CLR*OLR)
						call SpMMtrace('T',4*subM, & 
								coeffIF(:,j),coeffIFcolindex(:,j),coeffIFrowindex(:,j), &
								midmat2,midmatcol2,midmatrow2,ibondord(k,j))

						!=========================================================================
						! calculate transition density matrix
						! <ex|aR^+*aL|gs>=<gs|aL^+*aR|ex>
						if(nstate>1) then
							if(j/=1) then
								! trace(CLR*OLR)
								call SpMMtrace('T',4*subM, & 
										coeffIF(:,1),coeffIFcolindex(:,1),coeffIFrowindex(:,1), &
										midmat2,midmatcol2,midmatrow2,itransDM(2,k,j))
							else if(j==1) then
								! trace(CLR*OLR)
								do istate=2,nstate,1
									call SpMMtrace('T',4*subM, & 
											coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindex(:,istate), &
											midmat2,midmatcol2,midmatrow2,itransDM(1,k,istate))
								end do
							end if
						end if
						!=========================================================================
					end do
				end do
				call MPI_SEND(ibondord,nstate*2,mpi_real8,0,1,MPI_COMM_WORLD,ierr)
				if(nstate>1) then
					call MPI_SEND(itransDM,nstate*4,mpi_real8,0,1,MPI_COMM_WORLD,ierr)
				end if
				!=========================================================================

				deallocate(phase)
			else if(myid==0) then
				call MPI_RECV(ibondord,nstate*2,mpi_real8,orbid1(l,1),1,MPI_COMM_WORLD,status,ierr)
				do istate=1,nstate,1
				do k=1,2,1
					bondordmat(i,l,k,istate)=ibondord(k,istate)
					bondordmat(l,i,k,istate)=ibondord(k,istate)
				end do
				end do
				if(nstate>1) then
					call MPI_RECV(itransDM,nstate*4,mpi_real8,orbid1(l,1),1,MPI_COMM_WORLD,status,ierr)
					do istate=2,nstate,1
					do k=1,2,1
						transDM0(i,l,k,istate)=itransDM(2,k,istate)
						transDM0(l,i,k,istate)=itransDM(1,k,istate)
					end do
					end do
				end if
			end if
			end if
		end do

		! confirm that the pppVmat and hopmat can be used again without problem
		if(myid==orbid1(i,1) .and. ifhop==.true.) then
			do j=1,hopntouched,1
				call MPI_WAIT(hopsendrequest(j),status,ierr)
			end do
		end if
	end do

	if(allocated(hopmat)) deallocate(hopmat,hopmatcol,hopmatrow)
	if(allocated(midmat)) deallocate(midmat,midmatcol,midmatrow)
	if(allocated(midmat2)) deallocate(midmat2,midmatcol2,midmatrow2)
	if(allocated(hoppackbuf)) deallocate(hoppackbuf)

return
end subroutine Calc_BOmat_link

!===========================================================================
!===========================================================================
end module BondOrder_mod
