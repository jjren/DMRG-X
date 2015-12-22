module masterdiag
! master process gather all the operamatbig1 to do diagonalize
	use kinds_mod
	use variables
	use communicate
	use module_sparse
	use symmetry
	use mpi
	implicit none
	private

	character(len=1),allocatable :: packbuf(:)
	integer :: packsize
	public :: masterop,MasterGather
contains

!================================================================
!================================================================

subroutine masterop(bigdim,smadim,coeff,newcoeff)
	use mathlib
	implicit none
	include "mkl_spblas.fi"
	integer :: bigdim,smadim
	real(kind=r8) :: coeff(bigdim*smadim),newcoeff(bigdim*smadim)
	
	! local
	real(kind=r8),allocatable :: LRcoeffin(:,:),LRcoeffout(:,:),coeffnosymm(:)
	integer(kind=i4),allocatable :: LRcoeffincol(:,:),LRcoeffinrow(:,:),&
	LRcoeffoutcol(:,:),LRcoeffoutrow(:,:)
	real(kind=r8),allocatable :: hopmat(:),pppVmat(:),buffmat(:)
	integer(kind=i4),allocatable :: &
	hopmatcol(:),pppVmatcol(:),buffmatcol(:),&
	hopmatrow(:),pppVmatrow(:),buffmatrow(:)
	
	integer(kind=i4) :: pppnelement,hopnelement,LRoutnelement

	real(kind=r8),allocatable :: phase(:,:)
	integer :: error,i,j,k,l,m,l1,l2
	integer :: info
	logical :: ifhop
	
! allocate workspace
	! store nosymmetry coeff
	allocate(coeffnosymm(ngoodstates*smadim))

	! set the sparse matrix dim
	pppnelement=CEILING(DBLE(16*subM*subM)/pppmatratio)
	hopnelement=CEILING(DBLE(16*subM*subM)/hopmatratio)
	LRoutnelement=CEILING(DBLE(16*subM*subM)/LRoutratio)

	allocate(LRcoeffin(ngoodstates,smadim))   ! coeff to LR format
	allocate(LRcoeffincol(ngoodstates,smadim))   
	allocate(LRcoeffinrow(4*Lrealdim+1,smadim))   
	allocate(LRcoeffout(LRoutnelement,smadim))  ! newcoeff to LR format
	allocate(LRcoeffoutcol(LRoutnelement,smadim))  
	allocate(LRcoeffoutrow(4*Lrealdim+1,smadim)) 
	allocate(buffmat(LRoutnelement))
	allocate(buffmatcol(LRoutnelement))
	allocate(buffmatrow(4*Lrealdim+1))
	allocate(pppVmat(pppnelement)) ! store the pppV matrix
	allocate(pppVmatcol(pppnelement)) 
	allocate(pppVmatrow(4*Lrealdim+1)) 

	! unsymmetrize the coeff if needed
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

	! to transform the 16M*M coeff to 4M*4M(L*R) format ; coeff(16M^2,n) to coeff(4M,4M,n) 
	do k=1,smadim,1
		call coefftosparse(ngoodstates,&
			LRcoeffin(:,k),LRcoeffincol(:,k),LRcoeffinrow(:,k),&
			ngoodstates,coeffnosymm((k-1)*ngoodstates+1:k*ngoodstates),&
			Lrealdim,Rrealdim,quantabigL(1:4*Lrealdim,1:2),quantabigR(1:4*Rrealdim,1:2))
	end do

	LRcoeffoutrow=1  ! define the LRcoeffout matrix is 0
	
	!  calculate HL*1 and 1*HR 
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

	! PPP term
	do i=norbs,norbs-nright,-1
		! construct the pppVmat
		do j=1,smadim,1
			call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Rrealdim,4*Rrealdim, &
				LRcoeffin(:,j),LRcoeffincol(:,j),LRcoeffinrow(:,j), &
				moperamatbig1(:,i*3),mbigcolindex1(:,i*3),mbigrowindex1(:,i*3), &
				pppVmat,pppVmatcol,pppVmatrow,pppnelement,info)
			call checkinfo(info)

			do l=1,nleft+1,1
				! buffmat is to save the intermediate matrix
				call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Lrealdim,4*Rrealdim, &
					moperamatbig1(:,l*3),mbigcolindex1(:,l*3),mbigrowindex1(:,l*3), &
					pppVmat,pppVmatcol,pppVmatrow, &
					buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
				call checkinfo(info)
				! add LRcoeffout and buffmat
				call SpMatAdd(4*Rrealdim,4*Lrealdim,LRcoeffout(:,j),LRcoeffoutcol(:,j),LRcoeffoutrow(:,j),&
				'N',pppV(i,l),4*Rrealdim,4*Lrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
			end do
		end do
	end do

	deallocate(pppVmat,pppVmatrow,pppVmatcol)
	allocate(hopmat(hopnelement)) ! store the hopping matrix
	allocate(hopmatcol(hopnelement)) 
	allocate(hopmatrow(4*Lrealdim+1)) 

	!---------------------------------------------------------
	! the +1 -1 phase added to l' of hopmat
	allocate(phase(4*Lrealdim,2),stat=error)
	if(error/=0) stop

	do j=1,4*Lrealdim,1
		phase(j,1)=(-1.0D0)**(mod(quantabigL(j,1),2))
		phase(j,2)=(-1.0D0)*phase(j,1)
	end do

	!  hopping term
	do i=norbs,norbs-nright,-1
		! check if need hopping matrix
		ifhop=.false.
		do j=1,nleft+1,1
			if(bondlink(i,j)==1) then
				ifhop=.true.
				exit
			end if
		end do

		if(ifhop==.true.) then
			! construct hopmat
			do j=1,smadim,1
			do k=1,4,1
				if(k<=2) then
				! firstly needed to transfer from a(+) to a, so 'T' and 'T' is 'N', (ni-1)^+=(ni-1)
				! k=1 a up,k=2 a down,k=3 n,k=4 a+ up,k=5 a+ down;
					call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Rrealdim,4*Rrealdim, &
							LRcoeffin(:,j),LRcoeffincol(:,j),LRcoeffinrow(:,j), &
							moperamatbig1(:,(i-1)*3+k),mbigcolindex1(:,i*3-3+k),mbigrowindex1(:,i*3-3+k), &
							hopmat,hopmatcol,hopmatrow,hopnelement,info)
					call checkinfo(info)
					do l1=1,4*Lrealdim,1
					do l2=hopmatrow(l1),hopmatrow(l1+1)-1,1
						hopmat(l2)=hopmat(l2)*phase(l1,1)
					end do
					end do
				else
					call  SpMMtoSp('N','T',4*Lrealdim,4*Rrealdim,4*Rrealdim,4*Rrealdim,4*Lrealdim,&
							LRcoeffin(:,j),LRcoeffincol(:,j),LRcoeffinrow(:,j), &
							moperamatbig1(:,i*3-5+k),mbigcolindex1(:,i*3-5+k),mbigrowindex1(:,i*3-5+k), &
							hopmat,hopmatcol,hopmatrow,hopnelement)
					do l1=1,4*Lrealdim,1
					do l2=hopmatrow(l1),hopmatrow(l1+1)-1,1
						!transfer from al*ar^(+) to ar^(+)*al
						hopmat(l2)=hopmat(l2)*phase(l1,2)
					end do
					end do
				end if
			
				do l=1,nleft+1,1
				if(bondlink(i,l)==1) then
					!k<=2 al^+*ar,k>2 al*ar^(+) 
					if(k<=2) then
						call mkl_dcsrmultcsr('N',0,8,4*Lrealdim,4*Lrealdim,4*Rrealdim, &
							moperamatbig1(:,l*3-3+k),mbigcolindex1(:,l*3-3+k),mbigrowindex1(:,l*3-3+k), &
							hopmat,hopmatcol,hopmatrow, &
							buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
						call checkinfo(info)
					else
						call mkl_dcsrmultcsr('T',0,8,4*Lrealdim,4*Lrealdim,4*Rrealdim, &
							moperamatbig1(:,l*3-5+k),mbigcolindex1(:,l*3-5+k),mbigrowindex1(:,l*3-5+k), &
							hopmat,hopmatcol,hopmatrow, &
							buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
						call checkinfo(info)
					end if
					call SpMatAdd(4*Rrealdim,4*Lrealdim,LRcoeffout(:,j),LRcoeffoutcol(:,j),LRcoeffoutrow(:,j),&
					'N',t(i,l),4*Rrealdim,4*Lrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
				end if
				end do
			end do
			end do
		end if
	end do
	deallocate(phase)

	! transfer LRcoeffout to coeffnosymm
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

	newcoeff=0.0D0
	if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
		do i=1,smadim,1
			call symmetrizestate(ngoodstates,coeffnosymm(ngoodstates*(i-1)+1:i*ngoodstates),&
				newcoeff(bigdim*(i-1)+1:i*bigdim),'s')
		end do
	else
		call copy(coeffnosymm,newcoeff)
	end if

	deallocate(LRcoeffin,LRcoeffinrow,LRcoeffincol)
	deallocate(LRcoeffout,LRcoeffoutrow,LRcoeffoutcol)
	deallocate(buffmat,buffmatcol,buffmatrow)
	deallocate(hopmat,hopmatrow,hopmatcol)
	deallocate(coeffnosymm)
	
return
end subroutine masterop

!================================================================
!================================================================

subroutine MasterGather
	implicit none
	integer :: i

	call MasterGatherAllocateArray
	
	do i=1,nleft+1,1
		call ptopGather(i,'L')
	end do
	do i=norbs,norbs-nright,-1
		call ptopGather(i,'R')
	end do
	deallocate(packbuf)
return

end subroutine MasterGather

!================================================================
!================================================================

subroutine ptopGather(index1,subspace)
	implicit none
	character(len=1) :: subspace
	integer :: index1
	
	integer :: k
	integer :: dim1,position1,nonzero,location,mlocation
	integer :: status(MPI_STATUS_SIZE),ierr

	if(myid==orbid1(index1,1)) then
		if(subspace=="L") then
			dim1=Lrealdim
		else if(subspace=="R") then
			dim1=Rrealdim
		end if
		position1=0
		do k=1,3,1
			location=3*orbid1(index1,2)-3+k
			call MPI_PACK(bigrowindex1(1,location),(4*dim1+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			nonzero=bigrowindex1(4*dim1+1,location)-1
			call MPI_PACK(operamatbig1(1,location),nonzero,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(bigcolindex1(1,location),nonzero,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		end do
		call MPI_SEND(packbuf,position1,MPI_PACKED,0,index1,MPI_COMM_WORLD,ierr)
	else if(myid==0) then
		call MPI_RECV(packbuf,packsize,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)  ! get anyone is ok
		mlocation=status(MPI_TAG)
		if(mlocation<=nleft+1) then
			dim1=Lrealdim
		else
			dim1=Rrealdim
		end if
		position1=0
		do k=1,3,1
			call MPI_UNPACK(packbuf,packsize,position1,mbigrowindex1(1,3*mlocation-3+k),(4*dim1+1),MPI_integer4,MPI_COMM_WORLD,ierr)
			nonzero=mbigrowindex1(4*dim1+1,3*mlocation-3+k)-1
			call MPI_UNPACK(packbuf,packsize,position1,moperamatbig1(1,3*mlocation-3+k),nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
			call MPI_UNPACK(packbuf,packsize,position1,mbigcolindex1(1,3*mlocation-3+k),nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
		end do
	end if
return
end subroutine ptopGather

!================================================================
!================================================================

subroutine MasterGatherAllocateArray

	implicit none

	packsize=(bigdim1*12+4*(4*subM+1))*3
	allocate(packbuf(packsize))

return
end subroutine MasterGatherAllocateArray

!================================================================
!================================================================

end module masterdiag

