Module Renormalization_mod
! after diaganolizaiton we need to renormalizaiton the many body states
! in fact only renormalization all the operator matrix

	use variables
	use communicate
	use kinds_mod
	use module_sparse

	implicit none
	save
	private

	public :: Renormalization

	real(kind=r8),allocatable ::  &
	leftu(:,:) , &
	rightv(:,:)  , &
	singularvalue(:)

contains

!===================================================
!===================================================
Subroutine Renormalization(direction)
! direction=l means l block is the system
! direction=i means is the infinit MPS
! direction=r means r block is the system
	use mpi
	implicit none
	
	character(len=1) :: direction
	integer :: error,ierr

	call master_print_message("enter Renormalization subroutine")

	if(4*Lrealdim>subM .or. 4*Rrealdim>subM) then
		if(myid==0) then
			! allocate work array
			allocate(singularvalue(subM),stat=error)
			if(error/=0) stop
			allocate(leftu(4*Lrealdim,subM),stat=error)
			if(error/=0) stop
			allocate(rightv(subM,4*Rrealdim),stat=error)
			if(error/=0) stop
			
			! get rotate matrix leftu/rightv and singularvalue
			if(nstate==1 .or. exscheme==1) then
				! gs or state-average exscheme
				call splitsvdL(singularvalue,leftu,1,nstate,nleft+1)
				call splitsvdR(singularvalue,rightv,1,nstate,norbs-nright)
			else if(exscheme==2) then
				! my new exScheme
				call ExScheme2
			end if
			
			! store the wavefunction
			call StoreWaveFunction
		end if

		! bcast the updated quantasmaL/R because only the 0 process know it
		call MPI_BCAST(quantasmaL(1,1),subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(quantasmaR(1,1),subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
		
		! rotate the bigmat to smamat
		! L subspace
		if(direction=='i' .or. direction=='l') then
			call RotateBasis('L')
		end if
		! R subspace
		if(direction=='i' .or. direction=='r') then
			call RotateBasis('R')
		end if
	else
		! direct copy bigmat to smamat
		call DirectCopy('L')
		call DirectCopy('R')
	end if

	! destroy work array
	if(4*Lrealdim>subM .or. 4*Rrealdim>subM) then
		if(myid==0) then
			deallocate(singularvalue)
			deallocate(leftu)
			deallocate(rightv)
		end if
	end if
return

end subroutine Renormalization

!===================================================
!===================================================

subroutine ExScheme2
! my ExScheme2 method to target excited states

	implicit none

	real(kind=r8),allocatable :: &
	singualarvalue2(:) , &
	leftu2(:,:)        , &
	rightv2(:,:)        
	integer(kind=i4),allocatable :: &
	quantasmaL2(:,:) , &
	quantasmaR2(:,:) , &
	symmlinksma2(:,:,:)
	integer :: error
	integer :: i

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
	do i=1,nstate,1
		call splitsvdL(singularvalue((i-1)*subM+1:i*subM),leftu2(:,(i-1)*subM+1:i*subM),i,i,nleft+1)
		call splitsvdR(singularvalue((i-1)*subM+1:i*subM),rightv2((i-1)*subM+1:i*subM,:),i,i,norbs-nright)
		quantasmaL2((i-1)*subM+1:i*subM,:)=quantasmaL
		quantasmaR2((i-1)*subM+1:i*subM,:)=quantasmaR
		if(logic_spinreversal/=0) then
			symmlinksma2((i-1)*subM+1:i*subM,:,:)=symmlinksma
		end if
	end do
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

return
end Subroutine ExScheme2

!===================================================
!===================================================

subroutine StoreWaveFunction
! store the wavefunction in matrix product form
! and singularvalue

	implicit none

	integer :: reclength
	logical :: alive
	integer :: i
	
	! 4 byte as 1 direct file block
	reclength=2*subM*subM
	
	! wavefunction.tmp
	inquire(file="wavefunction.tmp",exist=alive)
	if(alive) then
		open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
	else
		open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
	end if

	! singularvalue.tmp
	open(unit=106,file="singularvalue.tmp",status="replace")

	! divid the leftu and rightv to small M*M matrix
	! leftu
	do i=1,4,1
		write(105,rec=4*nleft+i) leftu((i-1)*Lrealdim+1:i*Lrealdim,1:subM)
	end do
	! rightv
	do i=1,4,1
		write(105,rec=4*(norbs-nright-1)+i) rightv(1:subM,i:4*Rrealdim:4)
	end do

	! write singularvalue though only used in finit MPS
	! the singularvalue here is the exactly singlarvalue^2
	write(106,*) singularvalue
	close(105)
	close(106)
return

end subroutine StoreWaveFunction

!===================================================
!===================================================

subroutine RotateBasis(domain)
! Rotate Basis ; In fact transfer the operator matrix
! to new basis
	
	use exit_mod
	use mathlib
	use mpi
	implicit none
	include "mkl_spblas.fi"

	character(len=1) :: domain
	! local
	real(kind=r8),allocatable ::     &
			rotatematdens(:,:) , &     ! store leftu and rightv
			rotatemat(:)               ! store the CSR sparse format
	integer :: orbstart,orbend,Hindex,dim1
	integer :: arraylength,nrows,operaindex
	integer(kind=i4),allocatable :: rotcolindex(:),rotrowindex(:)
	integer :: job(8)
	integer :: error,ierr,info
	integer :: i,j

	character(len=1),allocatable :: packbuf(:)
	integer :: packsize
	integer :: position1
	
	
	if(domain=='L') then
		orbstart=1
		orbend=nleft+1
		Hindex=1
		dim1=Lrealdim
	else if(domain=='R') then
		orbstart=norbs-nright
		orbend=norbs
		Hindex=2
		dim1=Rrealdim
	else
		call exit_DMRG(sigAbort,"RotateBasis domain/=L/R")
	end if

	! define the U and V nonzero element numbers
	arraylength=NINT(DBLE(4*subM*subM)/UVmatratio,i4)

	! leftu rightv store in CSR format
	allocate(rotatemat(arraylength),stat=error)
	if(error/=0) stop
	allocate(rotcolindex(arraylength),stat=error)
	if(error/=0) stop
	allocate(rotrowindex(4*dim1+1),stat=error) 
	if(error/=0) stop

	packsize=arraylength*12+4*(4*dim1+1)+1000
	allocate(packbuf(packsize),stat=error) 
	if(error/=0) stop

	if(myid==0) then
		! store the dense U/V
		allocate(rotatematdens(4*dim1,subM),stat=error)
		if(error/=0) stop
		if(domain=='L') then
			rotatematdens=leftu
		else if(domain=='R') then
			rotatematdens=transpose(rightv)
		end if
		job(1)=0
		job(2)=1
		job(3)=1
		job(4)=2
		job(5)=arraylength
		job(6)=1
		nrows=4*dim1
		call mkl_ddnscsr(job,nrows,subM,rotatematdens,nrows,rotatemat,rotcolindex,rotrowindex,info)
		if(info==0) then
		!	call master_print_message(arraylength,"U/V maxnelement=")
		!	call master_print_message(rotrowindex(nrows+1)-1,"U/V nonzero=")
		else
			call master_print_message(info,"info/=0 dens to sparse in rotatebasis")
			stop
		end if
		deallocate(rotatematdens)

		position1=0
		call MPI_PACK(rotrowindex,4*dim1+1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(rotatemat,rotrowindex(4*dim1+1)-1,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(rotcolindex,rotrowindex(4*dim1+1)-1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
	end if
	
	! broadcast the rotate matrix
	call MPI_BCAST(position1,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(packbuf,position1,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
	
	if(myid/=0) then
		position1=0
		call MPI_UNPACK(packbuf,packsize,position1,rotrowindex,4*dim1+1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,rotatemat,rotrowindex(4*dim1+1)-1,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,rotcolindex,rotrowindex(4*dim1+1)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
	end if
	
	! rotate the matrix
	do i=orbstart,orbend,1
		if(myid==orbid1(i,1)) then
			do j=1,3,1
				operaindex=orbid1(i,2)*3-3+j
				call SpMatRotateBasis(subM,4*dim1,rotatemat,rotcolindex,rotrowindex, &
							4*dim1,4*dim1,operamatbig1(:,operaindex),bigcolindex1(:,operaindex),bigrowindex1(:,operaindex), &
							subM,operamatsma1(:,operaindex),smacolindex1(:,operaindex),smarowindex1(:,operaindex),smadim1)
			end do
		end if
	end do
	
	! rotate HL and HR
	if(myid==0) then
		call SpMatRotateBasis(subM,4*dim1,rotatemat,rotcolindex,rotrowindex, &
					4*dim1,4*dim1,Hbig(:,Hindex),Hbigcolindex(:,Hindex),Hbigrowindex(:,Hindex), &
					subM,Hsma(:,Hindex),Hsmacolindex(:,Hindex),Hsmarowindex(:,Hindex),Hsmadim)
	end if

	deallocate(rotatemat)
	deallocate(rotcolindex)
	deallocate(rotrowindex)
	deallocate(packbuf)
return

end subroutine RotateBasis

!===================================================
!===================================================

subroutine DirectCopy(domain)
! direct copy the operamatbig to operamatsma
! only in the infinit MPS process
! operamatbig -> operamatsma 
! Hbig -> Hsma
! symmlinkbig -> symmlinksma
	
	use exit_mod
	implicit none

	character(len=1) domain
	
	! local
	integer :: i,j
	integer :: operaindex
	integer :: orbstart,orbend,Hindex,dim1
	integer :: nelement

	if(domain=='L') then
		orbstart=1
		orbend=nleft+1
		Hindex=1
		dim1=Lrealdim
	else if(domain=='R') then
		orbstart=norbs-nright
		orbend=norbs
		Hindex=2
		dim1=Rrealdim
	else
		call exit_DMRG(sigAbort,"DirectCopy domain/=L/R")
	end if

	do i=orbstart,orbend,1
		if(myid==orbid1(i,1)) then
			do j=1,3,1
				operaindex=orbid1(i,2)*3-3+j
				smarowindex1(:,operaindex)=0
				smarowindex1(1:4*dim1+1,operaindex)=bigrowindex1(1:4*dim1+1,operaindex)
				nelement=smarowindex1(4*dim1+1,operaindex)-1
				smacolindex1(1:nelement,operaindex)=bigcolindex1(1:nelement,operaindex)
				operamatsma1(1:nelement,operaindex)=operamatbig1(1:nelement,operaindex)
			end do
		end if
	end do

	if(myid==0) then
		Hsmarowindex(:,Hindex)=0
		Hsmarowindex(1:4*dim1+1,Hindex)=Hbigrowindex(1:4*dim1+1,Hindex)
		nelement=Hsmarowindex(4*dim1+1,Hindex)-1
		Hsmacolindex(1:nelement,Hindex)=Hbigcolindex(1:nelement,Hindex)
		Hsma(1:nelement,Hindex)=Hbig(1:nelement,Hindex)
		if(logic_spinreversal/=0) then
			symmlinksma(1:4*dim1,1,Hindex)=symmlinkbig(1:4*dim1,1,Hindex)
		end if
	end if

	if(domain=='L') then
		quantasmaL(1:4*Lrealdim,:)=quantabigL(1:4*Lrealdim,:)
	else
		quantasmaR(1:4*Rrealdim,:)=quantabigR(1:4*Rrealdim,:)
	end if

return

end subroutine DirectCopy

!===================================================
!===================================================
end Module Renormalization_mod
