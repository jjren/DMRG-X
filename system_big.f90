Subroutine System_Big(domain)

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
!	include "mkl.fi"

	character(len=1) :: domain   !L/R
	
	! local
	integer :: orbstart,orbend,orbadd,Hindex,dim1
	! orbstart is from 1 or norbs-nright+1
	! orbend is nleft or norbs
	! orbadd is nleft+1 or norbs-nright
	! Hindex : HL/1 HR/2
	integer :: operaindex
	real(kind=r8) :: II(4),IM(subM)
	integer(kind=i4) :: IIrowindex(5),IIcolindex(4),IMrowindex(subM+1),IMcolindex(subM)

	integer(kind=i4),allocatable :: phase(:)

	real(kind=r8),allocatable :: Hbuf(:),Hbuf2(:),obuf(:,:)  ! Hbuffer and operator buffer
	integer(kind=i4),allocatable :: Hbufrowindex(:),Hbufcolindex(:),&
						Hbufrowindex2(:),Hbufcolindex2(:),&
						obufrowindex(:,:),obufcolindex(:,:)
	integer :: info
	integer :: i,j,k
	integer :: error

	integer :: status(MPI_STATUS_SIZE),recvtag,recvrequest ! MPI flag
	character(len=1),allocatable :: packbuf1(:),packbuf2(:)
	integer(kind=i4) :: packsize1,packsize2    
	integer(kind=i4) :: position1
	integer :: ierr

	call master_print_message("enter in subroutine system_big")

	packsize1=smadim1*72+12*subM+1000   ! store hopping matrix
	allocate(packbuf1(packsize1),stat=error)
	if(error/=0) stop
	packsize2=smadim1*24+4*subM+1000    ! store ppp matrix
	allocate(packbuf2(packsize2),stat=error)
	if(error/=0) stop
	call master_print_message(packsize1,"packsize1=")
	call master_print_message(packsize2,"packsize2=")


	if(domain=='L') then
		orbstart=1
		orbend=nleft
		orbadd=nleft+1
		dim1=Lrealdim
		Hindex=1
	else if(domain=='R') then
		orbstart=norbs-nright+1
		orbend=norbs
		orbadd=norbs-nright
		dim1=Rrealdim
		Hindex=2
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
	end if

!==================================================================
! set phase value
	allocate(phase(4*subM),stat=error)
	if(error/=0) stop
	if(domain=='R' .and. logic_C2==0) then
		phase(1:4*dim1:4)=1
		phase(2:4*dim1:4)=-1
		phase(3:4*dim1:4)=-1
		phase(4:4*dim1:4)=1
	else
		do j=1,4*dim1,1
			if(mod(j,dim1)==0) then
				k=dim1
			else
				k=mod(j,dim1)
			end if
			if(domain=='L') then
				phase(j)=(-1)**(mod(quantasmaL(k,1),2))
			else
				phase(j)=(-1)**(mod(quantasmaR(k,1),2))
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
	do i=1,subM,1
		IMrowindex(i)=i
		IMcolindex(i)=i
	end do
	IMrowindex(subM+1)=subM+1

!=========================================================================

! construct the L/R(without sigmaL/R) subspace operator matrix in 4M basis

	do i=orbstart,orbend,1
	if(myid==orbid1(i,1)) then
		! send the PPP operator small matrix
		if(bondlink(i,orbadd)==1) then
			position1=0
			call MPI_PACK(smarowindex1(1,(orbid1(i,2)*3-2)),(subM+1)*3,MPI_integer4,packbuf1,packsize1,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(smacolindex1(1,(orbid1(i,2)*3-2)),3*smadim1,MPI_integer4,packbuf1,packsize1,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(operamatsma1(1,(orbid1(i,2)*3-2)),3*smadim1,mpi_real8,packbuf1,packsize1,position1,MPI_COMM_WORLD,ierr)
			call MPI_SEND(packbuf1,packsize1,MPI_PACKED,0,i,MPI_COMM_WORLD,ierr)
		else
			position1=0
			call MPI_PACK(smarowindex1(1,(orbid1(i,2)*3)),(subM+1),MPI_integer4,packbuf2,packsize2,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(smacolindex1(1,(orbid1(i,2)*3)),smadim1,MPI_integer4,packbuf2,packsize2,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(operamatsma1(1,(orbid1(i,2)*3)),smadim1,mpi_real8,packbuf2,packsize2,position1,MPI_COMM_WORLD,ierr)
			call MPI_SEND(packbuf2,packsize2,MPI_PACKED,0,i,MPI_COMM_WORLD,ierr)
		end if
		write(*,*) "position1=",position1
		
		operamatbig1(:,orbid1(i,2)*3-2:orbid1(i,2)*3)=0.0D0
		bigrowindex1(:,orbid1(i,2)*3-2:orbid1(i,2)*3)=0
		bigcolindex1(:,orbid1(i,2)*3-2:orbid1(i,2)*3)=0

		if(domain=='R' .and. logic_C2==0 ) then
			do j=1,3,1
				operaindex=(orbid1(i,2)-1)*3+j
				call SparseDirectProduct(4,4,II,IIcolindex,IIrowindex,&
								dim1,dim1,operamatsma1(:,operaindex),smacolindex1(:,operaindex),smarowindex1(:,operaindex),&
								operamatbig1(:,operaindex),bigcolindex1(:,operaindex),bigrowindex1(:,operaindex),bigdim1)
				if(j<=2) then
					do k=1,bigrowindex1(4*dim1+1,operaindex)-1,1
						operamatbig1(k,operaindex)=operamatbig1(k,operaindex)*DBLE(phase(bigcolindex1(k,operaindex)))
					end do
				end if
			end do
		else
			do j=1,3,1
				operaindex=(orbid1(i,2)-1)*3+j
				call SparseDirectProduct(dim1,dim1,operamatsma1(:,operaindex),smacolindex1(:,operaindex),smarowindex1(:,operaindex),&
								4,4,II,IIcolindex,IIrowindex,&
								operamatbig1(:,operaindex),bigcolindex1(:,operaindex),bigrowindex1(:,operaindex),bigdim1)
			end do
		end if
	end if
	end do

!==============================================================================================================

! construct the sigmaL/R subspace operator matrix in 4M basis
	
	if(myid==orbid1(orbadd,1)) then
		operamatbig1(:,orbid1(orbadd,2)*3-2:orbid1(orbadd,2)*3)=0.0D0
		bigrowindex1(:,orbid1(orbadd,2)*3-2:orbid1(orbadd,2)*3)=0
		bigcolindex1(:,orbid1(orbadd,2)*3-2:orbid1(orbadd,2)*3)=0
		
		if(domain=='R' .and. logic_C2==0) then
			do j=1,3,1
				operaindex=(orbid1(orbadd,2)-1)*3+j
				call SparseDirectProduct(4,4,onesitemat(:,j),osmcolindex(:,j),osmrowindex(:,j),&
								dim1,dim1,IM,IMcolindex,IMrowindex,&
								operamatbig1(:,operaindex),bigcolindex1(:,operaindex),bigrowindex1(:,operaindex),bigdim1)
			end do
		else 
			do j=1,3,1
				operaindex=(orbid1(orbadd,2)-1)*3+j
				call SparseDirectProduct(dim1,dim1,IM,IMcolindex,IMrowindex,&
								4,4,onesitemat(:,j),osmcolindex(:,j),osmrowindex(:,j),&
								operamatbig1(:,operaindex),bigcolindex1(:,operaindex),bigrowindex1(:,operaindex),bigdim1)
				if(j<=2) then
					do k=1,bigrowindex1(4*dim1+1,operaindex)-1,1
						operamatbig1(k,operaindex)=operamatbig1(k,operaindex)*DBLE(phase(bigcolindex1(k,operaindex)))
					end do
				end if
			end do
		end if
	end if

! cosntruct the L+sigmaL/R+sigmaR Hamiltonian operator in 4M basis
	if(myid==0) then

		call MPI_IRECV(packbuf1,packsize1,mpi_packed,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,recvrequest,ierr)
		
		allocate(obuf(smadim1,3),stat=error)
		if(error/=0) stop
		allocate(obufrowindex(subM+1,3),stat=error)
		if(error/=0) stop
		allocate(obufcolindex(smadim1,3),stat=error)
		if(error/=0) stop
		
		allocate(Hbuf(Hbigdim),stat=error)
		if(error/=0) stop
		allocate(Hbufrowindex(4*subM+1),stat=error)
		if(error/=0) stop
		allocate(Hbufcolindex(Hbigdim),stat=error)
		if(error/=0) stop

		allocate(Hbuf2(Hbigdim),stat=error)
		if(error/=0) stop
		allocate(Hbufrowindex2(4*subM+1),stat=error)
		if(error/=0) stop
		allocate(Hbufcolindex2(Hbigdim),stat=error)
		if(error/=0) stop


!===========================================================

!     L/R Hamiltonian contribute
		Hbuf=0.0D0
		Hbufrowindex=0
		Hbufcolindex=0

		if(domain=='R' .and. logic_C2==0) then
			call SparseDirectProduct(4,4,II,IIcolindex,IIrowindex,&
							dim1,dim1,Hsma(:,Hindex),Hsmacolindex(:,Hindex),Hsmarowindex(:,Hindex),&
							Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
		else
			call SparseDirectProduct(dim1,dim1,Hsma(:,Hindex),Hsmacolindex(:,Hindex),Hsmarowindex(:,Hindex),&
							4,4,II,IIcolindex,IIrowindex,&
							Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
		end if

		Hbig(:,Hindex)=0.0D0
		Hbigcolindex(:,Hindex)=0
		Hbigrowindex(:,Hindex)=0
		call copy(Hbuf,Hbig(:,Hindex))
		call copy(Hbufrowindex,Hbigrowindex(:,Hindex))
		call copy(Hbufcolindex,Hbigcolindex(:,Hindex))

!===========================================================

!     sigmaL Hamiltonian contribute. site energy+HubbardU
		Hbuf=0.0D0
		Hbufrowindex=0
		Hbufcolindex=0

		if(domain=='R' .and. logic_C2==0) then
			call SparseDirectProduct(4,4,onesitemat(:,6),osmcolindex(:,6),osmrowindex(:,6),&
							dim1,dim1,IM,IMcolindex,IMrowindex,&
							Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
		else
			call SparseDirectProduct(dim1,dim1,IM,IMcolindex,IMrowindex,&
							4,4,onesitemat(:,6),osmcolindex(:,6),osmrowindex(:,6),&
							Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
		end if
		Hbuf2=0.0D0
		Hbufrowindex2=0
		Hbufcolindex2=0
		call mkl_dscradd('N',0,0,4*subM,4*subM,Hbig(:,Hindex),Hbigcolindex(:,Hindex),Hbigrowindex(:,Hindex),1.0D0,Hbuf,Hbufcolindex,Hbufrowindex,&
					Hbuf2,Hbuf2colindex,Hbuf2rowindex,Hbigdim,info)
		if(info/=0) then
			call master_print_message("info/=0 failed!")
			stop
		end if
		Hbig(:,Hindex)=0.0D0
		Hbigcolindex(:,Hindex)=0
		Hbigrowindex(:,Hindex)=0
		call copy(Hbuf2,Hbig(:,Hindex))
		call copy(Hbufrowindex2,Hbigrowindex(:,Hindex))
		call copy(Hbufcolindex2,Hbigcolindex(:,Hindex))

!===========================================================

!     L sigmaL/R sigmaR interaction operator contribute

		do i=orbstart,orbend,1
		! nonblock
			call MPI_WAIT(recvrequest,status,ierr)
			recvtag=status(MPI_TAG)

			if(bondlink(recvtag,orbadd)==1) then
				position1=0
				call MPI_UNPACK(packbuf1,packsize1,position1,obufrowindex(1,1),(subM+1)*3,MPI_integer4,MPI_COMM_WORLD,ierr)
				call MPI_UNPACK(packbuf1,packsize1,position1,obufcolindex(1,1),3*smadim1,MPI_integer4,MPI_COMM_WORLD,ierr)
				call MPI_UNPACK(packbuf1,packsize1,position1,obuf(1,1),3*smadim1,mpi_real8,MPI_COMM_WORLD,ierr)
			else
				position1=0
				call MPI_UNPACK(packbuf1,packsize1,position1,obufrowindex(1,3),(subM+1),MPI_integer4,MPI_COMM_WORLD,ierr)
				call MPI_UNPACK(packbuf1,packsize1,position1,obufcolindex(1,3),smadim1,MPI_integer4,MPI_COMM_WORLD,ierr)
				call MPI_UNPACK(packbuf1,packsize1,position1,obuf(1,3),smadim1,mpi_real8,MPI_COMM_WORLD,ierr)
			end if

			if(i<orbend) then
				call MPI_IRECV(packbuf1,packsize1,mpi_packed,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,recvrequest,ierr)
			end if
			
		!transfer integral term
			if(bondlink(recvtag,orbadd)==1) then
			do j=1,2,1

				Hbuf=0.0D0
				Hbufrowindex=0
				Hbufcolindex=0

				if(domain=='R' .and. logic_C2==0) then
					call SparseDirectProduct(4,4,onesitemat(:,j+3),osmcolindex(:,j+3),osmrowindex(:,j+3),&
									dim1,dim1,obuf(:,j),obufcolindex(:,j),obufrowindex(:,j),&
									Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
					do k=1,Hbufrowindex(4*dim1+1)-1,1
						Hbuf(k)=Hbuf(k)*DBLE(phase(Hbufcolindex(k)))*(-1.0D0)
					end do
				else
					call SparseDirectProduct(dim1,dim1,obuf(:,j),obufcolindex(:,j),obufrowindex(:,j),&
									4,4,onesitemat(:,j+3),osmcolindex(:,j+3),osmrowindex(:,j+3),&
									Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
					do k=1,Hbufrowindex(4*dim1+1)-1,1
						Hbuf(k)=Hbuf(k)*DBLE(phase(Hbufcolindex(k)))
					end do
				end if

				Hbuf2=0.0D0
				Hbufrowindex2=0
				Hbufcolindex2=0
				call mkl_dscradd('N',0,0,4*subM,4*subM,Hbig(:,Hindex),Hbigcolindex(:,Hindex),Hbigrowindex(:,Hindex),t(recvtag,orbadd),Hbuf,Hbufcolindex,Hbufrowindex,&
							Hbuf2,Hbuf2colindex,Hbuf2rowindex,Hbigdim,info)
				if(info/=0) then
					call master_print_message(info,"info/=0 failed!")
					stop
				end if

				Hbig(:,Hindex)=0.0D0
				Hbigcolindex(:,Hindex)=0
				Hbigrowindex(:,Hindex)=0
				call mkl_dscradd('T',0,0,4*subM,4*subM,Hbuf2,Hbuf2colindex,Hbuf2rowindex,t(recvtag,orbadd),Hbuf,Hbufcolindex,Hbufrowindex,&
							Hbig(:,Hindex),Hbigcolindex(:,Hindex),Hbigrowindex(:,Hindex),Hbigdim,info)
				if(info/=0) then
					call master_print_message(info,"info/=0 failed!")
					stop
				end if
			end do
			end if

		! ppp term
			Hbuf=0.0D0
			Hbufrowindex=0
			Hbufcolindex=0

			if(domain=='R' .and. logic_C2==0) then
				call SparseDirectProduct(4,4,onesitemat(:,3),osmcolindex(:,3),osmrowindex(:,3),&
								dim1,dim1,obuf(:,3),obufcolindex(:,3),obufrowindex(:,3),&
								Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
			else
				call SparseDirectProduct( dim1,dim1,obuf(:,3),obufcolindex(:,3),obufrowindex(:,3),&
								4,4,onesitemat(:,3),osmcolindex(:,3),osmrowindex(:,3),&
								Hbuf,Hbufcolindex,Hbufrowindex,Hbigdim)
			end if

			Hbuf2=0.0D0
			Hbufrowindex2=0
			Hbufcolindex2=0
			call mkl_dscradd('N',0,0,4*subM,4*subM,Hbig(:,Hindex),Hbigcolindex(:,Hindex),Hbigrowindex(:,Hindex),pppV(recvtag,orbadd),Hbuf,Hbufcolindex,Hbufrowindex,&
						Hbuf2,Hbuf2colindex,Hbuf2rowindex,Hbigdim,info)
			if(info/=0) then
				call master_print_message(info,"info/=0 failed!")
				stop
			end if

			Hbig(:,Hindex)=0.0D0
			Hbigcolindex(:,Hindex)=0
			Hbigrowindex(:,Hindex)=0

			call copy(Hbuf2,Hbig(:,Hindex))
			call copy(Hbufrowindex2,Hbigrowindex(:,Hindex))
			call copy(Hbufcolindex2,Hbigcolindex(:,Hindex))
		end do
!=========================================================================
	
	! construct the symmmlinkbig
		if(logic_spinreversal/=0) then
			call Creatsymmlinkbig(dim1,domain,Hindex)
		end if
		
		deallocate(Hbuf,Hbuf2)
		deallocate(obuf)
		deallocate(Hbufrowindex,Hbufcolindex)
		deallocate(Hbufrowindex2,Hbufcolindex2)
		deallocate(obufrowindex,obufcolindex)
	end if

	deallocate(packbuf1,packbuf2)
	deallocate(phase)
	return
end subroutine System_Big

