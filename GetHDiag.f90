module GetHdiag_mod
contains

Subroutine GetHDiag(HDIAGnosymm,num,&
cap_big,cap_bigcol,cap_bigrow,&
cap_Hbig,cap_Hbigcol,cap_Hbigrow,&
cap_quantabigL,cap_quantabigR,&
ifperturbation)
! This subroutine is to get the diagonal element of the Hamiltonian(no symmetry)
! from all the process which can be used in davidson diagonalization
! the hopping term did not contribute anyting to the diagnal term 
! because the number of electrons is not equal
! so the diagnol term only need the PPP term operator

! maxdim is the max{Lrealdim,Rrealdim}
	use mpi
	use variables
	use communicate
	use module_sparse

	implicit none
	
	integer,intent(in) :: num
	real(kind=r8),intent(out) :: HDIAGnosymm(num)
	real(kind=r8),intent(in) :: cap_big(:,:),cap_Hbig(:,:)
	integer(kind=i4),intent(in) :: &
		cap_bigcol(:,:),cap_bigrow(:,:),&
		cap_Hbigcol(:,:),cap_Hbigrow(:,:),&
		cap_quantabigL(:,:),cap_quantabigR(:,:)
	logical,intent(in) :: ifperturbation
	
	! local
	real(kind=8),allocatable :: buffermat(:),buffermat0(:,:),Hdiagdummy(:)
	integer :: operaindex
	integer :: status(MPI_STATUS_SIZE),ierr
	integer :: i,error,j,k,m
	integer :: iLrealdim,iRrealdim
	integer :: maxdim


	call master_print_message("enter in GetHDiag subroutine")
	
	if(ifperturbation==.true.) then
		iLrealdim=Lrealdimp
		iRrealdim=Rrealdimp
	else
		iLrealdim=Lrealdim
		iRrealdim=Rrealdim
	end if
	maxdim=max(iLrealdim,iRrealdim)
	
	if(myid/=0) then
		allocate(buffermat(4*maxdim),stat=error)
		if(error/=0) stop
		buffermat=0.0D0
	else
		allocate(buffermat0(4*maxdim,norbs),stat=error)
		if(error/=0) stop
		buffermat0=0.0D0
		allocate(Hdiagdummy(16*iLrealdim*iRrealdim),stat=error)
		if(error/=0) stop
	end if

	! L space
	do i=1,nleft+1,1
		if(myid==orbid1(i,1)) then
			operaindex=orbid1(i,2)
			! copy the diagonal element to the buffermat
			do j=1,4*iLrealdim,1
				do k=cap_bigrow(j,operaindex*3),cap_bigrow(j+1,operaindex*3)-1,1
					if(cap_bigcol(k,operaindex*3)==j) then
						buffermat(j)=cap_big(k,operaindex*3)
						exit
					end if
				end do
			end do
			! send the diagonal element to the 0 process
			call MPI_SEND(buffermat,4*maxdim,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(1,i),4*maxdim,mpi_real8,orbid1(i,1),i,MPI_COMM_WORLD,status,ierr)
		end if
	end do

	! R space
	do j=norbs,norbs-nright,-1
		if(myid==orbid1(j,1)) then
			operaindex=orbid1(j,2)
			! copy the diagonal element to the buffermat
			do i=1,4*iRrealdim,1
				do k=cap_bigrow(i,operaindex*3),cap_bigrow(i+1,operaindex*3)-1,1
					if(cap_bigcol(k,operaindex*3)==i) then
						buffermat(i)=cap_big(k,operaindex*3)
						exit
					end if
				end do
			end do
			! send the diagonal element to the 0 process
			call MPI_SEND(buffermat,4*maxdim,mpi_real8,0,j,MPI_COMM_WORLD,ierr)
		else if(myid==0) then
			call MPI_RECV(buffermat0(1,j),4*maxdim,mpi_real8,orbid1(j,1),j,MPI_COMM_WORLD,status,ierr)
		end if
	end do

	if(myid==0) then
		Hdiagdummy=0.0D0
		! HL contribution
		do i=1,4*iRrealdim,1
			do j=1,4*iLrealdim,1
				do k=cap_Hbigrow(j,1),cap_Hbigrow(j+1,1)-1,1
					if(cap_Hbigcol(k,1)==j) then
						Hdiagdummy((i-1)*4*iLrealdim+j)=cap_Hbig(k,1)
						exit
					end if
				end do
			end do
		end do
		! HR contribution
		do i=1,4*iRrealdim,1
			do k=cap_Hbigrow(i,2),cap_Hbigrow(i+1,2)-1,1
				if(cap_Hbigcol(k,2)==i) then
					Hdiagdummy((i-1)*4*iLrealdim+1:i*4*iLrealdim)=cap_Hbig(k,2)+Hdiagdummy((i-1)*4*iLrealdim+1:i*4*iLrealdim)
					exit
				end if
			end do
		end do
		! transfer integral contribution the contribute is zero
		! PPP term contribution
		if(logic_PPP==1) then
			do i=1,nleft+1,1
			do j=norbs,norbs-nright,-1
				do k=1,4*iRrealdim,1
					Hdiagdummy((k-1)*4*iLrealdim+1:k*4*iLrealdim)=buffermat0(1:4*iLrealdim,i)*buffermat0(k,j)*pppV(i,j)&
											+Hdiagdummy((k-1)*4*iLrealdim+1:k*4*iLrealdim)
				end do
			end do
			end do
		end if
		! copy Hdiagdummy to Hdiag
		! and ignore those diag element corresponding states without good quantum number
		m=1
		do i=1,4*iRrealdim,1
		do j=1,4*iLrealdim,1
			if((cap_quantabigL(j,1)+cap_quantabigR(i,1)==nelecs) .and. &
				cap_quantabigL(j,2)+cap_quantabigR(i,2)==totalSz) then
				HDIAGnosymm(m)=Hdiagdummy((i-1)*4*iLrealdim+j)
				m=m+1
			end if
		end do
		end do
		
		m=m-1
		if(m/=num) then
			call master_print_message(m,"HDIAGnosymm number wrong! failed!,m=")
		end if
	end if

	if(myid==0) then
		deallocate(buffermat0)
		deallocate(Hdiagdummy)
	else
		deallocate(buffermat)
	end if

return

end subroutine GetHDiag

end module GetHdiag_mod
