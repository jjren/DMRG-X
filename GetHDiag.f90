module GetHdiag_mod
contains

Subroutine GetHDiag(HDIAGnosymm,num,&
cap_big,cap_bigcol,cap_bigrow,&
cap_Hbig,cap_Hbigcol,cap_Hbigrow,&
cap_goodbasis,&
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
    USE BLAS95
    USE F95_PRECISION
    use mathlib
    implicit none
    
    integer,intent(in) :: num
    real(kind=r8),intent(out) :: HDIAGnosymm(num)
    real(kind=r8),intent(in) :: cap_big(:,:),cap_Hbig(:,:)
    integer(kind=i4),intent(in) :: &
        cap_bigcol(:,:),cap_bigrow(:,:),&
        cap_Hbigcol(:,:),cap_Hbigrow(:,:),&
        cap_goodbasis(:,:)
    logical,intent(in) :: ifperturbation
    
    ! local
    real(kind=r8),allocatable :: buffermat(:),buffermat0(:,:),Hdiagdummy(:)
    integer :: operaindex
    integer :: status(MPI_STATUS_SIZE),ierr
    integer :: i,error,j,k,m
    integer :: iLrealdim,iRrealdim
    integer :: maxdim
    real(kind=r8) :: alpha
    real(kind=r8),allocatable :: Houtput(:)
    
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
        allocate(Houtput(4*maxdim))
        do j=1,4*iLrealdim,1
            call SpMatIJ(4*iLrealdim,j,j,cap_Hbig(:,1),&
                cap_Hbigcol(:,1),cap_Hbigrow(:,1),Houtput(j))
        end do
        do i=1,4*iRrealdim,1
            call copy(Houtput(1:4*iLrealdim),Hdiagdummy((i-1)*4*iLrealdim+1:i*4*iLrealdim))
        end do

        ! HR contribution
        do i=1,4*iRrealdim,1
            call SpMatIJ(4*iRrealdim,i,i,cap_Hbig(:,2),&
                cap_Hbigcol(:,2),cap_Hbigrow(:,2),Houtput(i))
            Hdiagdummy((i-1)*4*iLrealdim+1:i*4*iLrealdim)=Houtput(i)+Hdiagdummy((i-1)*4*iLrealdim+1:i*4*iLrealdim)
        end do
        deallocate(Houtput)
        
        ! transfer integral contribution the contribute is zero
        ! PPP term contribution
        if(logic_PPP==1) then
            do i=1,nleft+1,1
            do j=norbs,norbs-nright,-1
                do k=1,4*iRrealdim,1
                   ! Hdiagdummy((k-1)*4*iLrealdim+1:k*4*iLrealdim)=buffermat0(1:4*iLrealdim,i)*buffermat0(k,j)*pppV(i,j)&
                   !                         +Hdiagdummy((k-1)*4*iLrealdim+1:k*4*iLrealdim)
                   alpha=buffermat0(k,j)*pppV(i,j)
                   call axpy(buffermat0(1:4*iLrealdim,i),Hdiagdummy((k-1)*4*iLrealdim+1:k*4*iLrealdim),alpha)
                end do
            end do
            end do
        end if
        ! copy Hdiagdummy to Hdiag
        ! and ignore those diag element corresponding states without good quantum number
        
        do i=1,num,1
            HDIAGnosymm(i)=Hdiagdummy((cap_goodbasis(i,2)-1)*4*iLrealdim+cap_goodbasis(i,1))
        end do
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
