module GetHdiag_mod
! This subroutine is to get the diagonal element of the Hamiltonian(no symmetry)
! from all the process which can be used in davidson diagonalization
! the hopping term did not contribute anyting to the diagnal term 
! because the number of electrons is not equal in the same L/R block
! so the diagnol term only need the PPP term operator
    USE MPI
    use kinds_mod
    use communicate
    use mathlib
!    use variables,only : Lrealdim,Rrealdim,Lrealdimp,Rrealdimp,opmethod,norbs,nleft,nright,orbid1,logic_PPP
    use variables
    use module_sparse

contains
!======================================================================
!======================================================================

Subroutine GetHDiag(Hdiagnosymm,nbasis,&
cap_big,cap_bigcol,cap_bigrow,&
cap_Hbig,cap_Hbigcol,cap_Hbigrow,&
cap_goodbasis,ifperturbation)

    use ABop
    
    implicit none
    
    integer,intent(in) :: nbasis
    real(kind=r8),intent(out) :: Hdiagnosymm(:)
    real(kind=r8),intent(in) :: cap_big(:,:),cap_Hbig(:,:)
    integer(kind=i4),intent(in) :: &
        cap_bigcol(:,:),cap_bigrow(:,:),&
        cap_Hbigcol(:,:),cap_Hbigrow(:,:),&
        cap_goodbasis(:,:)
    logical,intent(in) :: ifperturbation
    
    ! local
    integer :: iLrealdim,iRrealdim,ierr
    real(kind=r8) :: starttime,endtime
    
    call master_print_message("enter in GetHDiag subroutine")
    starttime=MPI_WTIME()
    
    if(myid==0) Hdiagnosymm=0.0D0

    if(ifperturbation==.true.) then
        iLrealdim=Lrealdimp
        iRrealdim=Rrealdimp
    else
        iLrealdim=Lrealdim
        iRrealdim=Rrealdim
    end if

    if(myid==0) then   
        call HDiagHLRcontribute(iLrealdim,iRrealdim,nbasis,&
                cap_Hbig,cap_Hbigcol,cap_Hbigrow,cap_goodbasis,Hdiagnosymm)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(opmethod=="comple") then
        call Com_HdiagNoSymm(iLrealdim,iRrealdim,nbasis,&
                cap_big,cap_bigcol,cap_bigrow,cap_goodbasis,&
                pppVCommat(:,1,:),pppVComcol(:,1,:),&
                pppVComrow(:,1,:),pppVComid,pppVComoperanum,Hdiagnosymm)
    else if(opmethod=="direct") then
        call Direct_HdiagNosymm(iLrealdim,iRrealdim,nbasis,&
                cap_big,cap_bigcol,cap_bigrow,cap_goodbasis,Hdiagnosymm)
    end if
    endtime=MPI_WTIME()
    call master_print_message(endtime-starttime,"HDIAGTIME:")
    return

end subroutine GetHDiag

!==============================================================
!==============================================================

subroutine Direct_HdiagNosymm(iLrealdim,iRrealdim,nbasis,&
                cap_big,cap_bigcol,cap_bigrow,cap_goodbasis,Hdiagnosymm)
    
    use CoeffTrans
    
    implicit none
    integer(kind=i4),intent(in) :: iLrealdim,iRrealdim,nbasis,&
        cap_bigcol(:,:),cap_bigrow(:,:),cap_goodbasis(:,:)
    real(kind=r8),intent(in) :: cap_big(:,:)
    real(kind=r8),intent(inout) :: Hdiagnosymm(:)
    
    ! local
    real(kind=r8),allocatable :: buffermat(:),buffermat0(:,:),Hdiagdummy(:),localdiag(:)
    real(kind=r8) :: alpha
    integer :: maxdim,operaindex
    integer :: i,j,k,m
    ! MPI_FLAG
    integer :: status(MPI_STATUS_SIZE),ierr
    
    maxdim=max(iLrealdim,iRrealdim)
    
    if(myid/=0) then
        allocate(buffermat(4*maxdim))
        buffermat=0.0D0
    else
        allocate(buffermat0(4*maxdim,norbs))
        buffermat0=0.0D0
        allocate(Hdiagdummy(16*iLrealdim*iRrealdim))
        Hdiagdummy=0.0D0
        allocate(localdiag(nbasis))
    end if

    ! L space
    do i=1,nleft+1,1
        if(myid==orbid1(i,1)) then
            operaindex=orbid1(i,2)*3
            ! copy the diagonal element to the buffermat
            do j=1,4*iLrealdim,1
                call SpMatIJ(4*iLrealdim,j,j,cap_big(:,operaindex),&
                        cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),buffermat(j))
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
            operaindex=orbid1(j,2)*3
            ! copy the diagonal element to the buffermat
            do i=1,4*iRrealdim,1
                call SpMatIJ(4*iRrealdim,i,i,cap_big(:,operaindex),&
                        cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),buffermat(i))
            end do
            ! send the diagonal element to the 0 process
            call MPI_SEND(buffermat,4*maxdim,mpi_real8,0,j,MPI_COMM_WORLD,ierr)
        else if(myid==0) then
            call MPI_RECV(buffermat0(1,j),4*maxdim,mpi_real8,orbid1(j,1),j,MPI_COMM_WORLD,status,ierr)
        end if
    end do

    if(myid==0) then   
        ! PPP term contribution
        if(logic_PPP==1) then
            do i=1,nleft+1,1
            do j=norbs,norbs-nright,-1
                do k=1,4*iRrealdim,1
                   alpha=buffermat0(k,j)*pppV(i,j)
                   call axpy(buffermat0(1:4*iLrealdim,i),Hdiagdummy((k-1)*4*iLrealdim+1:k*4*iLrealdim),alpha)
                end do
            end do
            end do
        end if
        call Dense16LRtoNgood(iLrealdim,iRrealdim,nbasis,cap_goodbasis,Hdiagdummy,localdiag)
        call dxpy(nbasis,localdiag,1,Hdiagnosymm,1)               
    end if

    if(myid==0) then
        deallocate(buffermat0)
        deallocate(Hdiagdummy)
        deallocate(localdiag)
    else
        deallocate(buffermat)
    end if

return

end subroutine Direct_HdiagNosymm

!======================================================================================    
!======================================================================================    

subroutine HDiagHLRcontribute(iLrealdim,iRrealdim,nbasis,&
                cap_Hbig,cap_Hbigcol,cap_Hbigrow,cap_goodbasis,Hdiagnosymm)
    
    use CoeffTrans
    implicit none
    integer(kind=i4),intent(in) :: iLrealdim,iRrealdim,nbasis,&
        cap_Hbigcol(:,:),cap_Hbigrow(:,:),cap_goodbasis(:,:)
    real(kind=r8),intent(in) :: cap_Hbig(:,:)
    real(kind=r8),intent(inout) :: Hdiagnosymm(:)
    ! local
    real(kind=r8),allocatable :: Houtput(:),Hdiagdummy(:),localdiag(:)
    integer :: i,j

    allocate(Houtput(4*iLrealdim))
    allocate(Hdiagdummy(16*iLrealdim*iRrealdim))
    allocate(localdiag(nbasis))
    
    ! HL contribution
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
            cap_Hbigcol(:,2),cap_Hbigrow(:,2),Houtput(1))
        call dapy(4*iLrealdim,Houtput(1),Hdiagdummy((i-1)*4*iLrealdim+1:i*4*iLrealdim),1)
    end do

    call Dense16LRtoNgood(iLrealdim,iRrealdim,nbasis,cap_goodbasis,Hdiagdummy,localdiag)
    call dxpy(nbasis,localdiag,1,Hdiagnosymm,1)               
    
    deallocate(Houtput,Hdiagdummy,localdiag)

end subroutine HDiagHLRcontribute

!======================================================================================    
!====2==================================================================================    

subroutine Com_HdiagNoSymm(iLrealdim,iRrealdim,nbasis,&
                cap_big,cap_bigcol,cap_bigrow,cap_goodbasis,&
                Commat,Comcol,Comrow,Comid,Comoperanum,Hdiagnosymm)
    implicit none
    integer(kind=i4),intent(in) :: iLrealdim,iRrealdim,nbasis,&
        cap_bigcol(:,:),cap_bigrow(:,:),cap_goodbasis(:,:),&
        Comcol(:,:),Comrow(:,:),Comid(:,:),Comoperanum(:)
    real(kind=r8),intent(in) :: cap_big(:,:),Commat(:,:)
    real(kind=r8),intent(inout) :: Hdiagnosymm(:)
    
    ! local
    real(kind=r8),allocatable :: localdiag(:)
    integer :: iproc
    ! MPI_FLAG
    integer :: status(MPI_STATUS_SIZE),ierr

    if(myid==0 .or. Comoperanum(myid)/=0) then
        allocate(localdiag(nbasis))
    end if

    do iproc=1,nprocs-1,1
        if(Comoperanum(iproc)/=0) then
            if(myid==iproc) then
                call HdiagPPPVcontribute(iLrealdim,iRrealdim,&
                            nbasis,cap_big,cap_bigcol,cap_bigrow,&
                            Commat,Comcol,Comrow,Comid,cap_goodbasis,localdiag)
                call MPI_SEND(localdiag,nbasis,MPI_REAL8,0,myid,MPI_COMM_WORLD,ierr)
            else if(myid==0) then
                call MPI_RECV(localdiag,nbasis,MPI_REAL8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                call dxpy(nbasis,localdiag,1,Hdiagnosymm,1)               
            end if
        end if
    end do

    if(allocated(localdiag)) deallocate(localdiag)
    
    return

end subroutine Com_HdiagNoSymm

!======================================================================================    
!======================================================================================    
subroutine HdiagPPPVcontribute(iLrealdim,iRrealdim,nbasis,&
        cap_big,cap_bigcol,cap_bigrow,&
        Commat,Comcol,Comrow,Comid,cap_goodbasis,localdiag)
    
    USE BLAS95
    USE F95_PRECISION
    use CoeffTrans
    implicit none
    integer(kind=i4),intent(in) :: iLrealdim,iRrealdim,nbasis,&
        cap_bigcol(:,:),cap_bigrow(:,:),&
        Comcol(:,:),Comrow(:,:),Comid(:,:),cap_goodbasis(:,:)
    real(kind=r8),intent(in) :: cap_big(:,:),Commat(:,:)
    real(kind=r8),intent(out) :: localdiag(:)

    ! local
    real(kind=r8),allocatable :: Ldiag(:),Rdiag(:)
    integer :: ileft,il,ir,ibasis
    integer :: Loperaindex,Roperaindex
    real(kind=r8),allocatable :: Hdiagdummy(:)
    
    allocate(Ldiag(4*iLrealdim))
    allocate(Rdiag(4*iRrealdim))
    allocate(Hdiagdummy(16*iLrealdim*iRrealdim))
    Hdiagdummy=0.0D0

    do ileft=1,nleft+1,1
        if(myid==Comid(ileft,1)) then
            Loperaindex=orbid1(ileft,2)*3
            Roperaindex=Comid(ileft,2)
            do il=1,4*iLrealdim,1
                call SpMatIJ(4*iLrealdim,il,il,cap_big(:,Loperaindex),&
                    cap_bigcol(:,Loperaindex),cap_bigrow(:,Loperaindex),Ldiag(il))
            end do
            do ir=1,4*iRrealdim,1
                call SpMatIJ(4*iRrealdim,ir,ir,Commat(:,Roperaindex),&
                    Comcol(:,Roperaindex),Comrow(:,Roperaindex),Rdiag(ir))
                call axpy(Ldiag,Hdiagdummy((ir-1)*4*iLrealdim+1:ir*4*iLrealdim),Rdiag(ir))
            end do
        end if
    end do

    call Dense16LRtoNgood(iLrealdim,iRrealdim,nbasis,cap_goodbasis,Hdiagdummy,localdiag)

    deallocate(Ldiag,Rdiag,Hdiagdummy)
    return

end subroutine HdiagPPPVcontribute

!======================================================================================    
!======================================================================================    

end module GetHdiag_mod
