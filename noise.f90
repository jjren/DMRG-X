module noise_mod
    use variables
    use module_sparse
    use communicate
    use mpi


    implicit none

    logical :: Ifnoise=.false.
    real(kind=r8),allocatable :: noiseweight(:)

contains
!=====================================================================
!=====================================================================

subroutine setvalue_noisemod
    implicit none
    integer :: i

    allocate(noiseweight(0:sweeps))
    do i=0,sweeps,1
        select case(i)
        case(0:1)
            noiseweight(i)=1.0D-1
        case(2:3)
            noiseweight(i)=1.0D-2
        case default
            noiseweight(i)=0.0D0
        end select
    end do
    return
end subroutine setvalue_noisemod

!=====================================================================
!=====================================================================

subroutine svd_noise (iLrealdim,iRrealdim,LRcoeff)
    
    implicit none

    integer,intent(in) :: iLrealdim,iRrealdim
    real(kind=r8),intent(inout) :: LRcoeff(:,:)

    call master_print_message("enter svd_noise subroutine")
    
    call setvalue_noisemod
    if(ifopenperturbation==.true.) then
        call svd_noise_wrapper(iLrealdim,iRrealdim,LRcoeff,&
        operamatbig1p,bigcolindex1p,bigrowindex1p,quantabigLp,quantabigRp,.true.)
    else
        call svd_noise_wrapper(iLrealdim,iRrealdim,LRcoeff,&
        operamatbig1,bigcolindex1,bigrowindex1,quantabigL,quantabigR,.false.)
    end if
    
    deallocate(noiseweight)
    return

end subroutine svd_noise

!=====================================================================
!=====================================================================

subroutine svd_noise_wrapper(iLrealdim,iRrealdim,LRcoeff,&
                    cap_big,cap_bigcol,cap_bigrow,&
                    cap_quantabigL,cap_quantabigR,ifperturbation)
    use BLAS95
    use F95_PRECISION
    implicit none
    integer,intent(in) :: iLrealdim,iRrealdim
    real(kind=r8),intent(inout) :: LRcoeff(:,:)
    real(kind=r8),intent(in) :: cap_big(:,:)
    integer(kind=i4),intent(in) :: cap_bigcol(:,:),cap_bigrow(:,:),&
                    cap_quantabigL(:,:),cap_quantabigR(:,:)
    logical :: ifperturbation
    ! local

    character(len=1),allocatable :: packbuf(:)
    integer :: position1,packsize
    real(kind=r8),allocatable :: leftmat(:,:),rightmat(:,:),midmat(:,:),dummymidmat(:,:),LRcoeffcorrect(:,:)
    integer(kind=i4),allocatable :: leftcol(:,:),leftrow(:,:),rightcol(:,:),rightrow(:,:)
    character(len=1) :: operation(6)
    integer :: ierr
    integer :: i,j,k,itrans
    integer :: ibigdim1,isubM,operaindex,nonzero,dim1
    character(len=1) :: trans1,trans2
    integer :: status(MPI_STATUS_SIZE)
    real(kind=r8) :: norm

    if(ifperturbation==.true.) then
        ibigdim1=bigdim1p
        isubM=subMp
    else
        ibigdim1=bigdim1
        isubM=subM
    end if
    packsize=(ibigdim1*12+(4*isubM+1)*4)*2

    if(myid==orbid1(norbs-nright,1) .or. myid==orbid1(nleft+1,1) .or. myid==0) then
        allocate(packbuf(packsize))
    end if
    if(myid==0) then
        allocate(leftmat(ibigdim1,2))
        allocate(leftcol(ibigdim1,2))
        allocate(leftrow(4*isubM+1,2))
        allocate(rightmat(ibigdim1,2))
        allocate(rightcol(ibigdim1,2))
        allocate(rightrow(4*isubM+1,2))
    end if
    
!   if(myid==0) then
!       norm=0.0D0
!       k=0
!       do i=1,4*iRrealdim,1
!       do j=1,4*iLrealdim,1
!           if((cap_quantabigL(j,1)+cap_quantabigR(i,1)==nelecs) .and. &
!               cap_quantabigL(j,2)+cap_quantabigR(i,2)==totalSz) then
!               norm=norm+LRcoeff(j,i)*LRcoeff(j,i)
!               k=k+1
!           end if
!       end do      
!       end do      
!       write(*,*) "noise norm1=",norm
!       write(*,*) "goodquantumnumber noise :" ,k
!   end if

    ! leftsend
    if(myid==orbid1(nleft+1,1)) then
        position1=0
        do k=1,2,1
            operaindex=orbid1(nleft+1,2)*3-3+k
            call MPI_PACK(cap_bigrow(1,operaindex),(4*iLrealdim+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            nonzero=cap_bigrow(4*iLrealdim+1,operaindex)-1
            call MPI_PACK(cap_big(1,operaindex),nonzero,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            call MPI_PACK(cap_bigcol(1,operaindex),nonzero,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
        end do
        call MPI_SEND(packbuf,position1,MPI_PACKED,0,nleft+1,MPI_COMM_WORLD,ierr)
    end if
    
    if(myid==0) then
        call MPI_RECV(packbuf,packsize,MPI_PACKED,orbid1(nleft+1,1),nleft+1,MPI_COMM_WORLD,status,ierr)
        position1=0
        do k=1,2,1
            call MPI_UNPACK(packbuf,packsize,position1,leftrow(1,k),(4*iLrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
            nonzero=leftrow(4*iLrealdim+1,k)-1
            call MPI_UNPACK(packbuf,packsize,position1,leftmat(1,k),nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position1,leftcol(1,k),nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
        end do
    end if

    ! rightsend
    if(myid==orbid1(norbs-nright,1)) then
        position1=0
        do k=1,2,1
            operaindex=orbid1(norbs-nright,2)*3-3+k
            call MPI_PACK(cap_bigrow(1,operaindex),(4*iRrealdim+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            nonzero=cap_bigrow(4*iRrealdim+1,operaindex)-1
            call MPI_PACK(cap_big(1,operaindex),nonzero,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            call MPI_PACK(cap_bigcol(1,operaindex),nonzero,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
        end do
        call MPI_SEND(packbuf,position1,MPI_PACKED,0,norbs-nright,MPI_COMM_WORLD,ierr)
    end if
    
    if(myid==0) then
        call MPI_RECV(packbuf,packsize,MPI_PACKED,orbid1(norbs-nright,1),norbs-nright,MPI_COMM_WORLD,status,ierr)
        position1=0
        do k=1,2,1
            call MPI_UNPACK(packbuf,packsize,position1,rightrow(1,k),(4*iRrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
            nonzero=rightrow(4*iRrealdim+1,k)-1
            call MPI_UNPACK(packbuf,packsize,position1,rightmat(1,k),nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position1,rightcol(1,k),nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
        end do
    end if
    
    if(myid==orbid1(norbs-nright,1) .or. myid==orbid1(nleft+1,1) .or. myid==0) then
        deallocate(packbuf)
    end if

    operation(1)='G'
    operation(2)='L'
    operation(3)='N'
    operation(4)='F'
    
    if(myid==0) then
        dim1=max(4*iLrealdim,4*iRrealdim)
        allocate(midmat(dim1,dim1))
        allocate(dummymidmat(dim1,dim1))
        allocate(LRcoeffcorrect(dim1,dim1))
        
        do i=1,2,1
            if(i==1) then
                trans1='N'
                trans2='N'
            else if(i==2) then
                trans1='T'
                trans2='T'
            end if
            do k=1,2,1
                midmat=0.0D0
                LRcoeffcorrect=0.0D0
                call mkl_dcsrmm(trans1,4*iLrealdim,4*iRrealdim,4*iLrealdim,1.0D0,operation,leftmat(:,k),&
                    leftcol(:,k),leftrow(1:4*iLrealdim,k),leftrow(2:4*iLrealdim+1,k),&
                    LRcoeff,4*iLrealdim,0.0D0,midmat,dim1)
                !midmat=transpose(midmat)
                dummymidmat=midmat
                do itrans=1,dim1,1
                    call copy(dummymidmat(:,itrans),midmat(itrans,:))
                end do

                call mkl_dcsrmm(trans2,4*iRrealdim,4*iLrealdim,4*iRrealdim,1.0D0,operation,rightmat(:,k),&
                    rightcol(:,k),rightrow(1:4*iRrealdim,k),rightrow(2:4*iRrealdim+1,k),&
                    midmat,dim1,0.0D0,LRcoeffcorrect,dim1)
                !LRcoeffcorrect=transpose(LRcoeffcorrect)
                dummymidmat=LRcoeffcorrect
                do itrans=1,dim1,1
                    call copy(dummymidmat(:,itrans),LRcoeffcorrect(itrans,:))
                end do

            !   if(i==2) then
                LRcoeff=LRcoeff+LRcoeffcorrect(1:4*iLrealdim,1:4*iRrealdim)*noiseweight(isweep)
            end do
        end do
        
        norm=0.0D0
        do i=1,4*iRrealdim,1
        do j=1,4*iLrealdim,1
            if((cap_quantabigL(j,1)+cap_quantabigR(i,1)==nelecs) .and. &
                cap_quantabigL(j,2)+cap_quantabigR(i,2)==totalSz) then
                norm=norm+LRcoeff(j,i)*LRcoeff(j,i)
            else
            !   write(*,*) LRcoeff(j,i)
                LRcoeff(j,i)=0.0D0
            end if
        end do      
        end do      
        write(*,*) "noise norm=",norm
        norm=sqrt(norm)
        LRcoeff=LRcoeff/norm

        deallocate(midmat,dummymidmat,LRcoeffcorrect)
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    if(myid==0) then
        deallocate(leftmat,leftcol,leftrow)
        deallocate(rightmat,rightcol,rightrow)
    end if

end subroutine svd_noise_wrapper 

!=====================================================================
!=====================================================================
end module noise_mod
