module OpExpec_mod
    use kinds_mod
    use communicate
    use mathlib
    USE MPI
    implicit none
contains
!================================================
!================================================

subroutine LinkOpExpec(expec0,ltrans,rtrans,lproc,rproc,leadproc,&
            iLrealdim,iRrealdim,isubM,lstate,rstate,loperaindex,roperaindex,&
            cap_big,cap_bigcol,cap_bigrow,cap_coeff,cap_coeffcol,cap_coeffrow,packratio,midratio,ifphase,phase)

    use SpMatTrans_mod
! this subroutine calculate the expectation value of opertors between L and R space
    implicit none
    
    character(len=1),intent(in) :: ltrans,rtrans
    integer(kind=i4),intent(in) :: lproc,rproc,leadproc,&
        iLrealdim,iRrealdim,isubM,&
        lstate,rstate,loperaindex,roperaindex
    real(kind=r8),intent(in) :: midratio,packratio
    real(kind=r8),intent(in) :: cap_big(:,:),cap_coeff(:,:),phase(:)
    integer(kind=i4),intent(in) :: cap_bigcol(:,:),cap_bigrow(:,:),&
                                cap_coeffcol(:,:),cap_coeffrow(:,:)
    
    logical,intent(in) :: ifphase
    real(kind=r8),intent(out) :: expec0

    ! local
    real(kind=r8),allocatable :: midmat1(:),midmat2(:)
    integer(kind=i4),allocatable :: midcol1(:),midrow1(:),&
                                    midcol2(:),midrow2(:)
    integer :: nmid,packsize,position1
    character(len=1),allocatable :: packbuf(:)
    integer :: operaindex,slaverproc,realdim
    character(len=1) :: ltranstmp
    real(kind=r8) :: expec
    integer :: i
    ! MPI 
    integer :: ierr,status(MPI_STATUS_SIZE)

    nmid=CEILING(DBLE(16*isubM*isubM)/midratio)

    if(myid==leadproc) then
        allocate(midmat1(nmid))
        allocate(midcol1(nmid))
        allocate(midrow1(4*isubM+1))
        allocate(midmat2(nmid))
        allocate(midcol2(nmid))
        allocate(midrow2(4*isubM+1))
    end if
    
    ! transter slaver matrix to midmat1
    if(myid/=0) then
        if(lproc/=rproc) then
            
            packsize=16*isubM*isubM*12/packratio+(4*isubM+1)*4
            allocate(packbuf(packsize)) ! packbuf to send local operator
            
            if(leadproc==lproc) then
                operaindex=roperaindex
                realdim=iRrealdim
                slaverproc=rproc
            else
                operaindex=loperaindex
                realdim=iLrealdim
                slaverproc=lproc
            end if

            if(myid==slaverproc) then
                position1=0
                call SpMatPack(cap_big(:,operaindex),cap_bigcol(:,operaindex),&
                        cap_bigrow(:,operaindex),4*realdim,position1,packbuf,packsize)
                call MPI_SEND(packbuf,position1,MPI_PACKED,leadproc,slaverproc,MPI_COMM_WORLD,ierr)
            else if(myid==leadproc) then
                call MPI_RECV(packbuf,packsize,MPI_PACKED,slaverproc,slaverproc,MPI_COMM_WORLD,status,ierr)
                position1=0
                call SpMatUnPack(midmat1(:),midcol1(:),midrow1(:),&
                    4*realdim,nmid,position1,packbuf,packsize)
            end if

            deallocate(packbuf)
        else
            ! if in the same process, copy the right operator to midcol
            call CopySpAtoB(4*iRrealdim,cap_big(:,roperaindex),cap_bigcol(:,roperaindex),&
                    cap_bigrow(:,roperaindex),midmat1(:),midcol1(:),midrow1(:),nmid)
        end if
    end if

    if(myid==leadproc) then
        if(myid==lproc) then
            call SpMMtoSp('N',rtrans,4*iLrealdim,4*iRrealdim,4*iRrealdim,4*iRrealdim,4*iLrealdim,&
                    cap_coeff(:,lstate),cap_coeffcol(:,lstate),cap_coeffrow(:,lstate),&
                    midmat1(:),midcol1(:),midrow1(:),midmat2(:),midcol2(:),midrow2(:),nmid)
            if(ltrans=="T") then
                ltranstmp="N"
            else
                ltranstmp="T"
            end if
            call SpMMtoSp(ltranstmp,'N',4*iLrealdim,4*iLrealdim,4*iLrealdim,4*iRrealdim,4*iLrealdim,&
                    cap_big(:,loperaindex),cap_bigcol(:,loperaindex),cap_bigrow(:,loperaindex),&
                    midmat2(:),midcol2(:),midrow2(:),midmat1(:),midcol1(:),midrow1(:),nmid)
        else if(myid==rproc) then
            if(ltrans=="T") then
                ltranstmp="N"
            else
                ltranstmp="T"
            end if
            call SpMMtoSp(ltranstmp,'N',4*iLrealdim,4*iLrealdim,4*iLrealdim,4*iRrealdim,4*iLrealdim,&
                    midmat1(:),midcol1(:),midrow1(:),&
                    cap_coeff(:,lstate),cap_coeffcol(:,lstate),cap_coeffrow(:,lstate),&
                    midmat2(:),midcol2(:),midrow2(:),nmid)
            call SpMMtoSp('N',rtrans,4*iLrealdim,4*iRrealdim,4*iRrealdim,4*iRrealdim,4*iLrealdim,&
                    midmat2(:),midcol2(:),midrow2(:),&
                    cap_big(:,roperaindex),cap_bigcol(:,roperaindex),cap_bigrow(:,roperaindex),&
                    midmat1(:),midcol1(:),midrow1(:),nmid)
        end if

        call CopySpAtoB(4*iLrealdim,cap_coeff(:,rstate),cap_coeffcol(:,rstate),cap_coeffrow(:,rstate),&
                midmat2(:),midcol2(:),midrow2(:),nmid)
        
        if(ifphase==.true.) then
            do i=1,4*iLrealdim,1
                midmat2(midrow2(i):midrow2(i+1)-1)=midmat2(midrow2(i):midrow2(i+1)-1)*phase(i)
            end do
        end if

        call SpMMtrace('T',4*iLrealdim,4*iRrealdim,midmat1(:),midcol1(:),midrow1(:),&
                4*iLrealdim,4*iRrealdim,&
                midmat2(:),midcol2(:),midrow2(:),expec)
        
        call MPI_SEND(expec,1,MPI_REAL8,0,0,MPI_COMM_WORLD,ierr)

    else if(myid==0) then
        expec0=0.0D0
        call MPI_RECV(expec0,1,MPI_REAL8,leadproc,0,MPI_COMM_WORLD,status,ierr)
    end if

    if(myid==leadproc) then
        deallocate(midmat1,midcol1,midrow1)
        deallocate(midmat2,midcol2,midrow2)
    end if
    
    return

end subroutine LinkOpExpec

!================================================
!================================================

subroutine SubSpaceOpExpec(expec0,iproc,domain,&
            iLrealdim,iRrealdim,isubM,lstate,rstate,operaindex,&
            cap_big,cap_bigcol,cap_bigrow,cap_coeff,cap_coeffcol,cap_coeffrow,midratio)

! this subroutine calculate the expectation value of opertors in L/R space
    implicit none
    
    character(len=1),intent(in) :: domain
    integer(kind=i4),intent(in) :: iproc,&
        iLrealdim,iRrealdim,isubM,&
        lstate,rstate,operaindex
    real(kind=r8),intent(in) :: midratio
    real(kind=r8),intent(in) :: cap_big(:,:),cap_coeff(:,:)
    integer(kind=i4),intent(in) :: cap_bigcol(:,:),cap_bigrow(:,:),&
                                cap_coeffcol(:,:),cap_coeffrow(:,:)
    real(kind=r8),intent(out) :: expec0
    
    ! local
    real(kind=r8),allocatable :: midmat(:)
    integer(kind=i4),allocatable :: midcol(:),midrow(:)
    integer :: nmid
    real(kind=r8) :: expec
    ! MPI 
    integer :: ierr,status(MPI_STATUS_SIZE)

    nmid=CEILING(DBLE(16*isubM*isubM)/midratio)
    
    if(myid==iproc) then
        
        allocate(midmat(nmid))
        allocate(midcol(nmid))
        allocate(midrow(4*isubM+1))
        
        if(domain=="L") then
            call SpMMtoSp('T','N',4*iLrealdim,4*iLrealdim,4*iLrealdim,4*iRrealdim,4*iRrealdim,&
                    cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),&
                    cap_coeff(:,lstate),cap_coeffcol(:,lstate),cap_coeffrow(:,lstate),&
                    midmat(:),midcol(:),midrow(:),nmid)
        else if(domain=="R") then
            call SpMMtoSp('N','N',4*iLrealdim,4*iRrealdim,4*iRrealdim,4*iRrealdim,4*iRrealdim,&
                    cap_coeff(:,lstate),cap_coeffcol(:,lstate),cap_coeffrow(:,lstate),&
                    cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),&
                    midmat(:),midcol(:),midrow(:),nmid)
        end if
        call SpMMtrace('T',4*iLrealdim,4*iRrealdim,midmat(:),midcol(:),midrow(:),&
                4*iLrealdim,4*iRrealdim,&
                cap_coeff(:,rstate),cap_coeffcol(:,rstate),cap_coeffrow(:,rstate),expec)
        call MPI_SEND(expec,1,MPI_REAL8,0,0,MPI_COMM_WORLD,ierr)

        deallocate(midmat,midcol,midrow)
    else if(myid==0) then
        expec0=0.0D0
        call MPI_RECV(expec0,1,MPI_REAL8,iproc,0,MPI_COMM_WORLD,status,ierr)
    end if

    return
    
end subroutine SubSpaceOpExpec   
    
!================================================
!================================================
end module OpExpec_mod
