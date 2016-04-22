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
    integer :: lsvddim,rsvddim,svdvaluedim

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
    integer :: iLrealdim,iRrealdim,isubM
    real(kind=r8),allocatable :: LRcoeff(:,:)
    integer :: isvd
    real(kind=r8),allocatable :: localweight(:)
    call master_print_message("enter Renormalization subroutine")

    if(ifopenperturbation==.true.) then
        iLrealdim=Lrealdimp
        iRrealdim=Rrealdimp
        isubM=subMp
    else
        iLrealdim=Lrealdim
        iRrealdim=Rrealdim
        isubM=subM
    end if
    lsvddim=min(4*iLrealdim,isubM)
    rsvddim=min(4*iRrealdim,isubM)
    svdvaluedim=min(lsvddim,rsvddim)


    if(4*Lrealdim>subM .or. 4*Rrealdim>subM) then  ! do Renormalization in this case
        
        if(myid==0) allocate(LRcoeff(4*iLrealdim,4*iRrealdim))
        
        call Prepare_LRcoeff(iLrealdim,iRrealdim,LRcoeff,1)
        
        if(myid==0) then
            ! allocate work array
            allocate(singularvalue(svdvaluedim))
            allocate(leftu(4*iLrealdim,lsvddim))
            allocate(rightv(rsvddim,4*iRrealdim))
            
            ! get rotate matrix leftu/rightv and singularvalue
            if(nstate==1 .or. exscheme==4) then
                ! the index 1 can be changed!
                if(ifopenperturbation==.false.) then
                    call splitsvd_direct(iLrealdim,iRrealdim,LRcoeff,lsvddim,rsvddim,svdvaluedim,leftu,rightv,singularvalue,&
                        quantabigL(1:4*iLrealdim,1:2),quantabigR(1:4*iRrealdim,1:2),&
                        quantasmaL(1:lsvddim,1:2),quantasmaR(1:rsvddim,1:2))
                else
                    call splitsvd_direct(iLrealdim,iRrealdim,LRcoeff,lsvddim,rsvddim,svdvaluedim,leftu,rightv,singularvalue,&
                        quantabigLp(1:4*iLrealdim,1:2),quantabigRp(1:4*iRrealdim,1:2),&
                        quantasmaLp(1:lsvddim,1:2),quantasmaRp(1:rsvddim,1:2))
                    write(*,*) "get lsvddim=",lsvddim
                    write(*,*) "get rsvddim=",rsvddim
                    if(lsvddim<subM .or. rsvddim<subM) then
                        call master_print_message("lsvddim/rsvddim<subM")
                    end if
                end if
            else if(exscheme==1) then
                if(C2method=="mix" .and. nleft==nright .and. direction/="i") then
                    allocate(localweight(C2state))
                    localweight=1.0D0/C2state
                    call splitsvd('L',iLrealdim,1,C2state,localweight)
                    call splitsvd('R',iRrealdim,1,C2state,localweight)
                    deallocate(localweight)
                else 
                    call splitsvd('L',iLrealdim,1,nstate,nweight)
                    call splitsvd('R',iRrealdim,1,nstate,nweight)
                end if
            else if(exscheme==2) then
                ! my new exScheme
                call ExScheme2
            end if
            write(*,*) "singular value:"
            do isvd=1,5,1
                write(*,*) singularvalue(isvd)
            end do
            
            ! store the wavefunction
            call StoreWaveFunction
        end if

        if(ifopenperturbation==.false.) then
            ! bcast the updated quantasmaL/R because only the 0 process know it
            call MPI_BCAST(quantasmaL(1,1),subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(quantasmaR(1,1),subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
        else
            call MPI_BCAST(lsvddim,1,MPI_integer4,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(rsvddim,1,MPI_integer4,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(quantasmaLp(1,1),subMp*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(quantasmaRp(1,1),subMp*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
            quantasmaL(:,:)=quantasmaLp(1:subM,:)
            quantasmaR(:,:)=quantasmaRp(1:subM,:)
        end if

        ! rotate the bigmat to smamat
        ! L subspace
        if(direction=='i' .or. direction=='l') then
            if(ifopenperturbation==.false.)  then
            !   if(myid==0) then
            !       write(*,*) "leftu1"
            !       write(*,*) leftu
            !   end if
                call RotateBasis('L',operamatbig1,bigcolindex1,bigrowindex1,operamatsma1,smacolindex1,smarowindex1,&
                Hbig,Hbigcolindex,Hbigrowindex,Hsma,Hsmacolindex,Hsmarowindex,Lrealdim,subM,Hsmadim,smadim1)
            !   if(myid==0) then
            !       write(*,*) "Hsma"
            !       write(*,*) Hsma(:,1)
            !       write(*,*) Hsmacolindex(:,1)
            !       write(*,*) Hsmarowindex(:,1)
            !   end if
            else
            !   if(myid==0) then
            !       write(*,*) "leftu1"
            !       write(*,*) leftu
            !   end if
                call RotateBasis('L',operamatbig1p,bigcolindex1p,bigrowindex1p,operamatsma1p,smacolindex1p,smarowindex1p,&
                    Hbigp,Hbigcolindexp,Hbigrowindexp,Hsmap,Hsmacolindexp,Hsmarowindexp,Lrealdimp,lsvddim,Hsmadimp,smadim1p)
                call RotateBasis('L',operamatbig1p,bigcolindex1p,bigrowindex1p,operamatsma1,smacolindex1,smarowindex1,&
                    Hbigp,Hbigcolindexp,Hbigrowindexp,Hsma,Hsmacolindex,Hsmarowindex,Lrealdimp,subM,Hsmadim,smadim1)
                ! update the Lrealdimp value
                Lrealdimp=lsvddim
            !   if(myid==0) then
            !       write(*,*) "Hsma"
            !       write(*,*) Hsma(:,1)
            !       write(*,*) Hsmacolindex(:,1)
            !       write(*,*) Hsmarowindex(:,1)
            !   end if
            end if
        end if
        ! R subspace
        if(direction=='i' .or. direction=='r') then
            if(ifopenperturbation==.false.)  then
                call RotateBasis('R',operamatbig1,bigcolindex1,bigrowindex1,operamatsma1,smacolindex1,smarowindex1,&
                    Hbig,Hbigcolindex,Hbigrowindex,Hsma,Hsmacolindex,Hsmarowindex,Rrealdim,subM,Hsmadim,smadim1)
            else
                call RotateBasis('R',operamatbig1p,bigcolindex1p,bigrowindex1p,operamatsma1p,smacolindex1p,smarowindex1p,&
                    Hbigp,Hbigcolindexp,Hbigrowindexp,Hsmap,Hsmacolindexp,Hsmarowindexp,Rrealdimp,rsvddim,Hsmadimp,smadim1p)
                call RotateBasis('R',operamatbig1p,bigcolindex1p,bigrowindex1p,operamatsma1,smacolindex1,smarowindex1,&
                    Hbigp,Hbigcolindexp,Hbigrowindexp,Hsma,Hsmacolindex,Hsmarowindex,Rrealdimp,subM,Hsmadim,smadim1)
                ! update the Rrealdimp value
                Rrealdimp=rsvddim
            end if
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
            deallocate(LRcoeff)
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
    !   call splitsvdL(singularvalue((i-1)*subM+1:i*subM),leftu2(:,(i-1)*subM+1:i*subM),i,i,nleft+1)
    !   call splitsvdR(singularvalue((i-1)*subM+1:i*subM),rightv2((i-1)*subM+1:i*subM,:),i,i,norbs-nright)
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
    integer :: iLrealdim,iRrealdim
    
    ! 4 byte as 1 direct file block
    ! matrix plus two integer number
    reclength=2*subMp*subMp+2

    if(ifopenperturbation==.true.) then
        iLrealdim=Lrealdimp
        iRrealdim=Rrealdimp
    else
        iLrealdim=Lrealdim
        iRrealdim=Rrealdim
    end if
    
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
        write(105,rec=4*nleft+i) iLrealdim,lsvddim,leftu((i-1)*iLrealdim+1:i*iLrealdim,1:lsvddim)
    end do
    ! rightv
    do i=1,4,1
        write(105,rec=4*(norbs-nright-1)+i) iRrealdim,rsvddim,rightv(1:rsvddim,i:4*iRrealdim:4)
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

subroutine RotateBasis(domain,cap_big,cap_bigcol,cap_bigrow,cap_sma,cap_smacol,cap_smarow,&
cap_Hbig,cap_Hbigcol,cap_Hbigrow,cap_Hsma,cap_Hsmacol,cap_Hsmarow,&
bigLRdim,smaLRdim,dummyHsmadim,dummysmadim)
! Rotate Basis ; In fact transfer the operator matrix
! to new basis
    
    use exit_mod
    use mathlib
    use mpi
    use checkmem_mod
    implicit none
    include "mkl_spblas.fi"

    character(len=1),intent(in) :: domain
    real(kind=r8),intent(in) :: cap_big(:,:),cap_Hbig(:,:)
    real(kind=r8),intent(out) :: cap_sma(:,:),cap_Hsma(:,:)
    integer(kind=i4),intent(in) :: cap_bigcol(:,:),cap_bigrow(:,:),&
        cap_Hbigcol(:,:),cap_Hbigrow(:,:)
    integer(kind=i4),intent(out) :: cap_smacol(:,:),cap_smarow(:,:),&
        cap_Hsmacol(:,:),cap_Hsmarow(:,:)
    integer,intent(in) :: bigLRdim,smaLRdim,dummyHsmadim,dummysmadim
    ! local
    real(kind=r8),allocatable ::     &
            rotatematdens(:,:) , &     ! store leftu and rightv
            rotatemat(:)               ! store the CSR sparse format
    integer :: orbstart,orbend,Hindex
    integer :: arraylength,nrows,operaindex
    integer(kind=i4),allocatable :: rotcolindex(:),rotrowindex(:)
    integer :: job(8)
    integer :: error,ierr,info
    integer :: i,j,k,itrans
    logical :: ifbondord,iflocalspin

    character(len=1),allocatable :: packbuf(:)
    integer :: packsize
    integer :: position1
    real(kind=r8) :: tmpratio(1)
    
    if(domain=='L') then
        orbstart=1
        orbend=nleft+1
        Hindex=1
        if(nleft<=(norbs-1)/2) then
            ifbondord=.true.
            iflocalspin=.true.
        else
            ifbondord=.false.
            iflocalspin=.false.
        end if
    else if(domain=='R') then
        orbstart=norbs-nright
        orbend=norbs
        Hindex=2
        if(nright<=(norbs-2)/2) then
            ifbondord=.true.
            iflocalspin=.true.
        else
            ifbondord=.false.
            iflocalspin=.false.
        end if
    else
        call exit_DMRG(sigAbort,"RotateBasis domain/=L/R")
    end if

    ! define the U and V nonzero element numbers
    arraylength=CEILING(DBLE(4*bigLRdim*smaLRdim)/UVmatratio)

    ! leftu rightv store in CSR format
    allocate(rotatemat(arraylength),stat=error)
    if(error/=0) stop
    allocate(rotcolindex(arraylength),stat=error)
    if(error/=0) stop
    allocate(rotrowindex(4*bigLRdim+1),stat=error) 
    if(error/=0) stop

    packsize=arraylength*12+4*(4*bigLRdim+1)+1000
    allocate(packbuf(packsize),stat=error) 
    if(error/=0) stop

    if(myid==0) then
        ! store the dense U/V
        allocate(rotatematdens(4*bigLRdim,smaLRdim),stat=error)
        if(error/=0) stop
        if(domain=='L') then
            rotatematdens=leftu(1:4*bigLRdim,1:smaLRdim)
        else if(domain=='R') then
            ! be careful about the transpose
            ! rotatematdens=transpose(rightv(1:smaLRdim,1:4*bigLRdim))
            do itrans=1,smaLRdim,1
                call copy(rightv(itrans,1:4*bigLRdim),rotatematdens(:,itrans))
            end do
        end if
        job(1)=0
        job(2)=1
        job(3)=1
        job(4)=2
        job(5)=arraylength
        job(6)=1
        nrows=4*bigLRdim
        call mkl_ddnscsr(job,nrows,smaLRdim,rotatematdens,nrows,rotatemat,rotcolindex,rotrowindex,info)
        if(info==0) then
            tmpratio(1)=DBLE(rotrowindex(4*bigLRdim+1)-1)/DBLE(4*bigLRdim*smaLRdim)
            call checkmem_OPmodMat("UVratio",tmpratio(1),1)
        else
            call master_print_message(info,"info/=0 dens to sparse in rotatebasis")
            stop
        end if
        deallocate(rotatematdens)

        position1=0
        call MPI_PACK(rotrowindex,4*bigLRdim+1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
        call MPI_PACK(rotatemat,rotrowindex(4*bigLRdim+1)-1,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
        call MPI_PACK(rotcolindex,rotrowindex(4*bigLRdim+1)-1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
    end if
    
    ! broadcast the rotate matrix
    call MPI_BCAST(position1,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(packbuf,position1,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    
    if(myid/=0) then
        position1=0
        call MPI_UNPACK(packbuf,packsize,position1,rotrowindex,4*bigLRdim+1,MPI_integer4,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position1,rotatemat,rotrowindex(4*bigLRdim+1)-1,MPI_real8,MPI_COMM_WORLD,ierr)
        call MPI_UNPACK(packbuf,packsize,position1,rotcolindex,rotrowindex(4*bigLRdim+1)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
    end if
    
    ! rotate the matrix
    do i=orbstart,orbend,1
        if(myid==orbid1(i,1)) then
            do j=1,3,1
                if(logic_PPP==0 .and. j==3) exit
                operaindex=orbid1(i,2)*3-3+j
                call SpMatRotateBasis(smaLRdim,4*bigLRdim,rotatemat,rotcolindex,rotrowindex, &
                            4*bigLRdim,4*bigLRdim,cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex), &
                            smaLRdim,cap_sma(:,operaindex),cap_smacol(:,operaindex),cap_smarow(:,operaindex),dummysmadim)
            end do
        end if
    end do

    ! rotate HL and HR
    if(myid==0) then
    !   write(*,*) cap_Hbigrow(1:4*bigLRdim+1,Hindex)
        call SpMatRotateBasis(smaLRdim,4*bigLRdim,rotatemat,rotcolindex,rotrowindex, &
                    4*bigLRdim,4*bigLRdim,cap_Hbig(:,Hindex),cap_Hbigcol(:,Hindex),cap_Hbigrow(:,Hindex), &
                    smaLRdim,cap_Hsma(:,Hindex),cap_Hsmacol(:,Hindex),cap_Hsmarow(:,Hindex),dummyHsmadim)
    !   write(*,*) cap_Hsmarow(1:smaLRdim+1,Hindex)
    end if
    
    ! rotate the bond order matrix
    if(logic_bondorder/=0 .and. ifbondord==.true.) then
        do i=orbstart,orbend,1
        do j=i,orbend,1
            if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
                if(myid==orbid2(i,j,1)) then
                    do k=1,2,1
                        operaindex=orbid2(i,j,2)*2-2+k
                        call SpMatRotateBasis(subM,4*bigLRdim,rotatemat,rotcolindex,rotrowindex, &
                                4*bigLRdim,4*bigLRdim,operamatbig2(:,operaindex),bigcolindex2(:,operaindex),bigrowindex2(:,operaindex), &
                                subM,operamatsma2(:,operaindex),smacolindex2(:,operaindex),smarowindex2(:,operaindex),smadim2)
                    end do
                end if
            end if
        end do
        end do
    end if
    ! rotate the localspin matrix
    if(logic_localspin==1 .and. iflocalspin==.true.) then
        do i=orbstart,orbend,1
        do j=i,orbend,1
            if(myid==orbid3(i,j,1)) then
                do k=1,2,1
                    operaindex=orbid3(i,j,2)-2+k
                    call SpMatRotateBasis(subM,4*bigLRdim,rotatemat,rotcolindex,rotrowindex, &
                            4*bigLRdim,4*bigLRdim,operamatbig3(:,operaindex),bigcolindex3(:,operaindex),bigrowindex3(:,operaindex), &
                            subM,operamatsma3(:,operaindex),smacolindex3(:,operaindex),smarowindex3(:,operaindex),smadim3)
                end do
            end if
        end do
        end do
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
    use blas95
    use f95_precision
    use exit_mod
    implicit none

    character(len=1) domain
    
    ! local
    integer :: i,j,k
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
                call copy(operamatbig1(1:nelement,operaindex),operamatsma1(1:nelement,operaindex))
            end do
        end if
    end do

    ! bond order matrix
    if(logic_bondorder/=0) then
        do i=orbstart,orbend,1
        do j=i,orbend,1
            if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
                if(myid==orbid2(i,j,1)) then
                    do k=1,2,1
                        operaindex=orbid2(i,j,2)*2-2+k
                        smarowindex2(:,operaindex)=0
                        smarowindex2(1:4*dim1+1,operaindex)=bigrowindex2(1:4*dim1+1,operaindex)
                        nelement=smarowindex2(4*dim1+1,operaindex)-1
                        smacolindex2(1:nelement,operaindex)=bigcolindex2(1:nelement,operaindex)
                        call copy(operamatbig2(1:nelement,operaindex),operamatsma2(1:nelement,operaindex))
                    end do
                end if
            end if
        end do
        end do
    end if

    ! local spin matrix
    if(logic_localspin==1) then
        do i=orbstart,orbend,1
        do j=i,orbend,1
            if(myid==orbid3(i,j,1)) then
                do k=1,2,1
                    operaindex=orbid3(i,j,2)-2+k
                    smarowindex3(:,operaindex)=0
                    smarowindex3(1:4*dim1+1,operaindex)=bigrowindex3(1:4*dim1+1,operaindex)
                    nelement=smarowindex3(4*dim1+1,operaindex)-1
                    smacolindex3(1:nelement,operaindex)=bigcolindex3(1:nelement,operaindex)
                    call copy(operamatbig3(1:nelement,operaindex),operamatsma3(1:nelement,operaindex))
                end do
            end if
        end do
        end do
    end if

    if(myid==0) then
        Hsmarowindex(:,Hindex)=0
        Hsmarowindex(1:4*dim1+1,Hindex)=Hbigrowindex(1:4*dim1+1,Hindex)
        nelement=Hsmarowindex(4*dim1+1,Hindex)-1
        Hsmacolindex(1:nelement,Hindex)=Hbigcolindex(1:nelement,Hindex)
        call copy(Hbig(1:nelement,Hindex),Hsma(1:nelement,Hindex))
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

subroutine splitsvd(domain,dim1,statebegin,stateend,iweight)
! this subroutine is used to split the reduced density matrix
! to different subspace according to good quantum number
! and diagonalizaiton it to get the renormalized vector
    
    use variables
    use mathlib
    use module_sparse
    USE blas95
    use lapack95
    use F95_PRECISION

    implicit none
    
    character(len=1) :: domain
    integer :: dim1,statebegin,stateend
    real(kind=r8),intent(in) :: iweight(stateend-statebegin+1)
    ! statebegin is the index of the begin state
    ! stateend is the index of the end state

    ! local
    real(kind=r8),allocatable :: &
        valuework(:)     , &     ! store eigenvalue
        coeffwork(:,:)   , &     ! workarray
        coeffbuf(:,:)    , &     ! store initial reduced density matrix
        coeffresult(:,:) , &     ! store eigenvector
        transform(:,:)   , &     ! Sz=0 basis transform matrix  ! only used in spin symmetry
        coeffdummy(:,:)          ! intemediate array in Sz=0 ! only used in spin symmetry

    integer,allocatable :: &
        valueindex(:)     , &    ! store the eigenvalue descending basis index
        szzeroindex(:)    , &    ! store the Sz=0 index
        quantabigbuf(:,:) , &    ! store the new quantabigL/R
        quantabiginit(:,:), &    ! store the old quantabigL/R
        symmlinkbigbuf(:) , &    ! store the new symmlinkbig
        subspacenum(:)           ! store the number of basis in every good quantum number subspace
    integer :: szzero,szl0,nsuborbs,Hindex
    ! szl0 means the number of sz>0 basis
    ! szzero means the number of sz=0 basis
    real(kind=r8) :: checktrace,discard
    integer :: i,j,k,l,m,n,p,q,l1,himp1
    logical :: done
    integer :: error,info

    call master_print_message("enter in subroutine splitsvd")
    
    ! store old quantabigL/R 
    allocate(quantabiginit(4*subM,2),stat=error)
    if(error/=0) stop

    if(domain=='L') then
        if(dim1/=Lrealdim) then
            call master_print_message(Lrealdim,"domain=L dim1/=Lrealdim")
            stop
        end if
        quantabiginit=quantabigL
        nsuborbs=nleft+1
        Hindex=1
    else if(domain=='R') then
        if(dim1/=Rrealdim) then
            call master_print_message(Rrealdim,"domain=R dim1/=Rrealdim")
            stop
        end if
        quantabiginit=quantabigR
        nsuborbs=nright+1
        Hindex=2
    end if

    allocate(coeffbuf(4*dim1,4*dim1),stat=error)
    if(error/=0) stop
    allocate(coeffwork(4*dim1,4*dim1),stat=error)
    if(error/=0) stop
    allocate(coeffresult(4*dim1,4*dim1),stat=error)
    if(error/=0) stop
    coeffwork=0.0D0
    coeffbuf=0.0D0
    coeffresult=0.0D0

    ! construct L/R subspace reduced density matrix
    do i=statebegin,stateend,1
        if(domain=='L') then
            call SpMMtoSptodens('N','T',4*Lrealdim,4*Rrealdim,4*Lrealdim,4*Rrealdim,4*dim1,4*dim1,&
                coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
                coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
                coeffbuf)
        !   call SpMMtoDens('N','T',4*Lrealdim,4*Rrealdim,4*Lrealdim,4*Rrealdim,4*dim1,4*dim1,&
        !       coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
        !       coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
        !       coeffbuf,ldc)
        else if(domain=='R') then
            call SpMMtoSptodens('T','N',4*Lrealdim,4*Rrealdim,4*Lrealdim,4*Rrealdim,4*dim1,4*dim1,&
                coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
                coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
                coeffbuf)
        !   call SpMMtoDens('T','N',4*Lrealdim,4*Rrealdim,4*Lrealdim,4*Rrealdim,4*dim1,4*dim1,&
        !       coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
        !       coeffIF(:,i),coeffIFcolindex(:,i),coeffIFrowindex(:,i),&
        !       coeffbuf,4*dim1)
        end if
        
        if(exscheme==1) then  ! state-average method
            coeffwork=coeffwork+coeffbuf*iweight(i)
        else  
            coeffwork=coeffbuf
        end if
    end do
    ! coeffbuf store the initial coeffwork
    coeffbuf=coeffwork

!------------------------------------------------------------------------------
    ! check if the trace == 1
    checktrace=0.0D0
    do i=1,4*dim1,1
        checktrace=checktrace+coeffwork(i,i)
    end do
    call master_print_message(checktrace,"reduced density matrix checktrace=")
!------------------------------------------------------------------------------

    allocate(valuework(4*dim1),stat=error)
    if(error/=0) stop
    valuework=0.0D0
    allocate(valueindex(subM),stat=error)
    if(error/=0) stop
    allocate(quantabigbuf(4*subM,2),stat=error)
    if(error/=0) stop
    allocate(subspacenum((2*nsuborbs+1)*(2*nsuborbs+3)+1),stat=error)
    if(error/=0) stop
    subspacenum=0
    
    if(logic_spinreversal/=0) then
        allocate(szzeroindex(4*dim1),stat=error)
        if(error/=0) stop
        allocate(symmlinkbigbuf(4*dim1),stat=error)
        if(error/=0) stop
        symmlinkbigbuf=0
    end if
    
    ! n is the total number of basis
    n=0

    ! Sz loop
    do j=nsuborbs,-nsuborbs,-1
        
        if(logic_spinreversal/=0 .and. j<0) then
            ! the Sz>0 rotate matrix is the same as Sz<0 rotatematrix; only need copy
            exit
        end if
        
        ! partical number loop
        do i=0,2*nsuborbs,1
            
            ! m is the number of basis in each good quantum number subspace
            m=0
            
            ! every loop initiate the coeffwork
            do k=1,4*dim1,1
                call copy(coeffbuf(:,k),coeffwork(:,k))
            end do
            
            ! initiate the szzeroindex
            if(logic_spinreversal/=0) then
                szzeroindex=0
            end if

            do k=1,4*dim1,1
                if(quantabiginit(k,1)==i .and. quantabiginit(k,2)==j) then
                    m=m+1
                    ! swap the good quantum number subspace basis to
                    ! the first few rows/cols
                    call swap(coeffwork(:,m),coeffwork(:,k))
                    call swap(coeffwork(m,:),coeffwork(k,:))

                    if(logic_spinreversal/=0 .and. j==0) then
                        szzeroindex(m)=k   ! store the index of Sz=0 basis
                    end if
                end if
            end do

            if(m/=0) then
                ! Sz=0 condition
                ! when j==0 we first transform the basis to the new basis
                ! which the symmlink is him self
                if(logic_spinreversal/=0 .and. j==0 ) then
                    allocate(transform(m,m),stat=error)
                    if(error/=0) stop
                    transform=0.0D0
                
                    do l=1,m,1
                        if(abs(symmlinkbig(szzeroindex(l),1,Hindex))==szzeroindex(l)) then
                            transform(l,l)=1.0D0  ! the symmlink is himself
                        else
                            do q=l+1,m,1
                                if(symmlinkbig(szzeroindex(l),1,Hindex)==szzeroindex(q)) then
                                ! the fisrt column the parity is 1
                                ! the second column the parity is -1
                                    transform(l,l)=sqrt(2.0D0)/2.0D0
                                    transform(q,l)=sqrt(2.0D0)/2.0D0
                                    transform(l,q)=sqrt(2.0D0)/2.0D0
                                    transform(q,q)=-sqrt(2.0D0)/2.0D0
                                else if(symmlinkbig(szzeroindex(l),1,Hindex)==-szzeroindex(q)) then
                                    transform(l,l)=sqrt(2.0D0)/2.0D0
                                    transform(q,l)=-sqrt(2.0D0)/2.0D0
                                    transform(l,q)=sqrt(2.0D0)/2.0D0
                                    transform(q,q)=sqrt(2.0D0)/2.0D0
                                end if
                            end do
                        end if
                    end do
                    
                    allocate(coeffdummy(m,m),stat=error)
                    if(error/=0) stop
                    
                    ! transfer to the new basis
                    call gemm(coeffwork(1:m,1:m),transform,coeffdummy(1:m,1:m),'N','N',1.0D0,0.0D0)
                    call gemm(transform,coeffdummy(1:m,1:m),coeffwork(1:m,1:m),'T','N',1.0D0,0.0D0)
                    
                    ! swap the symmlink==1 to the first few columns 
                    ! and the symmlink==-1 to the last few columns
                    himp1=0  ! himself plus +1
                    do l=1,m,1
                        if(symmlinkbig(szzeroindex(l),1,Hindex)==szzeroindex(l)) then
                            himp1=himp1+1
                            call swap(coeffwork(:,himp1),coeffwork(:,l))
                            call swap(coeffwork(himp1,:),coeffwork(l,:))
                            call swap(transform(:,l),transform(:,himp1))  ! the transform matrix need to swap too
                        else if(symmlinkbig(szzeroindex(l),1,Hindex)==-szzeroindex(l)) then
                            cycle
                        else
                            done=.true.
                            ! check if have swaped
                            do l1=1,l-1,1
                                if(abs(symmlinkbig(szzeroindex(l),1,Hindex))==szzeroindex(l1)) then
                                    done=.false.
                                    exit
                                end if
                            end do
                            if(done==.true.) then
                                himp1=himp1+1
                                call swap(coeffwork(:,himp1),coeffwork(:,l))
                                call swap(coeffwork(himp1,:),coeffwork(l,:))
                                call swap(transform(:,l),transform(:,himp1))
                            end if
                        end if
                    end do
                    
                    ! diagonalize this matrix
                    call syevd(coeffwork(1:himp1,1:himp1),valuework(n+1:n+himp1),'V','U',info)
                    if(info/=0) then
                        write(*,*) "Sz=0 diagnolization failed! himp1"
                        stop
                    end if
                    symmlinkbigbuf(n+1:n+himp1)=1
                    ! symmlinkbigbuf==1 means the symmlink is himself
                    ! symmlinkbigbuf==-1 means the symmlink is -himself
                    call syevd(coeffwork(himp1+1:m,himp1+1:m),valuework(n+himp1+1:n+m),'V','U',info)
                    if(info/=0) then
                        write(*,*) "Sz=0 diagnolization failed! himm1"
                        stop
                    end if
                    symmlinkbigbuf(n+himp1+1:n+m)=-1
                    
                    ! transfer to its' inital basis representation
                    call gemm(transform,coeffwork(1:m,1:m),coeffdummy(1:m,1:m),'N','N',1.0D0,0.0D0)
                    coeffwork(1:m,1:m)=coeffdummy(1:m,1:m)
                    
                    deallocate(transform)
                    deallocate(coeffdummy)
                else  
                    ! Sz>0 condition
                    call syevd(coeffwork(1:m,1:m),valuework(n+1:n+m),'V','U',info)
                    if(info/=0) then
                        write(*,*) "Sz>0 diagnolization failed!"
                        stop
                    end if
                end if
                
                ! copy the U/V matrix to coeffresult/ it's inital location
                p=0
                do k=1,4*dim1,1
                    if(quantabiginit(k,1)==i .and. quantabiginit(k,2)==j) then
                        p=p+1
                        call copy(coeffwork(p,1:m),coeffresult(k,n+1:n+m))
                    end if
                end do

                if(p/=m) then
                    write(*,*) "p/=m failed!",p,m
                    stop
                end if
                
                ! update the quantabig
                quantabigbuf(n+1:n+m,1)=i
                quantabigbuf(n+1:n+m,2)=j
                
                ! subspacenum(1) stores the number of good quantum number subspace
                ! subspacenum(x) stores the number of basis in each subspace
                subspacenum(1)=subspacenum(1)+1
                subspacenum(subspacenum(1)+1)=m
            end if

        ! update the total basis
        n=m+n

        end do

        if(j>0) then
            szl0=n   ! store the Sz>0 basis
        else if(j==0 .and. logic_spinreversal/=0) then
            szzero=n-szl0   ! store the Sz=0 basis if needed
        end if
    end do
    

    if(logic_spinreversal/=0) then
        write(*,*) "szl0=",szl0,"szzero=",szzero
        if((szl0*2+szzero)/=4*dim1) then
            write(*,*) "(szl0*2+szzero)/=4*dim1 failed!","szl0",szl0,"szzero",szzero
            stop
        end if
        if(sum(subspacenum(2:subspacenum(1)+1))/=szl0+szzero) then
            write(*,*) "------------------------"
            write(*,*) "subspacenum=szl0+szzero,failed!",subspacenum
            write(*,*) "------------------------"
            stop
        end if
    else
        if(n/=4*dim1) then
            write(*,*) "------------------------"
            write(*,*) "n/=4*dim1,failed!",n
            write(*,*) "------------------------"
            stop
        end if
        if(sum(subspacenum(2:subspacenum(1)+1))/=4*dim1) then
            write(*,*) "------------------------"
            write(*,*) "subspacenum=4*dim1,failed!",subspacenum
            write(*,*) "------------------------"
            stop
        end if
    end if
    


! copy the Sz>0 part and Sz=0(symmetry pair is not himself) to the symmetry pair

    if(logic_spinreversal/=0) then
        do i=1,4*dim1,1
            if(quantabiginit(i,2)>0) then
                ! transfer every -1 link to 1 link
                coeffresult(abs(symmlinkbig(i,1,Hindex)),szl0+szzero+1:4*dim1)=coeffresult(i,1:szl0)*DBLE(sign(1,symmlinkbig(i,1,Hindex)))
            end if
        end do
        quantabigbuf(szl0+szzero+1:4*dim1,1)=quantabigbuf(1:szl0,1)
        quantabigbuf(szl0+szzero+1:4*dim1,2)=-1*quantabigbuf(1:szl0,2)
        valuework(szl0+szzero+1:4*dim1)=valuework(1:szl0)
    end if
    
    ! check if valuework<0
    do i=1,4*dim1,1
        if(valuework(i)<0.0D0) then
            if(abs(valuework(i))>1.0D-10) then
                write(*,*) "-----------------------------"
                write(*,*) "caution valuework<0.0D0",valuework(i)
                write(*,*) "-----------------------------"
            end if
            valuework(i)=0.0D0
        end if
    end do


    if(domain=='R' .and. (nstate==1 .or. exscheme==4)) then
        call RspaceCorrespond(valuework,subspacenum,quantabigbuf,szzero,szl0,valueindex)
    else 
        call selectstates(valuework,4*dim1,valueindex,singularvalue,subM,subspacenum,nsuborbs,szzero,szl0)
    end if

!------------------------------------------------------------------------
    ! check if the valueindex is right
    do i=1,subM,1
        if(valueindex(i)==0) then
            write(*,*) "----------------------------------"
            write(*,*) "splitsvd valueindex(i)==0",i
            write(*,*) "----------------------------------"
            stop
        end if
    end do
    do i=1,subM,1
        do j=i+1,subM,1
            if(valueindex(i)==valueindex(j)) then
                write(*,*) "----------------------------------"
                write(*,*) "splitsvd valueindex(i)=valueindex(j)",i,j,valueindex(i)
                write(*,*) "----------------------------------"
                stop
            end if
        end do
    end do
!------------------------------------------------------------------------
    
    ! copy the coeffresult to the U/V according to the valueindex
    do i=1,subM,1
        if(domain=='L') then
            call copy(coeffresult(:,valueindex(i)),leftu(:,i))
            quantasmaL(i,:)=quantabigbuf(valueindex(i),:)
        else if (domain=='R') then
            call copy(coeffresult(:,valueindex(i)),rightv(i,:))
            quantasmaR(i,:)=quantabigbuf(valueindex(i),:)
        end if

        if(logic_spinreversal/=0) then
            if(valueindex(i)<=szl0) then
                if((nstate==1 .or. exscheme==4) .and. domain=='R') then
                    symmlinksma(i,1,Hindex)=i-1
                else
                    symmlinksma(i,1,Hindex)=i+1   ! Sz>0
                end if
            else if(valueindex(i)>szl0 .and. valueindex(i)<=szl0+szzero) then
                if(symmlinkbigbuf(valueindex(i))==0) then
                    write(*,*) "----------------------------------------"
                    write(*,*) "symmlinkbigbuf(valueindex(i))==0 failed!"
                    write(*,*) "----------------------------------------"
                    stop
                end if
                symmlinksma(i,1,Hindex)=i*symmlinkbigbuf(valueindex(i))  ! Sz=0
            else 
                if((nstate==1 .or. exscheme==4) .and. domain=='R') then
                    symmlinksma(i,1,Hindex)=i+1
                else 
                    symmlinksma(i,1,Hindex)=i-1  ! Sz<0
                end if
            end if
        end if
    end do

    discard=1.0D0-sum(valuework(valueindex(1:subM)))
    write(*,'(A20,A1,D12.5)') "totaldiscard=",domain,discard

    deallocate(valuework)
    deallocate(coeffwork)
    deallocate(coeffbuf)
    deallocate(coeffresult)

    deallocate(valueindex)
    deallocate(quantabigbuf)
    deallocate(quantabiginit)
    deallocate(subspacenum)
    if(logic_spinreversal/=0) then
        deallocate(szzeroindex)
        deallocate(symmlinkbigbuf)
    end if

return
end subroutine splitsvd

!===================================================
!===================================================

subroutine RspaceCorrespond(valuework,subspacenum,quantabigbuf,&
        szzero,szl0,valueindex)
! this subroutine is to select R space basis
! which is corresponding the selected L space basis
! R space select states should be corresponse to the L space states
! when nstate==1 we can find the corresponding state in the L and R state
! with the same eigenvalue , but when nstate/=1 there is no such condition
    
    implicit none
    
    integer :: szzero,szl0
    integer :: subspacenum((2*(nright+1)+1)**2+1),quantabigbuf(4*subM,2)
    real(kind=r8) :: valuework(4*Rrealdim)
    integer :: valueindex(subM)   ! output
    
    ! local
    integer,parameter :: scalebound=15
    logical :: ifexist,iffind
    real(kind=r8) :: diffzero,scale1
    integer :: scalenum
    integer :: i,j,k,p
    
    
    valueindex=0
    
    do i=1,subM,1
        
        ! scale the eigenvalue to be exponential format
        scale1=1.0D0
        scalenum=0
        do while(.true.)
            if(singularvalue(i)>scale1) exit
            scale1=scale1*0.1D0
            scalenum=scalenum+1
            if(scalenum>=scalebound) then
                write(*,*) "---------------------"
                write(*,*) "caution! scalenum>15"
                write(*,*) "---------------------"
                exit
            end if
        end do
        
        if(logic_spinreversal==0) then
            if(scalenum<scalebound) then
                do j=1,subspacenum(1),1 
                    ! find the right good quantum number subspace
                    if((quantasmaL(i,1)+quantabigbuf(sum(subspacenum(2:j+1)),1)==nelecs) .and. &
                        (quantasmaL(i,2)+quantabigbuf(sum(subspacenum(2:j+1)),2)==totalSz)) then
                        diffzero=1.0D-10*(0.1D0**scalenum)
                        ! find the exact eigenvalue
                        do while(.true.)
                            do k=sum(subspacenum(2:j+1)),sum(subspacenum(2:j+1))-subspacenum(j+1)+1,-1
                                if(abs(valuework(k)-singularvalue(i))<diffzero) then
                                    ! check if have exist
                                    Ifexist=.false.
                                    do p=1,i-1,1
                                        if(valueindex(p)==k) then
                                            Ifexist=.true.
                                            exit
                                        end if
                                    end do
                                    if(Ifexist==.false.) then
                                        valueindex(i)=k
                                        exit
                                    end if
                                end if
                            end do
                            if(valueindex(i)/=0) exit
                            diffzero=diffzero*10.0D0
                        end do
                        exit
                    end if
                end do
            else
                ! find a random basis 
                ! singularvalue from big to small
                do j=1,4*Rrealdim,1
                    ifexist=.false.
                    do k=1,i-1,1
                        if(valueindex(k)==j) then
                            ifexist=.true.
                            exit
                        end if
                    end do
                    if(ifexist==.false.) then
                        valueindex(i)=j
                        exit
                    end if
                end do
            end if
        else 
            if(scalenum<scalebound) then
                ! only find the quantasmaL(i,2)<0 case
                if(quantasmaL(i,2)<=0) then
                    do j=1,subspacenum(1),1
                        if((quantasmaL(i,1)+quantabigbuf(sum(subspacenum(2:j+1)),1)==nelecs) .and. &
                            (quantasmaL(i,2)+quantabigbuf(sum(subspacenum(2:j+1)),2)==totalSz)) then
                            diffzero=1.0D-10*(0.1D0**scalenum)
                            do while(.true.)
                                do k=sum(subspacenum(2:j+1)),sum(subspacenum(2:j+1))-subspacenum(j+1)+1,-1
                                    if(abs(valuework(k)-singularvalue(i))<diffzero) then
                                        Ifexist=.false.
                                        do p=1,i-1,1
                                            if(valueindex(p)==k) then
                                                Ifexist=.true.
                                                exit
                                            end if
                                        end do
                                        if(Ifexist==.false.) then
                                            if(quantasmaL(i,2)<0) then
                                                valueindex(i)=k
                                                valueindex(i-1)=k+szl0+szzero
                                                exit
                                            else if(quantasmaL(i,2)==0) then
                                                valueindex(i)=k
                                                exit
                                            end if
                                        end if
                                    end if
                                end do
                                if(valueindex(i)/=0) exit
                                diffzero=diffzero*10.0D0
                            end do
                            exit
                        end if
                    end do
                end if

                if(quantasmaL(i,2)<=0 .and. valueindex(i)==0) then
                    write(*,*) "----------------------------"
                    write(*,*) "R space valueindex==0",i
                    write(*,*) "----------------------------"
                    stop
                end if
            else
                ! singularvalue small condition; random find
                if(i<subM .and. valueindex(i)==0) then
                    do j=szl0+1,2*szl0+szzero,1
                        ifexist=.false.
                        do k=1,i-1,1
                            if(valueindex(k)==j) then
                                ifexist=.true.
                                exit
                            end if
                        end do
                        if(ifexist==.false.) then
                            valueindex(i)=j
                            if(j>szl0+szzero) then
                                valueindex(i+1)=j-(szzero+szl0)
                            end if
                            exit
                        end if
                    end do
                else if(i==subM .and. valueindex(i)==0) then
                    ! find sz=0 basis
                    iffind=.false.
                    do  j=szl0+1,szl0+szzero,1
                        ifexist=.false.
                        do k=1,i-1,1
                            if(valueindex(k)==j) then
                                ifexist=.true.
                                exit
                            end if
                        end do
                        if(ifexist==.false.) then
                            valueindex(i)=j
                            iffind=.true.
                            exit
                        end if
                    end do
                    if(iffind==.false.) then
                        write(*,*) "-----------------------------------"
                        write(*,*) "R correspond did not find the last index valueindex=0"
                        write(*,*) "-----------------------------------"
                        stop
                    end if
                end if
            end if
        end if
    end do

return

end subroutine RspaceCorrespond

!===================================================
!===================================================

subroutine Prepare_LRcoeff(iLrealdim,iRrealdim,LRcoeff,istate)
    use blas95
    use noise_mod
    implicit none
    integer,intent(in) :: iLrealdim,iRrealdim,istate
    real(kind=r8),intent(out) :: LRcoeff(:,:)
    real(kind=r8) :: norm

    ! local
    integer :: job(8),info
    integer :: i

    if(myid==0) then
        ! recover coeffIF to its dense format
        LRcoeff=0.0D0
        if(ifopenperturbation==.false.) then
            job(1)=1
            job(2)=1
            job(3)=1
            job(4)=2
            job(5)=0
            job(6)=1
            call mkl_ddnscsr(job,4*iLrealdim,4*iRrealdim,LRcoeff(:,:),4*iLrealdim,coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindex(:,istate),info)
            if(info/=0) then
                call master_print_message(info,"coeffIF to dense format info/=")
                stop
            end if
        else
        ! recover the coeffIFp in coordinate format to dense format
            norm=dot(coeffIFp(1:coeffIFplast,1),coeffIFp(1:coeffIFplast,1))
            write(*,*) "Perturbation after normalization3:",1,norm
            do i=1,coeffIFplast,1
                LRcoeff(coeffIFrowindexp(i),coeffIFcolindexp(i))=coeffIFp(i,istate)
                !write(*,*) coeffIFrowindexp(i),coeffIFcolindexp(i)
                if(quantabigLp(coeffIFrowindexp(i),1)+quantabigRp(coeffIFcolindexp(i),1)/=nelecs &
                .or. quantabigLp(coeffIFrowindexp(i),2)+quantabigRp(coeffIFcolindexp(i),2)/=totalSz) then
                    write(*,*) "!!!!" 
                    write(*,*) i,coeffIFrowindexp(i),coeffIFcolindexp(i)
                end if
            end do
            write(*,*) "goodquantum number",coeffIfplast
        end if
    end if

    if(Ifnoise==.true.) then
        call svd_noise(iLrealdim,iRrealdim,LRcoeff)
    end if

    return
end subroutine Prepare_LRcoeff

!===================================================
!===================================================
end Module Renormalization_mod
