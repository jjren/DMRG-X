module module_sparse
! this module contains the sparse format matrix in 3 array CSR format
! operamatbig and Hbig
! the core workarray of this program

    use kinds_mod
    use variables
    
    implicit none
    private
    save

    public :: AllocateArray

    ! sparse form in 3 array CSR format
    real(kind=r8),allocatable,public :: &
    operamatbig1(:,:)  , &        ! sparse form 1 electron operamatbig
    operamatsma1(:,:)  , &        ! sparse form 1 electron operamatsma
    operamatbig2(:,:)  , &        ! sparse form 2 electron operamatbig
    operamatsma2(:,:)  , &        ! sparse form 2 electron operamatsma
    operamatbig3(:,:)  , &        ! sparse form local spin operamatbig
    operamatsma3(:,:)  , &        ! sparse form local spin operamatsma
    Hbig(:,:)          , &        ! Hbig in sparse form
    Hsma(:,:)          , &        ! Hsma in sparse form
    coeffIF(:,:)       , &        ! coeffIF is the inital and final wavefunction coefficient 
    operamatbig1p(:,:) , &        ! sparse form 1 electron operamatbig in perturbation mode
    operamatsma1p(:,:) , &        ! sparse form 1 electron operamatsma in perturbation mode
    Hbigp(:,:)         , &        ! Hbig in sparse form in perturbation mode
    Hsmap(:,:)         , &        ! Hsma in sparse form in perturbation mode
    coeffIFp(:,:)                 ! coeffIF is the inital and final wavefunction coefficient in perturbation mode
    
    integer(kind=i4),allocatable,public :: &
    bigrowindex1(:,:) , &         ! 1 electron operamatbig rowindex
    bigcolindex1(:,:) , &         ! 1 electron oepramatbig columnindex
    smarowindex1(:,:) , &         ! 1 electron operamatsma rowindex
    smacolindex1(:,:) , &         ! 1 electron operamatsma columnindex
    bigrowindex2(:,:) , &         ! 2 electron operamatbig rowindex
    bigcolindex2(:,:) , &         ! 2 electron oepramatbig columnindex
    smarowindex2(:,:) , &         ! 2 electron operamatsma rowindex
    smacolindex2(:,:) , &         ! 2 electron operamatsma columnindex
    bigrowindex3(:,:) , &         ! local spin operamatbig rowindex
    bigcolindex3(:,:) , &         ! local spin oepramatbig columnindex
    smarowindex3(:,:) , &         ! local spin operamatsma rowindex
    smacolindex3(:,:) , &         ! local spin operamatsma columnindex
    Hbigcolindex(:,:) , &         ! Hbig colindex
    Hbigrowindex(:,:) , &         ! Hbig rowindex
    Hsmacolindex(:,:) , &         ! Hsma colindex
    Hsmarowindex(:,:) , &         ! Hsma rowindex
    coeffIFcolindex(:,:) ,&       ! coeffIF colindex
    coeffIFrowindex(:,:) ,&       ! coeffIF rowindex
    bigrowindex1p(:,:) , &         ! 1 electron operamatbig rowindex in perturbation mode
    bigcolindex1p(:,:) , &         ! 1 electron oepramatbig columnindex in perturbation mode
    smarowindex1p(:,:) , &         ! 1 electron operamatsma rowindex in perturbation mode
    smacolindex1p(:,:) , &         ! 1 electron operamatsma columnindex in perturbation mode
    Hbigcolindexp(:,:) , &         ! Hbig colindex in perturbation mode
    Hbigrowindexp(:,:) , &         ! Hbig rowindex in perturbation mode
    Hsmacolindexp(:,:) , &         ! Hsma colindex in perturbation mode
    Hsmarowindexp(:,:) , &         ! Hsma rowindex in perturbation mode
    coeffIFcolindexp(:) ,&         ! coeffIF colindex in perturbation mode
    coeffIFrowindexp(:)            ! coeffIF rowindex in perturbation mode
    
    integer(kind=i4),public :: coeffIFplast   ! the number of nonzero element in the coeffIFp
    
    ! in sparse form operamatbig/operamatsma,Hbig/Hsma dim
    integer(kind=i4),public :: &
    bigdim1,smadim1   , & 
    bigdim2,smadim2   , &
    bigdim3,smadim3   , & 
    Hbigdim,Hsmadim   , &
    coeffIFdim        , &
    bigdim1p,smadim1p , &
    Hbigdimp,Hsmadimp , &
    coeffIFdimp

    ! sparse parameter
    real(kind=r8),public :: pppmatratio,hopmatratio,LRoutratio,UVmatratio,coeffIFratio
    real(kind=r8),public :: bigratio1,smaratio1,bigratio2,smaratio2,bigratio3,smaratio3,Hbigratio,Hsmaratio  ! sparse radio
    
    integer,allocatable,public :: operanum1(:),operanum2(:),operanum3(:)
    ! store the number of operators on every process
    ! operanum1 is the max site operator every process have
    
    real(kind=r8),allocatable,public:: moperamatbig1(:,:)
    integer(kind=i4),allocatable,public :: mbigrowindex1(:,:),mbigcolindex1(:,:)
    contains

!=========================================================================================================
!=========================================================================================================

subroutine AllocateArray
    
    use communicate
    implicit none
    
    
    ! local
    integer :: error
    integer :: C2state
    
    call sparse_default

! set the sparse mat dim
    bigdim1=CEILING(DBLE(16*subM*subM)/bigratio1)
    bigdim2=CEILING(DBLE(16*subM*subM)/bigratio2)
    bigdim3=CEILING(DBLE(16*subM*subM)/bigratio3)
    smadim1=CEILING(DBLE(subM*subM)/smaratio1)
    smadim2=CEILING(DBLE(subM*subM)/smaratio2)
    smadim3=CEILING(DBLE(subM*subM)/smaratio3)
    Hbigdim=CEILING(DBLE(16*subM*subM)/Hbigratio)
    Hsmadim=CEILING(DBLE(subM*subM)/Hsmaratio)
    coeffIFdim=CEILING(DBLE(16*subM*subM)/coeffIFratio)
    
    ! perturbation mode
    bigdim1p=CEILING(DBLE(16*subMp*subMp)/bigratio1)
    smadim1p=CEILING(DBLE(subMp*subMp)/smaratio1)
    Hbigdimp=CEILING(DBLE(16*subMp*subMp)/Hbigratio)
    Hsmadimp=CEILING(DBLE(subMp*subMp)/Hsmaratio)
    coeffIFdimp=CEILING(DBLE(16*subMp*subMp)/coeffIFratio)

! allocate memory 
    if(myid/=0) then
        allocate(operamatbig1(bigdim1,3*operanum1(myid)),stat=error)
        if(error/=0) stop
        allocate(bigcolindex1(bigdim1,3*operanum1(myid)),stat=error)
        if(error/=0) stop
        allocate(bigrowindex1(4*subM+1,3*operanum1(myid)),stat=error)
        if(error/=0) stop

        allocate(operamatsma1(smadim1,3*operanum1(myid)),stat=error)
        if(error/=0) stop
        allocate(smacolindex1(smadim1,3*operanum1(myid)),stat=error)
        if(error/=0) stop
        allocate(smarowindex1(subM+1,3*operanum1(myid)),stat=error)
        if(error/=0) stop
        bigrowindex1=1   ! set the matrix to be 0
        smarowindex1=1

        if(logic_perturbation/=0) then
            allocate(operamatbig1p(bigdim1p,3*operanum1(myid)))
            allocate(bigcolindex1p(bigdim1p,3*operanum1(myid)))
            allocate(bigrowindex1p(4*subMp+1,3*operanum1(myid)))

            allocate(operamatsma1p(smadim1p,3*operanum1(myid)))
            allocate(smacolindex1p(smadim1p,3*operanum1(myid)))
            allocate(smarowindex1p(subMp+1,3*operanum1(myid)))
            bigrowindex1p=1   ! set the matrix to be 0
            smarowindex1p=1
        end if

        if(logic_bondorder/=0) then
            allocate(operamatbig2(bigdim2,2*operanum2(myid)),stat=error)
            if(error/=0) stop
            allocate(bigcolindex2(bigdim2,2*operanum2(myid)),stat=error)
            if(error/=0) stop
            allocate(bigrowindex2(4*subM+1,2*operanum2(myid)),stat=error)
            if(error/=0) stop

            allocate(operamatsma2(smadim2,2*operanum2(myid)),stat=error)
            if(error/=0) stop
            allocate(smacolindex2(smadim2,2*operanum2(myid)),stat=error)
            if(error/=0) stop
            allocate(smarowindex2(subM+1,2*operanum2(myid)),stat=error)
            if(error/=0) stop
            bigrowindex2=1
            smarowindex2=1
        end if
        
        ! local spin operator matrix
        if(logic_localspin==1) then
            allocate(operamatbig3(bigdim3,operanum3(myid)),stat=error)
            if(error/=0) stop
            allocate(bigcolindex3(bigdim3,operanum3(myid)),stat=error)
            if(error/=0) stop
            allocate(bigrowindex3(4*subM+1,operanum3(myid)),stat=error)
            if(error/=0) stop

            allocate(operamatsma3(smadim3,operanum3(myid)),stat=error)
            if(error/=0) stop
            allocate(smacolindex3(smadim3,operanum3(myid)),stat=error)
            if(error/=0) stop
            allocate(smarowindex3(subM+1,operanum3(myid)),stat=error)
            if(error/=0) stop
            bigrowindex3=1
            smarowindex3=1
        end if
    else
    ! 2 means the R space ;1 means the L space
    
        allocate(Hbig(Hbigdim,2),stat=error)
        if(error/=0) stop
        allocate(Hbigcolindex(Hbigdim,2),stat=error)
        if(error/=0) stop
        allocate(Hbigrowindex(4*subM+1,2),stat=error)
        if(error/=0) stop
        
        allocate(Hsma(Hsmadim,2),stat=error)
        if(error/=0) stop
        allocate(Hsmacolindex(Hsmadim,2),stat=error)
        if(error/=0) stop
        allocate(Hsmarowindex(subM+1,2),stat=error)
        if(error/=0) stop
        
        ! in C2 mode we can calculate the two subspace together
        ! without loss of accuracy
        if(logic_C2==0) then
            C2state=nstate
        else
            C2state=2*nstate
        end if

        allocate(coeffIF(coeffIFdim,C2state),stat=error)
        if(error/=0) stop
        allocate(coeffIFcolindex(coeffIFdim,C2state),stat=error)
        if(error/=0) stop
        allocate(coeffIFrowindex(4*subM+1,C2state),stat=error)
        if(error/=0) stop
        
        Hbigrowindex=1
        Hsmarowindex=1
        coeffIFrowindex=1
        
        if(logic_perturbation/=0) then
            allocate(Hbigp(Hbigdimp,2))
            allocate(Hbigcolindexp(Hbigdimp,2))
            allocate(Hbigrowindexp(4*subMp+1,2))
            
            allocate(Hsmap(Hsmadimp,2))
            allocate(Hsmacolindexp(Hsmadimp,2))
            allocate(Hsmarowindexp(subMp+1,2))
            
            allocate(coeffIFp(coeffIFdimp,C2state))
            allocate(coeffIFcolindexp(coeffIFdimp))
            allocate(coeffIFrowindexp(coeffIFdimp))
            Hbigrowindexp=1
            Hsmarowindexp=1
            coeffIFrowindexp=1
        end if

        if(diagmethod=="MD") then
            allocate(moperamatbig1(bigdim1,3*norbs))
            allocate(mbigcolindex1(bigdim1,3*norbs))
            allocate(mbigrowindex1(4*subM+1,3*norbs))
        end if
    end if

return

end subroutine AllocateArray

!=========================================================================================================
!=========================================================================================================

subroutine sparse_default
! set the default ratio according to the subM
    use communicate
    implicit none

    if(subM<=128) then
        bigratio1=1.0
        smaratio1=1.0
        bigratio2=1.0
        smaratio2=1.0
        bigratio3=1.0
        smaratio3=1.0
        Hbigratio=1.0
        Hsmaratio=1.0
        pppmatratio=1.0
        hopmatratio=1.0
        LRoutratio=1.0
        UVmatratio=1.0
        coeffIFratio=1.0
    else if(abs(subM-120)<20) then
        bigratio1=32.0
        smaratio1=8.0
        bigratio2=32.0
        smaratio2=8.0
        bigratio3=1.0
        smaratio3=1.0
        Hbigratio=12.0
        Hsmaratio=8.0
        pppmatratio=10.0
        hopmatratio=12.0
        LRoutratio=10.0
        UVmatratio=5.0
        coeffIFratio=10.0
    else if (abs(subM-256)<50) then
        bigratio1=20.0
        smaratio1=5.0
        bigratio2=20.0
        smaratio2=5.0
        bigratio3=1.0
        smaratio3=1.0
        Hbigratio=10.0
        Hsmaratio=5.0
        pppmatratio=10.0
        hopmatratio=12.0
        LRoutratio=10.0
        UVmatratio=5.0
        coeffIFratio=10.0
    else if (abs(subM-512)<50) then
        bigratio1=15.0
        smaratio1=4.0
        bigratio2=15.0
        smaratio2=4.0
        bigratio3=15.0
        smaratio3=4.0
        Hbigratio=10.0
        Hsmaratio=5.0
        pppmatratio=8.0
        hopmatratio=10.0
        LRoutratio=5.0
        UVmatratio=5.0
        coeffIFratio=8.0
    else if (abs(subM-1024)<100) then
        bigratio1=25.0
        smaratio1=5.0
        bigratio2=25.0
        smaratio2=5.0
        bigratio3=15.0
        smaratio3=4.0
        Hbigratio=10.0
        Hsmaratio=5.0
        pppmatratio=10.0
        hopmatratio=12.0
        LRoutratio=10.0
        UVmatratio=5.0
        coeffIFratio=10.0
    else if (abs(subM-2048)<500) then
        bigratio1=27.0
        smaratio1=5.0
        bigratio2=27.0
        smaratio2=5.0
        bigratio3=27.0
        smaratio3=5.0
        Hbigratio=10.0
        Hsmaratio=5.0
        pppmatratio=10.0
        hopmatratio=12.0
        LRoutratio=10.0
        UVmatratio=5.0
        coeffIFratio=10.0
    else
        bigratio1=15.0
        smaratio1=4.0
        bigratio2=15.0
        smaratio2=4.0
        bigratio3=15.0
        smaratio3=4.0
        Hbigratio=10.0
        Hsmaratio=5.0
        pppmatratio=8.0
        hopmatratio=10.0
        LRoutratio=5.0
        UVmatratio=5.0
        coeffIFratio=8.0
    end if

    if(myid==0) then
        write(*,*) "bigratio1=",    bigratio1
        write(*,*) "smaratio1=",    smaratio1
        write(*,*) "bigratio2=",    bigratio2
        write(*,*) "smaratio2=",    smaratio2
        write(*,*) "bigratio3=",    bigratio3
        write(*,*) "smaratio3=",    smaratio3
        write(*,*) "Hbigratio=",    Hbigratio
        write(*,*) "Hsmaratio=",    Hsmaratio
        write(*,*) "pppmatratio=",  pppmatratio
        write(*,*) "hopmatratio=",  hopmatratio
        write(*,*) "LRoutratio=" ,  LRoutratio
        write(*,*) "UVmatratio=" ,  UVmatratio
        write(*,*) "coeffIFratio=", coeffIFratio
    end if
return

end subroutine sparse_default

!=========================================================================================================
!=========================================================================================================

end module module_sparse
