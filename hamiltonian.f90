Module Hamiltonian_mod
! this is the core Module
! contruct the hamiltonian on the fly and muliply it with coefficient of
! the wavefunction and transfter it to the 0 process
! use two different diagonalization method : Davidson and Jacobi-Davidson

    USE variables
    use communicate
    use mpi
    USE InitialGuess
    USE symmetry
    use module_sparse
    use masterdiag
    use GetHdiag_mod
    use perturbation_mod
    use ABop

    implicit none
    integer :: dimN
    real(kind=r8),allocatable :: HDIAG(:)

contains
!====================================================================
!====================================================================

subroutine Hamiltonian(direction)
    implicit none
    character(len=1) :: direction  ! i,l,r direction l is L space to R space sweep
    real(kind=r8) :: starttime,endtime
    
    starttime=MPI_WTIME()
    call master_print_message("enter hamiltonian subroutine")
    
    write(*,*) "myid=",myid
    write(*,*) "Lrealdim",Lrealdim
    write(*,*) "Rrealdim",Rrealdim
    write(*,*) "Lrealdimp",Lrealdimp
    write(*,*) "Rrealdimp",Rrealdimp

    if(diagmethod=="Davidson" .or. diagmethod=="D" .or. diagmethod=="MD") then
        call Davidson_wrapper(direction)
    else if(diagmethod=="JacobiDavidson" .or. diagmethod=="JD") then
        call JacobiDavidson_wrapper(direction)
    end if
    endtime=MPI_WTIME()
    call master_print_message(endtime-starttime,"DVDTIME:")
return

end subroutine Hamiltonian

!====================================================================
!====================================================================

subroutine JacobiDavidson_Wrapper(direction)

    implicit none
    character(len=1) :: direction  ! i,l,r direction l is L space to R space sweep
    ! local
    real(kind=r8),allocatable :: EIGS(:),RES(:),eigenvector(:)
!   integer,allocatable :: HDIAGrowindex(:),HDIAGcolindex(:)
    integer ::        NEIG         , &  ! number of eigenstate
                MADSPACE     , &  ! subspace dimension
                iter         , &
                NINIT        , &  ! initial Guess
                Leigenvector , &
                isearch      , &
                info         , &
                iprint       , &
                ijob         , &
                ICNTL(5)     , &
                NDX1         , &
                NDX2         
    real(kind=r8) ::  Tol          , &
                sigma        , &
                shift        , &
                mem          , &
                droptol      , &
                gap
    integer :: error,ierr
    integer :: i
    real(kind=r8) :: starttime,endtime
    

!--------------------------------------------------------------------
    ! initial value
    NEIG=nstate
    MADSPACE=20
    iter=2000
    mem=20.0
    droptol=1.0D-3
    ICNTL(1)=0
    ICNTL(2)=0
    ICNTL(3)=0
    ICNTL(4)=0
    ICNTL(5)=0
    IPRINT=6
    NINIT=nstate
    call GetDimSym
    Leigenvector=dimN*(3*MADSPACE+NEIG+1)+3*MADSPACE**2+MAX(MADSPACE**2,NEIG)+100

    if(myid==0) then
        allocate(EIGS(nstate))
        allocate(RES(nstate))
        allocate(eigenvector(Leigenvector))
    else
        allocate(eigenvector(1))
    end if
    
    if(myid==0) then
        if(isweep/=0 .and. nelecs==realnelecs) then
            isearch=1
            sigma=dmrgenergy(1)
            if(nstate/=1) then
                shift=dmrgenergy(1)*2.0D0-dmrgenergy(2)
            else 
                shift=dmrgenergy(1)-1.0D0
            end if
            EIGS(1:nstate)=dmrgenergy(1:nstate)
        else
            isearch=0
            EIGS(1:nstate)=0.0D0
        end if
    end if
        
    ! in the initial few sweeps
!   if(1.0D-3*(1.0D-1)**isweep>1.1D-6) then
!       Tol=1.0D-4*(1.0D-1)**isweep
!   else
!       Tol=1.0D-6
!   end if
!   ! in the last 3 sweeps
!   if(isweep==sweeps .or. isweep==sweeps-1 .or. &
!   isweep==sweeps-2) then
!       Tol=1.0D-4
!   end if
    if(isweep==0) then
        Tol=5.0D-3
    else
        Tol=1.0D-4
    end if
!--------------------------------------------------------------------

    if(myid==0) then
        allocate(HDIAG(dimN),stat=error)
        if(error/=0) stop
!       allocate(HDIAGrowindex(dimN+1))  ! only used in non-diagonal precondition method
!       allocate(HDIAGcolindex(dimN))
!       do i=1,dimN,1
!           HDIAGrowindex(i)=i
!           HDIAGcolindex(i)=i
!       end do
!       HDIAGrowindex(dimN+1)=dimN+1
    end if

    ! Get the diagonal element of hamiltonian
    starttime=MPI_WTIME()
    if(logic_spinreversal/=0 .or. &
        (logic_C2/=0 .and. nleft==nright)) then
        call SymmHDiag(HDIAG)
    else 
        call GetHDiag(HDIAG,dimN,&
        operamatbig1,bigcolindex1,bigrowindex1,&
        Hbig,Hbigcolindex,Hbigrowindex,&
        quantabigL,quantabigR,&
        .false.)
    end if
    endtime=MPI_WTIME()
    call master_print_message(endtime-starttime,"HDIAGTIME:")

    ! Get the Initialcoeff Guess
    if(myid==0) then
        call InitialStarter(direction,dimN,NINIT,eigenvector)
    end if

    if(myid/=0) then  ! these parameter nouse in slaver process
        dimN=1
        NDX1=1
        NDX2=1
    endif
    if(myid==0) then
        write(*,*) "isweep=",isweep,"site=",nleft+1,norbs-nright
    end if
    IJOB=0
    do while(.true.)
        if(myid==0) then
            call DPJDREVCOM(dimN,HDIAG,-1,-1,EIGS,RES,eigenvector,Leigenvector,NEIG,sigma,&
                isearch,NINIT,MADSPACE,ITER,TOL,SHIFT,DROPTOL,MEM,ICNTL,IJOB, &
                NDX1,NDX2,IPRINT,INFO,GAP)
        end if
        call MPI_BCAST(IJOB,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
        if(IJOB/=1) then
            exit
        else
            call op(dimN,1,eigenvector(NDX1),eigenvector(NDX2),&
                Lrealdim,Rrealdim,subM,ngoodstates,&
                operamatbig1,bigcolindex1,bigrowindex1,&
                Hbig,Hbigcolindex,Hbigrowindex,&
                quantabigL,quantabigR)
        end if
    end do

    if(myid==0) then
        if(INFO/=0) then
            call master_print_message(info,"INFO/=0")
            stop
        end if
        call DavidOutput(EIGS,eigenvector)
        write(*,*) "low state energy"
        do i=1,nstate,1
            write(*,*) nleft+1,norbs-nright,i,"th energy=",EIGS(i)
        end do
    end if
    
    if(ifopenperturbation==.true.) then
        call perturbation(EIGS,nstate)
    end if
    

    ! deallocate symmetry workarray
    if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
        call DestroySymm
    end if
    if(myid==0) then
        call DPJDCLEANUP
        deallocate(HDIAG)
    !   deallocate(HDIAGrowindex,HDIAGcolindex)
        deallocate(EIGS,RES)
    end if
    deallocate(eigenvector)

    return
end subroutine JacobiDavidson_Wrapper

!====================================================================
!====================================================================

Subroutine Davidson_Wrapper(direction)
! this is the davidson wrapper to call DVDSON written by Andreas
! aim to allocate memory and set variables in the global array
! mainly allocate memory on 0 process
    use MKL_SERVICE
    implicit none

    character(len=1) :: direction  ! i,l,r direction l is L space to R space sweep
    integer :: &
        lim   ,  &       ! the expanding small space in davidson iteration
        ilow  ,  &       ! index of lowest eigenpair 
        ihigh ,  &       ! index of the highest eigenpair
        niv   ,  &       ! number of initial vector
        mblock,  &       ! number of vector to be targeted in each iteration
        maxiter          ! maxiter iterations
    real(kind=8) :: &
        crite ,   &      ! convergence eigenvalue
        critc ,   &      ! convergence coefficiency
        critr ,   &      ! convergence residual vector norm
        ortho            ! orthogonality threshold
    integer,allocatable :: iselec(:)   ! selected eigenpair not used
    ! davidson parameter
    integer :: nloops,nmv,ierror,smadim,IWRSZ
    logical :: hiend
!   external op

    real(kind=r8),allocatable :: DavidWORK(:)
    real(kind=r8),allocatable :: dummycoeff(:),dummynewcoeff(:) ! have no use in fact
    integer :: error,i,j,k,m

    integer :: ierr ! MPI flag
    real(kind=r8) :: starttime,endtime
    
    starttime=MPI_WTIME()

    call master_print_message("enter in davidson diagonalization subroutine")
    
    ! default value
    ilow=1
    ihigh=nstate
    lim=nstate+35    ! this number 20 can be changed consider the efficiency
    niv=nstate
    mblock=nstate
    maxiter=400
    allocate(iselec(lim),stat=error)
    if(error/=0) stop
    iselec=-1

    ! in the initial few sweeps
    if(1.0D-5*(1.0D-1)**isweep>1.1D-10) then
        crite=1.0D-5*(1.0D-1)**isweep
        critc=1.0D-5*(1.0D-1)**isweep
        critr=1.0D-5*(1.0D-1)**isweep
        ortho=1.0D-4*(1.0D-1)**isweep
    else
        crite=1.0D-9
        critc=1.0D-9
        critr=1.0D-9
        ortho=1.0D-8
    end if

    ! in the last 3 sweeps
    if(isweep==sweeps .or. isweep==sweeps-1 .or. &
    isweep==sweeps-2) then
        crite=1.0D-10
        critc=1.0D-10
        critr=1.0D-10
        ortho=1.0D-9
    end if
    
    if(myid==0) then
        write(*,*) "criteria"
        write(*,*)  crite,critc,critr
    end if

    call GetDimSym
! allocate the davidson workarray needed by DVDSON
    IWRSZ=lim*(2*dimN+lim+9)+lim*(lim+1)/2+nstate+100
    if(myid==0) then
        allocate(HDIAG(dimN),stat=error)
        if(error/=0) stop
        allocate(DavidWORK(IWRSZ),stat=error)
        if(error/=0) stop
    end if
    
! Get the diagonal element of hamiltonian
    starttime=MPI_WTIME()
    if(logic_spinreversal/=0 .or. &
        (logic_C2/=0 .and. nleft==nright)) then
        call SymmHDiag(HDIAG)
    else 
        call GetHDiag(HDIAG,dimN,&
        operamatbig1,bigcolindex1,bigrowindex1,&
        Hbig,Hbigcolindex,Hbigrowindex,&
        quantabigL,quantabigR,&
        .false.)
    !   write(*,*) HDIAG
    end if
    endtime=MPI_WTIME()
    call master_print_message(endtime-starttime,"HDIAGTIME:")

! Get the Initialcoeff Guess
    if(myid==0) then
        call InitialStarter(direction,dimN,niv,DavidWORK)
    end if

!--------------------------------------------------------------------
! The core part of davidson diagnolization
    if(diagmethod=="D" .or. diagmethod=="Davidson") then
        if(myid==0) then
            call DVDSON(dimN,lim,HDIAG,ilow,ihigh,iselec &
                ,niv,mblock,crite,critc,critr,ortho,maxiter,DavidWORK,&
                IWRSZ,hiend,nloops,nmv,ierror)
                smadim=0
            call MPI_BCAST(smadim,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
        end if

        if(myid/=0) then
            do while(.true.)
                call MPI_BCAST(smadim,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
                if(smadim>0) then
                    allocate(dummycoeff(smadim),stat=error)
                    if(error/=0) stop
                    allocate(dummynewcoeff(smadim),stat=error)
                    if(error/=0) stop
                    call op(1,smadim,dummycoeff,dummynewcoeff,&
                            Lrealdim,Rrealdim,subM,ngoodstates,&
                            operamatbig1,bigcolindex1,bigrowindex1,&
                            Hbig,Hbigcolindex,Hbigrowindex,&
                            quantabigL,quantabigR)
                    deallocate(dummycoeff)
                    deallocate(dummynewcoeff)
                else
                    exit
                end if
            end do
        end if
    else if(diagmethod=="MD") then
        if(myid==0) then
            call MKL_SET_NUM_THREADS(nthreads(3))
        else
            call MKL_SET_NUM_THREADS(nthreads(4))
        end if

        call MasterGather
        if(myid==0) then
            call DVDSON(masterop,dimN,lim,HDIAG,ilow,ihigh,iselec &
                ,niv,mblock,crite,critc,critr,ortho,maxiter,DavidWORK,&
                IWRSZ,hiend,nloops,nmv,ierror)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        if(myid==0) then
            call MKL_SET_NUM_THREADS(nthreads(1))
        else
            call MKL_SET_NUM_THREADS(nthreads(2))
        end if
    end if

!-----------------------------------------------------------------------------
    
    if(myid==0) then
        if(hiend/=.false.) then
            call master_print_message("didn't get the lowest state")
            stop
        end if
        if(ierror/=0) then
            call master_print_message(ierror,"caution! IERROR/=0")
            if(ierror/=2048) then
                call master_print_message("failed!")
                stop
            end if
        end if
        write(*,*) "low state energy"
        do i=1,ihigh,1
            write(*,*) nleft+1,norbs-nright,i,"th energy=",DavidWORK(IHIGH*dimN+i)
            write(*,*) "energy converge:",DavidWORK(IHIGH*dimN+IHIGH+i)
            write(*,*) "residual norm:",DavidWORK(IHIGH*dimN+2*IHIGH+i)
        end do
        write(*,*) "NLOOPS=",nloops
        write(*,*) "NMV=",nmv
        call DavidOutput(DavidWORK(IHIGH*dimN+1),DavidWORK(1))
    end if

    if(ifopenperturbation==.true.) then
        if(myid==0) then
            call perturbation(DavidWORK(IHIGH*dimN+1:IHIGH*dimN+nstate),nstate)
        else
            call perturbation(DavidWORK,nstate)
        end if
    end if
    
    if(myid==0) then
        deallocate(HDIAG,DavidWORK)
    end if
    
    ! deallocate symmetry workarray
    if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
        call DestroySymm
    end if

    deallocate(iselec)
return

end subroutine Davidson_Wrapper

!====================================================================
!====================================================================

subroutine DavidOutput(eigenvalue,eigenvector)
    use perturbation_mod
    implicit none

    real(kind=r8) :: eigenvalue(nstate),eigenvector(nstate*dimN)
    real(kind=r8),allocatable :: nosymmout(:)
    integer :: i,k,error
    integer :: nonzero

    if(myid==0) then
        ! transfer the symmetry state to the non-symmetry state S*fai
        allocate(nosymmout(ngoodstates*nstate),stat=error)
        if(error/=0) stop

        if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
            do i=1,nstate,1
                call SymmetrizeState(ngoodstates,&
                    nosymmout((i-1)*ngoodstates+1:i*ngoodstates),eigenvector((i-1)*nsymmstate+1:i*nsymmstate),'u')
            end do
        else
            nosymmout=eigenvector(1:ngoodstates*nstate)
        end if
        !  no need do GramSchmit; the result fullfill orthnormal
        !  call GramSchmit(niv,ngoodstates,nosymmout,norm)
        !  write(*,*) "davidson nosymmout norm=",norm
         
        ! transfer the nosymmout to coeffIF
        do k=1,nstate,1
            call coefftosparse(&
                coeffIFdim,coeffIF(:,k),coeffIFcolindex(:,k),coeffIFrowindex(:,k),&
                ngoodstates,nosymmout((k-1)*ngoodstates+1:k*ngoodstates),&
                Lrealdim,Rrealdim,quantabigL(1:4*Lrealdim,1:2),quantabigR(1:4*Rrealdim,1:2))
            coeffIFrowindex(4*Lrealdim+1:4*subM+1,k)=coeffIFrowindex(4*Lrealdim+1,k)
        end do
    !   write the coeffIF in two partical density matrix calculation

        open(unit=109,file="coeffIF.tmp",form="unformatted",status="replace")
        do i=1,nstate,1
            write(109) coeffIFrowindex(:,i)
            nonzero=coeffIFrowindex(4*subM+1,i)-1
            write(109) coeffIF(1:nonzero,i)
            write(109) coeffIFcolindex(1:nonzero,i)
        end do
        close(109)
    end if

!   if(ifopenperturbation==.true.) then
!       call perturbation(eigenvalue(1:nstate),nstate)
!   end if
    
    if(ifopenperturbation==.false.) then
    if(myid==0) then
    ! update the sweepenergy
    ! use the middle site as the sweepenergy
        if(nleft==(norbs+1)/2-1) then
            do i=1,nstate,1
                sweepenergy(isweep,i)=eigenvalue(i)
            end do
        end if
        ! update the energy
        dmrgenergy(1:nstate)=eigenvalue(1:nstate)
    end if
    end if

    if(myid==0) then
        deallocate(nosymmout)
    end if

    return
end subroutine DavidOutput

!====================================================================
!====================================================================

subroutine GetDimSym
! Get dimension of H matrix and Symmetry matrix

    implicit none
    integer :: i,j
    ! check how many states fullfill good quantum number
    ! every process do it
    ! ngoodstates is the number of good quantum number states
    
    ! allocate the symmetry work array
    if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
        call SymmAllocateArray
    end if

    ngoodstates=0
    do i=1,4*Rrealdim,1
    do j=1,4*Lrealdim,1
        if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
        quantabigL(j,2)+quantabigR(i,2)==totalSz) then
            ngoodstates=ngoodstates+1
            ! construct the symmlinkgood 
            if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
                call CreatSymmlinkgood(ngoodstates,j,i)
            end if
        end if
    end do
        ! construct the symmlinkcol
        ! the nonzero LRcoeff element of every column
        if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
            symmlinkcol(i+1)=ngoodstates+1
        end if
    end do

    dimN=ngoodstates
    
    if((logic_spinreversal/=0 .or. &
        (logic_C2/=0 .and. nleft==nright))) then
        ! construct the symmetry matrix S in sparse format
        call SymmetryMat
        dimN=nsymmstate
    end if

    call master_print_message(ngoodstates,"ngoodstates:")
    call master_print_message(dimN,"total Hamiltonian dimension:")
return
end subroutine GetDimSym

!====================================================================
!====================================================================

end Module Hamiltonian_mod

