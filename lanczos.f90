subroutine LANCZOS_OP(localN,localP,localQ,localR)
    USE MPI
    use variables
    use kinds_mod
    use ABop
    integer :: localN
    real(kind=r8) :: localP(localN),localQ(localN),localR(localN)

    ! local 
    integer :: IJOB,ierr
    
    IJOB = 1
    call MPI_BCAST(IJOB,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
    
    !write(*,*) "localQ=",localQ
    !write(*,*) "localR=",localR

    if(opmethod=="comple") then
        call op(localN,1,localQ,localR,&
            Lrealdim,Rrealdim,subM,ngoodstates,&
            operamatbig1,bigcolindex1,bigrowindex1,&
            Hbig,Hbigcolindex,Hbigrowindex,&
            quantabigL,quantabigR,goodbasis,goodbasiscol)
    else if(opmethod=="direct") then
        call opdirect(localN,1,localQ,localR,&
            Lrealdim,Rrealdim,subM,ngoodstates,&
            operamatbig1,bigcolindex1,bigrowindex1,&
            Hbig,Hbigcolindex,Hbigrowindex,&
            quantabigL,quantabigR,goodbasis,goodbasiscol)
    else
        stop
    end if

    return

end subroutine LANCZOS_OP

!====================================================================
!====================================================================

subroutine LANCZOS_OPM(localN,localA,localB)
    
    use kinds_mod
    integer :: localN
    real(kind=r8) :: localA(localN),localB(localN)

    localB = localA
    
    return

end subroutine LANCZOS_OPM


