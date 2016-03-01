subroutine Sweep(direction)
! this subroutine is finit MPS procedure, sweep process

    USE mpi
    USE variables
    use communicate
    use exit_mod
    use Renormalization_mod
    use OnesiteMatrix
    use hamiltonian_mod
    use module_sparse
    use construct_system_big
    use checkmat_mod

    implicit none
    
    character(len=1) :: direction
    
    ! local 
    character(len=1) :: domain,envirodomain
    integer :: nsysorb, &      ! nsysorb=nleft or nright
               orbnow          ! orbnow=nleft+1 or norbs-nright
    integer :: i
    integer :: LRdim,LRdimp
    
    call master_print_message("enter in subroutine Sweep")

    if(direction=='l') then
        domain="L"
        envirodomain="R"
        nsysorb=nleft
        orbnow=nleft+1
        LRdim=Lrealdim
        LRdimp=Lrealdimp
    else if(direction=='r') then
        domain="R"
        envirodomain="L"
        nsysorb=nright
        orbnow=norbs-nright
        LRdim=Rrealdim
        LRdimp=Rrealdimp
    else
        call exit_DMRG(sigAbort,"Sweep subroutine direction/=l .and. /=r failed!")
    end if

    if(nsysorb/=exactsite) then
    !   if(myid==0) then
    !       write(*,*) Lrealdim,Lrealdimp
    !       write(*,*)  "Hsma-Hsmap"
    !       write(*,*) Hsmarowindex(Lrealdim+1,1)-1,Hsmarowindexp(Lrealdimp+1,1)-1
    !       write(*,*) Hsma(:,1)-Hsmap(:,1)
    !       write(*,*)  "Hbig-Hbigp"
    !       write(*,*) Hbigrowindex(4*Lrealdim+1,1)-1,Hbigrowindexp(4*Lrealdimp+1,1)-1
    !       write(*,*)  Hbig(:,1)-Hbigp(:,1)
    !   end if
    !   if(myid/=0) then
    !       write(*,*) Lrealdim,Lrealdimp
    !       write(*,*) "operabig"
    !       write(*,*) smarowindex1(Lrealdim+1,1:3)-1,smarowindex1p(Lrealdimp+1,1:3)-1
    !       write(*,*) operamatsma1(:,1:3)-operamatsma1p(:,1:3)
    !   end if
        call ConstructOnesiteMatrix(orbnow)
        
!       call checkmat(Hsma,Hsmacolindex,Hsmarowindex,&
!           operamatsma1,smacolindex1,smarowindex1,quantasmaL,quantasmaR,"sma",domain)
        
        call System_Big(domain,operamatsma1,smacolindex1,smarowindex1,&
            operamatbig1,bigcolindex1,bigrowindex1,&
            Hsma,Hsmacolindex,Hsmarowindex,&
            Hbig,Hbigcolindex,Hbigrowindex,&
            quantasmaL,quantasmaR,&
            LRdim,subM,.false.)
        if(ifOpenperturbation==.true.) then
            call System_Big(domain,operamatsma1p,smacolindex1p,smarowindex1p,&
                operamatbig1p,bigcolindex1p,bigrowindex1p,&
                Hsmap,Hsmacolindexp,Hsmarowindexp,&
                Hbigp,Hbigcolindexp,Hbigrowindexp,&
                quantasmaLp,quantasmaRp,&
                LRdimp,subMp,.true.)
        end if
    !   if(myid==0) then
    !       write(*,*) Lrealdim,Lrealdimp
    !       write(*,*)  "Hsma-Hsmap"
    !       write(*,*) Hsmarowindex(Lrealdim+1,1)-1,Hsmarowindexp(Lrealdimp+1,1)-1
    !       write(*,*) Hsma(:,1)-Hsmap(:,1)
    !       write(*,*)  "Hbig-Hbigp"
    !       write(*,*) Hbigrowindex(4*Lrealdim+1,1)-1,Hbigrowindexp(4*Lrealdimp+1,1)-1
    !       write(*,*)  Hbig(:,1)-Hbigp(:,1)
    !   end if
        
        if(domain=='L') then
            call System_Constructquanta(domain,Lrealdim,quantabigL(1:4*Lrealdim,1:2),quantasmaL(1:Lrealdim,1:2))
            if(ifopenperturbation==.true.) then
                call System_Constructquanta(domain,Lrealdimp,quantabigLp(1:4*Lrealdimp,1:2),quantasmaLp(1:Lrealdimp,1:2))
            end if
        else
            call System_Constructquanta(domain,Rrealdim,quantabigR(1:4*Rrealdim,1:2),quantasmaR(1:Rrealdim,1:2))
            if(ifopenperturbation==.true.) then
                call System_Constructquanta(domain,Rrealdimp,quantabigRp(1:4*Rrealdimp,1:2),quantasmaRp(1:Rrealdimp,1:2))
            end if
        end if

        if(logic_perturbation/=0 .and. ifopenperturbation==.false.) then
            call pre_perturbation(domain,'big')
        end if
        call Store_Operator(domain)
    else
        call Enviro_Big(domain)
    end if
    
    if(nleft==nright .and. logic_C2/=0) then
        call C2_Copy(direction)
        ! in the C2 symmetry, we calculate the two subspace
        ! first the reverse space; then the input space
        do i=1,2,1
            logic_C2=logic_C2*(-1)
            call master_print_message(logic_C2,"logic_C2=")
            call Hamiltonian(direction)
        end do
    else
        call Enviro_Big(envirodomain)
        call Hamiltonian(direction)
    end if

    if(isweep==sweeps .and. direction=='l' .and. nleft==(norbs-1)/2) then
        ! this is to do the chan proposed trace excited algrithom
        call master_print_message("In the last step, did not do Renormalization")
    else
        call Renormalization(direction)
    end if

return

end subroutine Sweep
