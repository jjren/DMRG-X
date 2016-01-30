subroutine pre_perturbation(domain,matchar)
! this subroutine used in the pre perturbation stage
! copy the 4m*4m to 4M*4M and big quanta array
! copy the m*m to M*M and small quanta array

    use variables
    use communicate
    use exit_mod
    use BLAS95
    use f95_precision
    use module_sparse
    implicit none
    
    character(len=1) :: domain
    character(len=3) :: matchar
    integer :: orbstart,orbend,Hindex,operaindex,nonzero
    integer :: dimsma,dimbig
    integer :: i,j
    
    if(matchar/="big" .and. matchar/="sma") then
        call master_print_message("matchar check wrong")
        stop
    end if

    if(domain=='L') then
        orbstart=1
        orbend=nleft+1
        Hindex=1
        dimbig=4*Lrealdim
        dimsma=Lrealdim
    else if(domain=='R') then
        orbstart=norbs-nright
        orbend=norbs
        Hindex=2
        dimbig=4*Rrealdim
        dimsma=Rrealdim
    else
        call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
    end if
    
    ! operamat
    do i=orbstart,orbend,1
        if(myid==orbid1(i,1)) then
            do j=1,3,1
                operaindex=orbid1(i,2)*3-3+j
                ! big matrix
                if(matchar=="big") then
                    bigrowindex1p(:,operaindex)=1  ! initialize
                    bigrowindex1p(1:dimbig+1,operaindex)=bigrowindex1(1:dimbig+1,operaindex)
                    nonzero=bigrowindex1(dimbig+1,operaindex)-1
                    bigcolindex1p(1:nonzero,operaindex)=bigcolindex1(1:nonzero,operaindex)
                    call copy(operamatbig1(1:nonzero,operaindex),operamatbig1p(1:nonzero,operaindex))
                else if(matchar=="sma") then
                    smarowindex1p(:,operaindex)=1
                    smarowindex1p(1:dimsma+1,operaindex)=smarowindex1(1:dimsma+1,operaindex)
                    nonzero=smarowindex1(dimsma+1,operaindex)-1
                    smacolindex1p(1:nonzero,operaindex)=smacolindex1(1:nonzero,operaindex)
                    call copy(operamatsma1(1:nonzero,operaindex),operamatsma1p(1:nonzero,operaindex))
                end if
            end do
        end if
    end do

    ! Hmat
    if(myid==0) then
        if(matchar=="big") then
            Hbigrowindexp(:,Hindex)=1
            Hbigrowindexp(1:dimbig+1,Hindex)=Hbigrowindex(1:dimbig+1,Hindex)
            nonzero=Hbigrowindex(dimbig+1,Hindex)-1
            Hbigcolindexp(1:nonzero,Hindex)=Hbigcolindex(1:nonzero,Hindex)
            call copy(Hbig(1:nonzero,Hindex),Hbigp(1:nonzero,Hindex))
        else if(matchar=="sma") then
            Hsmarowindexp(:,Hindex)=1
            Hsmarowindexp(1:dimsma+1,Hindex)=Hsmarowindex(1:dimsma+1,Hindex)
            nonzero=Hsmarowindex(dimsma+1,Hindex)-1
            Hsmacolindexp(1:nonzero,Hindex)=Hsmacolindex(1:nonzero,Hindex)
            call copy(Hsma(1:nonzero,Hindex),Hsmap(1:nonzero,Hindex))
        end if
    end if

    ! quanta array
    if(domain=='L') then
        if(matchar=="big") then
            quantabigLp=0
            quantabigLp(1:dimbig,1:2)=quantabigL(1:dimbig,1:2)
        else if(matchar=="sma") then
            quantasmaLp=0
            quantasmaLp(1:dimsma,1:2)=quantasmaL(1:dimsma,1:2)
        end if
    else if(domain=='R') then
        if(matchar=="big") then
            quantabigRp=0
            quantabigRp(1:dimbig,1:2)=quantabigR(1:dimbig,1:2)
        else if(matchar=="sma") then
            quantasmaRp=0
            quantasmaRp(1:dimsma,1:2)=quantasmaR(1:dimsma,1:2)
        end if
    end if
    
    ! Lrealdimp/Rrealdimp equals Lrealdim/Rrealdim in pre_perturbation case
    if(domain=='L') then
        Lrealdimp=Lrealdim
    else if(domain=="R") then
        Rrealdimp=Rrealdim
    end if

    return
end subroutine pre_perturbation






