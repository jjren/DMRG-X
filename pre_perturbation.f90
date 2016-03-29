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
    use mathlib
    implicit none
    
    character(len=1) :: domain
    character(len=3) :: matchar
    integer :: orbstart,orbend,Hindex,operaindex
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
                    call CopySpAtoB(dimbig,operamatbig1(:,operaindex),bigcolindex1(:,operaindex),&
                            bigrowindex1(:,operaindex),operamatbig1p(:,operaindex),&
                            bigcolindex1p(:,operaindex),bigrowindex1p(:,operaindex),bigdim1p) 
                else if(matchar=="sma") then
                    smarowindex1p(:,operaindex)=1
                    call CopySpAtoB(dimsma,operamatsma1(:,operaindex),smacolindex1(:,operaindex),&
                            smarowindex1(:,operaindex),operamatsma1p(:,operaindex),&
                            smacolindex1p(:,operaindex),smarowindex1p(:,operaindex),smadim1p) 
                end if
            end do
        end if
    end do

    ! Hmat
    if(myid==0) then
        if(matchar=="big") then
            Hbigrowindexp(:,Hindex)=1
            call CopySpAtoB(dimbig,Hbig(:,Hindex),Hbigcolindex(:,Hindex),Hbigrowindex(:,Hindex),&
                    Hbigp(:,Hindex),Hbigcolindexp(:,Hindex),Hbigrowindexp(:,Hindex),Hbigdimp)
        else if(matchar=="sma") then
            Hsmarowindexp(:,Hindex)=1
            call CopySpAtoB(dimsma,Hsma(:,Hindex),Hsmacolindex(:,Hindex),Hsmarowindex(:,Hindex),&
                    Hsmap(:,Hindex),Hsmacolindexp(:,Hindex),Hsmarowindexp(:,Hindex),Hsmadimp)
        end if
    end if

    ! quanta array
    if(domain=='L') then
        if(matchar=="big") then
            quantabigLp=0
            do i=1,2,1
                call scopy(dimbig,quantabigL(:,i),1,quantabigLp(:,i),1)
            end do
        else if(matchar=="sma") then
            quantasmaLp=0
            do i=1,2,1
                call scopy(dimsma,quantasmaL(:,i),1,quantasmaLp(:,i),1)
            end do
        end if
    else if(domain=='R') then
        if(matchar=="big") then
            quantabigRp=0
            do i=1,2,1
                call scopy(dimbig,quantabigR(:,i),1,quantabigRp(:,i),1)
            end do
        else if(matchar=="sma") then
            quantasmaRp=0
            do i=1,2,1
                call scopy(dimsma,quantasmaR(:,i),1,quantasmaRp(:,i),1)
            end do
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






