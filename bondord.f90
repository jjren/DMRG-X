module BondOrder_mod
! this module calculate the bond order matrix;
! is the same as the one body density matrix a(i,sigma)^+ a(j,sigma)

    use module_sparse
    use variables
    use communicate
    use mathlib

    implicit none

    real(kind=r8),allocatable ::  transDM0(:,:,:,:,:),&
        transDMMO(:,:,:,:,:)
    
    contains
!===========================================================================
!===========================================================================

subroutine init_BOmat(orbindex)
! initiate the on site niup,nidown matrix
    use onesitematrix

    implicit none

    integer :: orbindex
    ! local
    integer :: operaindex2
    
    if(myid==orbid2(orbindex,orbindex,1)) then
        call ConstructOnesiteMatrix(orbindex)
        operaindex2=orbid2(orbindex,orbindex,2)
        smarowindex2(:,operaindex2*2-1:operaindex2*2)=0
        call CopySpAtoB(4,onesitemat(:,11),osmcolindex(:,11),osmrowindex(:,11)&
            ,operamatsma2(:,2*operaindex2-1),smacolindex2(:,2*operaindex2-1),&
            smarowindex2(:,2*operaindex2-1),smadim2)
        call CopySpAtoB(4,onesitemat(:,8),osmcolindex(:,8),osmrowindex(:,8)&
            ,operamatsma2(:,2*operaindex2),smacolindex2(:,2*operaindex2),&
            smarowindex2(:,2*operaindex2),smadim2)
    end if
return

end subroutine init_BOmat

!===========================================================================
!===========================================================================

subroutine BondOrder
! this subroutine calculate the BondOrder Matrix in the last step
    implicit none

    integer :: error
    integer :: i,j,k,istate
    real(kind=8) :: tmp

    call master_print_message("enter BondOrder subroutine")
    
    if(myid==0) then
        allocate(transDM0(norbs,norbs,2,nstate,nstate),stat=error)
        if(error/=0) stop
        transDM0=0.0D0
    end if
    
    call Calc_BOmat_link
    call Calc_BOmat_subspace('L')
    call Calc_BOmat_subspace('R')
    
    if(myid==0) then

        !  transition density matrix
        open(unit=398,file="AO-transOpdm.out",form="unformatted",status="replace")
        write(398) norbs,nstate
        do istate=1,nstate,1
        do k=istate,nstate,1  ! including the one partical density matrix
            do i=1,norbs,1
            do j=1,norbs,1
                if(i==j) then    ! recover niup+nidown,niup-nidown to niup,nidown
                    tmp=transDM0(i,i,1,k,istate)-transDM0(i,i,2,k,istate)
                    transDM0(i,i,1,k,istate)=(transDM0(i,i,1,k,istate)+transDM0(i,i,2,k,istate))/2.0D0
                    transDM0(i,i,2,k,istate)=tmp/2.0D0
                end if
            end do
            end do
            write(398) transDM0(:,:,:,k,istate)
        end do
        end do
        close(398)

        ! bondorder/one partical density matrix
        call master_print_message("bondorder matrix")
        open(unit=399,file="bondord.out",status="replace")
        do k=1,nstate,1
            write(399,*) k
            do i=1,norbs,1
            do j=i,norbs,1
                if(bondlink(i,j)==1) then
                    write(399,*) i,j,transDM0(i,j,:,k,k)
                end if
            end do
            end do
        end do
        
        ! electron density
        do k=1,nstate,1
            write(399,*) k
            do i=1,norbs,1
                write(399,*) i,i,transDM0(i,i,:,k,k)
            end do
        end do
        close(399)

        if(logic_bondorder==2 .and. mod(nelecs,2)==0 .and. logic_meanfield==1) then  ! need the Opdm as basis
            ! transform the one partical reduced matrix to MO basis
            allocate(transDMMO(norbs,norbs,2,nstate,nstate))
            call transDMAO2MO
            ! transform it to the natural orbital basis
            call NatOrbAnalysis
            if(nstate>1) then
                ! NTO analysis
                call NatTraOrb
            end if
            call proxyNAC
            deallocate(transDMMO)
        end if
    end if

    if(myid==0) then
        deallocate(transDM0)
    end if
return
end subroutine BondOrder

!===========================================================================
!===========================================================================

subroutine proxyNAC
! calculate the approximate Nonadiabatic coupling used in 
! Anna I. Krylov JPCL 2013 4 3845

    use blas95
    use f95_precision

    implicit none
    real(kind=r8) :: pNAC(nstate,nstate)
    integer :: i,j,k,l

    pNAC=0.0D0
    call master_print_message("proxNAC")
    do i=1,nstate,1
        write(*,*) "istate=",i
        do j=i+1,nstate,1
            do l=1,2,1
            do k=1,norbs,1
                pNAC(j,i)=pNAC(j,i)+nrm2(transDM0(:,k,l,j,i))**2
            end do
            end do
            write(*,*) pNAC(j,i)
        end do
    end do
return

end subroutine proxyNAC

!===========================================================================
!===========================================================================

subroutine transDMAO2MO
! the transition density matrix or one partical density matrix from PPP-AO to MO
! MOij=CDC^+   C is MO*AO column format
    use BLAS95
    use F95_PRECISION
    implicit none
    integer :: i,j,k,l,istate2
    real(kind=r8),allocatable :: midmat(:,:),coeffC(:,:)
    real(kind=r8) tmp
    integer :: itmp

    allocate(midmat(norbs,norbs))
    allocate(coeffC(norbs,norbs))

    ! transition density matrix
    call master_print_message("MO transition density matrix")
    ! read MO coefficient
    open(unit=157,file="MO.out",status="old")
    do i=1,norbs,1
        read(157,*) itmp,tmp
        read(157,*) coeffC(:,i) 
    end do
    close(157)

    open(unit=151,file="MO-transOpdm.out",form="unformatted",status="replace")
    write(151) norbs,nstate
    do i=1,nstate,1
    do istate2=i,nstate,1
        do j=1,2,1   ! spin up down
            call gemm(coeffC,transDM0(:,:,j,istate2,i),midmat,'T','N',1.0D0,0.0D0)
            call gemm(midmat,coeffC,transDMMO(:,:,j,istate2,i),'N','N',1.0D0,0.0D0)
        end do
        write(151) transDMMO(:,:,:,istate2,i)
        
        if(i==1 .and. istate2/=i) then
            write(*,*) "istate=",istate2
            do k=1,norbs,1
            do l=1,norbs,1
                tmp=abs(transDMMO(k,l,1,istate2,1)+transDMMO(k,l,2,istate2,1))
                if(tmp>0.44) then
                    write(*,*) k,"<<--",l,sqrt(2.0D0)/2.0D0*tmp
                end if
            end do
            end do
        end if
    end do
    end do
    close(151)

    deallocate(midmat,coeffC)
return
end subroutine transDMAO2MO

!===========================================================================
!===========================================================================

subroutine NatOrbAnalysis
! this subroutine aims to transfer the one partical reduced density matrix
! to MO basis and Natural orbital basis
! analyse the N odd using the Anna I.Krylov proposed method 
! JPCL 2013,4,3845
    use meanfield
    use LAPACK95
    use BLAS95
    USE F95_PRECISION
    implicit none

    real(kind=r8),allocatable :: MOonepdm(:,:,:),AOonepdm(:,:,:),midmat(:,:),NOcoeff(:,:),&
    NOeigenvalue(:,:)
    integer :: i,j
    ! all are the spatial one partical RDM
    allocate(AOonepdm(norbs,norbs,nstate))
    allocate(NOcoeff(norbs,norbs))
    allocate(NOeigenvalue(norbs,nstate))
    allocate(midmat(norbs,norbs))

    call master_print_message("MO diagonal Opdm")
    do i=1,nstate,1
        write(*,*) "istate=",i
        AOonepdm(:,:,i)=transDM0(:,:,1,i,i)+transDM0(:,:,2,i,i)
        ! MO occupation
        do j=1,norbs,1
            write(*,*) transDMMO(j,j,1,i,i)+transDMMO(j,j,2,i,i)
        end do
    end do

    call master_print_message("NO diagonal Opdm")
    open(unit=153,file="NO.out",status="replace")
    do i=1,nstate,1
        call  syevr(AOonepdm(:,:,i),NOeigenvalue(:,i),'U',NOcoeff)
        write(*,*) "istate=",i
        do j=1,norbs,1
            write(*,*) NOeigenvalue(j,i)
        end do
        ! output NO orbital
        write(153,*) "istate=",i
        do j=1,norbs,1
            write(153,*) j
            write(153,*) NOcoeff(:,j)
        end do
    end do
    close(153)

    call master_print_message("odd electron numbers")
    do i=1,nstate,1
        write(*,*) 2*realnelecs-(sum(NOeigenvalue(:,i)**2))
    end do

    deallocate(AOonepdm,NOcoeff,NOeigenvalue,midmat)

return
end subroutine NatOrbAnalysis

!===========================================================================
!===========================================================================

subroutine NatTraOrb
! this subroutine do natural transition orbital(NTO) analysis if the ex is
! single excitation
! <ex|a^+a*ai|gs>
    
    use LAPACK95
    USE F95_PRECISION
    use meanfield
    implicit none
    real(kind=r8),allocatable :: Tai(:,:,:),svdvalue(:),&
    rightv(:,:),leftu(:,:),ww(:)
    integer :: info,tmp
    integer :: i

    allocate(Tai(norbs-nocc,nocc,nstate))
    
    ! the operator is the singlet excitation operator 
    ! sqrt(2)/2*(apup^+*aqup+apdown^+*aqdown)
    Tai(:,:,:)=sqrt(2.0D0)/2.0D0*(TransDMMO(nocc+1:norbs,1:nocc,1,:,1)+TransDMMO(nocc+1:norbs,1:nocc,2,:,1))
    
    tmp=min(nocc,norbs-nocc)
    allocate(svdvalue(tmp))
    allocate(leftu(norbs-nocc,tmp))
    allocate(rightv(tmp,nocc))
    allocate(ww(tmp-1))
    
    call master_print_message("NTO analysis result")
    open(unit=150,file="NTO.out",status="replace")
    do i=2,nstate,1
        call gesvd(Tai(:,:,i),svdvalue,leftu,rightv,ww,'N',info)
        if(info/=0) then
            write(*,*) "NatTraOrb info/=0",info
            stop
        end if
        write(*,*) "stateindex",i
        write(*,*) "svdvalue",svdvalue
        write(150,*) i
        write(150,*) svdvalue
        write(150,*) leftu
        write(150,*) rightv
    end do
    close(150)
    
    deallocate(Tai,leftu,rightv,ww,svdvalue)

return
end subroutine NatTraOrb

!===========================================================================
!===========================================================================

subroutine Calc_BOmat_subspace(domain)
    
    use OpExpec_mod
    implicit none

    character(len=1),intent(in) :: domain
    ! local
    real(kind=r8) :: midratio
    integer :: orbstart,orbend,iproc,operaindex
    integer :: lstate,rstate,iorb,jorb,ispin
    real(kind=r8),allocatable :: midtrans(:,:)

    midratio=pppVmidratio

    if(domain=="L") then
        orbstart=1
        orbend=nleft+1
    else
        orbstart=norbs-nright
        orbend=norbs
    end if

    do iorb=orbstart,orbend,1
    do jorb=iorb,orbend,1
    do ispin=1,2,1
        if(iorb==jorb .and. ispin==1) then
            ! L space <n_{iup}+n_{idown}-nuclQ_i>
            iproc=orbid1(iorb,1)
            if(myid==0 .or. myid==iproc) then
                do lstate=1,nstate,1
                do rstate=1,nstate,1
                    operaindex=orbid1(iorb,2)*3
                    call SubSpaceOpExpec(transDM0(iorb,iorb,ispin,lstate,rstate),iproc,domain,&
                            Lrealdim,Rrealdim,subM,lstate,rstate,operaindex,&
                            operamatbig1,bigcolindex1,bigrowindex1,&
                            coeffIF,coeffIfcolindex,coeffIFrowindex,midratio)
                end do
                end do
            end if
        else
            
            iproc=orbid2(iorb,jorb,1)
            if(myid==0 .or. myid==iproc) then
                do lstate=1,nstate,1
                do rstate=1,nstate,1
                    operaindex=orbid2(iorb,jorb,2)*2-2+ispin
                    call SubSpaceOpExpec(transDM0(iorb,jorb,ispin,lstate,rstate),iproc,domain,&
                            Lrealdim,Rrealdim,subM,lstate,rstate,operaindex,&
                            operamatbig2,bigcolindex2,bigrowindex2,&
                            coeffIF,coeffIfcolindex,coeffIFrowindex,midratio)
                end do
                end do
            end if
        end if

        if(myid==0) then
            do lstate=1,nstate,1
            do rstate=1,nstate,1
                if(iorb==jorb .and. ispin==1 .and. lstate==rstate) then
                    transDM0(iorb,iorb,ispin,lstate,rstate)=transDM0(iorb,iorb,ispin,lstate,rstate)+nuclQ(iorb)
                end if
            end do
            end do
            
            
            if(iorb/=jorb) then
                allocate(midtrans(nstate,nstate))
                midtrans=transDM0(iorb,jorb,ispin,:,:)
                
                do lstate=1,nstate,1
                do rstate=1,nstate,1
                    if(domain=="L") then
                        transDM0(jorb,iorb,ispin,lstate,rstate)=midtrans(rstate,lstate)
                    else
                        transDM0(iorb,jorb,ispin,lstate,rstate)=midtrans(rstate,lstate)
                        transDM0(jorb,iorb,ispin,lstate,rstate)=midtrans(lstate,rstate)
                    end if
                end do
                end do

                deallocate(midtrans)
            end if
        end if
    end do
    end do
    end do

    return

end subroutine Calc_BOmat_subspace

!===========================================================================
!===========================================================================

subroutine Calc_BOmat_link
    use OpExpec_mod
    implicit none
    
    real(kind=r8),allocatable :: phase(:)
    integer :: lstate,rstate,ileft,iright,ispin
    integer :: lproc,rproc,loperaindex,roperaindex,leadproc,i
    real(kind=r8) :: midratio

    midratio=LRoutratio

    allocate(phase(4*Lrealdim))
    do i=1,4*Lrealdim,1
        phase(i)=(-1.0D0)**(mod(quantabigL(i,1),2))
    end do

    do ileft=1,nleft+1,1
    do iright=norbs-nright,norbs,1
    if(bondlink(ileft,iright)/=0 .or. logic_bondorder==2) then
        do ispin=1,2,1
            lproc=orbid1(ileft,1)
            rproc=orbid1(iright,1)
            loperaindex=orbid1(ileft,2)*3-3+ispin
            roperaindex=orbid1(iright,2)*3-3+ispin
            leadproc=lproc
            
            if(myid==lproc .or. myid==rproc .or. myid==0) then
                do lstate=1,nstate,1
                do rstate=1,nstate,1
                    call LinkOpExpec(transDM0(ileft,iright,ispin,lstate,rstate),'N','T',lproc,rproc,leadproc,&
                            Lrealdim,Rrealdim,subM,lstate,rstate,&
                            loperaindex,roperaindex,operamatbig1,bigcolindex1,bigrowindex1,&
                            coeffIF,coeffIfcolindex,coeffIFrowindex,bigratio1,midratio,.true.,phase)
                end do
                end do
                
                if(myid==0) then
                    do lstate=1,nstate,1
                    do rstate=1,nstate,1
                        transDM0(iright,ileft,ispin,lstate,rstate) = transDM0(ileft,iright,ispin,rstate,lstate)
                    end do
                    end do
                end if
            end if
        end do
    end if
    end do
    end do

    deallocate(phase)

    return
end subroutine Calc_BOmat_link

!===========================================================================
!===========================================================================

subroutine LocalMagMoment(domain)
! local spin
    use OpExpec_mod
    implicit none

    character(len=1),intent(in) :: domain
    
    ! local
    real(kind=r8) :: midratio
    integer :: orbstart,orbend,iproc,operaindex
    integer :: istate,iorb
    real(kind=r8),allocatable ::  localmagmoment0(:,:)
    character(len=50) :: filename

    midratio=pppVmidratio

    if(domain=="L") then
        orbstart=1
        orbend=nleft+1
    else
        orbstart=norbs-nright
        orbend=norbs
    end if

    if(myid==0) allocate(localmagmoment0(norbs,nstate))

    do iorb=orbstart,orbend,1
        iproc=orbid2(iorb,iorb,1)

        if(myid==0 .or. myid==iproc) then
            do istate=1,nstate,1
                operaindex=orbid2(iorb,iorb,2)*2-1
                call SubSpaceOpExpec(localmagmoment0(iorb,istate),iproc,domain,&
                        Lrealdim,Rrealdim,subM,istate,istate,operaindex,&
                        operamatbig2,bigcolindex2,bigrowindex2,&
                        coeffIF,coeffIfcolindex,coeffIFrowindex,midratio)
            end do
        end if
    end do

    if(myid==0) then
        filename=trim(domain)//"localmagmoment.out"
        open(unit=666,file=filename,status="replace")
        do istate=1,nstate,1
            write(666,*)  istate
            do iorb=orbstart,orbend,1
                write(666,*) iorb,localmagmoment0(iorb,istate)
            end do
        end do
        deallocate(localmagmoment0)
    end if
            
    return

end subroutine LocalMagMoment

!===========================================================================
!===========================================================================
end module BondOrder_mod
