module BondOrder_mod
! this module calculate the bond order matrix;
! is the same as the one body density matrix a(i,sigma)^+ a(j,sigma)

    use module_sparse
    use variables
    use communicate
    
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

        operamatsma2(1:4,2*operaindex2-1:2*operaindex2)=onesitemat(1:4,7:8)
        smacolindex2(1:4,2*operaindex2-1:2*operaindex2)=osmcolindex(1:4,7:8)
        smarowindex2(1:5,2*operaindex2-1:2*operaindex2)=osmrowindex(1:5,7:8)
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
            ! electron density
            do i=1,norbs,1
                write(399,*) i,i,transDM0(i,i,:,k,k)
            end do
        end do
        close(399)

        if(logic_bondorder==2) then  ! need the Opdm as basis
            ! transform the one partical reduced matrix to MO basis
            call transDMAO2MO
            ! transform it to the natural orbital basis
            call NatOrbAnalysis
            if(nstate>1) then
                ! NTO analysis
                call NatTraOrb
            end if
            call proxyNAC
        end if
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
    use MeanField
    use BLAS95
    use F95_PRECISION
    implicit none
    integer :: i,j,k,l,istate2
    real(kind=r8),allocatable :: midmat(:,:)
    real(kind=r8) tmp

    allocate(transDMMO(norbs,norbs,2,nstate,nstate))
    allocate(midmat(norbs,norbs))

    ! transition density matrix
    call master_print_message("MO transition density matrix")
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

    deallocate(midmat)
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
    
    deallocate(Tai,leftu,rightv,ww,svdvalue)

return
end subroutine NatTraOrb

!===========================================================================
!===========================================================================

subroutine Calc_BOmat_subspace(domain)
! calculate operator in the L/R subspace i,j<=nleft+1,or i,j>=norbs-nright
! <R|<L|CLR ai^+aj CL'R'|L'>|R'>
    use exit_mod
    use mpi
    use mathlib
    implicit none
    
    character(len=1) :: domain
    
    ! local
    integer :: i,j,istate,k,l,istate2
    integer :: operaindex2,orbstart,orbend
    real(kind=r8) :: itransDM(2,nstate,nstate)
    ! itransDM :: the first 2 means up/down, the second and third is the transition pair
    integer :: error,ierr
    integer :: nmid
    real(kind=r8),allocatable :: midmat(:)
    integer(kind=i4),allocatable :: midcolindex(:),midrowindex(:),coeffIFrowindexdummy(:,:)
    
    integer :: status(MPI_STATUS_SIZE) ! MPI flag

    if(myid/=0) then
        nmid=CEILING(DBLE(16*subM*subM)/pppVmidratio)
        allocate(midmat(nmid),stat=error)
        if(error/=0) stop
        allocate(midcolindex(nmid),stat=error)
        if(error/=0) stop
        allocate(midrowindex(4*subM+1),stat=error)
        if(error/=0) stop
        allocate(coeffIFrowindexdummy(4*subM+1,nstate),stat=error)
        if(error/=0) stop
    end if

    ! set the parameters
    if(domain=='L') then
        orbstart=1
        orbend=nleft+1
    else if(domain=='R') then
        orbstart=norbs-nright
        orbend=norbs
    else
        call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
    end if
    
    ! two operator matrix => no phase
    
    ! for example :: CLR*CL'R'OLL'IRR'=CLR*(OLL'CL'R') the same as transition moment
    ! in the R domain need to tranpose the coeffIF
    if(myid/=0) then
        if(domain=='R') then
            do istate=1,nstate,1
                call CSCtoCSR('RC',4*Rrealdim,4*Lrealdim,&
                coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindex(:,istate),&
                coeffIFrowindexdummy(:,istate))
            end do
        else
            coeffIFrowindexdummy=coeffIFrowindex
        end if
    end if

    do i=orbstart,orbend,1
    do j=i,orbend,1
        if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
            if(myid==orbid2(i,j,1)) then
                do istate=1,nstate,1
                do k=1,2,1
                    operaindex2=orbid2(i,j,2)*2-2+k
                    call SpMMtoSp('N','N',4*subM,4*subM,4*subM,4*subM,4*subM,&
                        operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2), &
                        coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
                        midmat,midcolindex,midrowindex,nmid)
                    !=====================================================================================
                    ! calculate the transition density matrix
                    ! including bondorder matrix
                    ! <psai1|ai^+*aj|psai2>/=<psai1|aj^+*ai|psai2>
                    ! <ex|ai^+*aj|gs>=<gs|aj^+*ai|ex>
                    do istate2=1,nstate,1
                        ! trace(CLR*QLR) or trace(CRL*QRL)
                        call SpMMtrace('T',4*subM,&
                            coeffIF(:,istate2),coeffIFcolindex(:,istate2),coeffIFrowindexdummy(:,istate2),&
                            midmat,midcolindex,midrowindex,itransDM(k,istate2,istate))
                    end do
                    !=====================================================================================
                end do
                end do
                call MPI_SEND(itransDM,2*nstate**2,mpi_real8,0,orbid2(i,j,2),MPI_COMM_WORLD,ierr)
            
            else if(myid==0) then
            ! transition density matrix 
                call MPI_RECV(itransDM,2*nstate**2,mpi_real8,orbid2(i,j,1),orbid2(i,j,2),MPI_COMM_WORLD,status,ierr)
                do istate=1,nstate,1
                do istate2=istate,nstate,1
                do k=1,2,1
                    if(domain=='L') then   ! in the L space l=1 means (i,j) pair
                        transDM0(i,j,k,istate2,istate)=itransDM(k,istate2,istate)  ! i<j
                        transDM0(j,i,k,istate2,istate)=itransDM(k,istate,istate2)
                    else if(domain=='R') then   ! in the R space l=1 means (j,i) pair
                        ! <psai1|ai^+*aj|psai2>=<psai2|aj^+*ai|psai1>
                        transDM0(i,j,k,istate2,istate)=itransDM(k,istate,istate2)  ! i<j  in calculation we store in the R space (i>j,ai^+aj)
                        transDM0(j,i,k,istate2,istate)=itransDM(k,istate2,istate)
                    end if
                end do
                end do
                end do
            end if
        end if
    end do
    end do

    ! recovery the coeffIF
    if(myid/=0) then
        if(domain=='R') then
            do istate=1,nstate,1
                call CSCtoCSR('CR',4*Lrealdim,4*Rrealdim,&
                coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
                coeffIFrowindex(:,istate))
            end do
        end if
    end if

    if(myid/=0) deallocate(midmat,midcolindex,midrowindex,coeffIFrowindexdummy)

return

end subroutine Calc_BOmat_subspace

!===========================================================================
!===========================================================================

subroutine Calc_BOmat_link
! this subroutine calculate bond order matrix belong to different subspace
! i<=nleft+1,j>=norbs-nright 
! the same as op subrouitne transfer integral algorithm
    use mathlib
    use mpi
    use blas95
    use F95_PRECISION
    
    implicit none
    include "mkl_spblas.fi"

    integer :: hopnelement,midnelement
    real(kind=r8),allocatable :: hopmat(:,:),midmat(:),midmat2(:)
    integer(kind=i4),allocatable :: &
    hopmatcol(:,:),hopmatrow(:,:),&
    midmatcol(:),midmatrow(:),&
    midmatcol2(:),midmatrow2(:),&
    phase(:)

    character(len=1),allocatable :: hoppackbuf(:)
    integer :: position1,hoppacksize
    integer :: operaindex,nnonzero
    real(kind=r8) :: itransDM(2,nstate,nstate)
    integer :: i,j,k,l,m,istate,p,istate2
    logical :: ifhop,ifhopsend
    integer :: hoptouched(nprocs-1),hopntouched
    integer :: error,info
    
    integer :: status(MPI_STATUS_SIZE),hopsendrequest(nprocs-1)
    integer :: ierr

    hopnelement=CEILING(DBLE(16*subM*subM)/bigratio1)
    midnelement=CEILING(DBLE(16*subM*subM)/hopmidratio)
    
    do i=1,nleft+1,1
    do j=norbs,norbs-nright,-1
        if(bondlink(i,j)==1 .or. logic_bondorder==2) then
            if(myid==orbid1(i,1) .or. myid==orbid1(j,1)) then
                if(.not. allocated(hopmat)) then
                    allocate(hopmat(hopnelement,2),stat=error) ! store the hopping matrix
                    if(error/=0) stop
                    allocate(hopmatcol(hopnelement,2),stat=error) 
                    if(error/=0) stop
                    allocate(hopmatrow(4*subM+1,2),stat=error) 
                    if(error/=0) stop

                    hoppacksize=(hopnelement*12+4*(4*subM+1))*2+1000  ! 1000 is redundant
                    allocate(hoppackbuf(hoppacksize),stat=error) ! packbuf to send hopping matrix
                    if(error/=0) stop
                end if
            end if

            if(myid==orbid1(i,1)) then
                if(.not. allocated(midmat)) then
                    allocate(midmat(midnelement),stat=error) ! store the intermediate matrix
                    if(error/=0) stop
                    allocate(midmatcol(midnelement),stat=error) 
                    if(error/=0) stop
                    allocate(midmatrow(4*subM+1),stat=error) 
                    if(error/=0) stop
                    allocate(midmat2(midnelement),stat=error) ! store the intermediate matrix
                    if(error/=0) stop
                    allocate(midmatcol2(midnelement),stat=error) 
                    if(error/=0) stop
                    allocate(midmatrow2(4*subM+1),stat=error) 
                    if(error/=0) stop
                end if
            end if
        end if
    end do
    end do
    
    do i=norbs,norbs-nright,-1
        if(myid==orbid1(i,1)) then
            ! check if need hopping matrix
            ifhop=.false.
            do j=1,nleft+1,1
                if(bondlink(i,j)==1 .or. logic_bondorder==2) then
                    ifhop=.true.
                    exit
                end if
            end do
            
            if(ifhop==.true.) then
                operaindex=orbid1(i,2)
                
                ! copy operamatbig to hopmat
                ! send operamatbig, not R*C
                do l=1,2,1
                    ! integer can not call copy
                    hopmatrow(:,l)=bigrowindex1(:,operaindex*3-3+l)
                    nnonzero=bigrowindex1(4*subM+1,operaindex*3-3+l)-1
                    hopmatcol(1:nnonzero,l)=bigcolindex1(1:nnonzero,operaindex*3-3+l)
                    call copy(operamatbig1(1:nnonzero,operaindex*3-3+l),hopmat(1:nnonzero,l))
                end do

                ! pack hopmatrix
                position1=0
                do l=1,2,1
                    call MPI_PACK(hopmatrow(1,l),(4*subM+1),MPI_integer4,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
                    call MPI_PACK(hopmat(1,l),(hopmatrow(4*subM+1,l)-1),MPI_real8,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
                    call MPI_PACK(hopmatcol(1,l),(hopmatrow(4*subM+1,l)-1),MPI_integer4,hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
                end do

                ! send the hopmat
                hoptouched=0
                hopntouched=0
                do l=1,nleft+1,1
                    if((bondlink(i,l)==1 .or. logic_bondorder==2) .and. orbid1(l,1)/=myid) then ! if orbid(l)==myid need not send hopmat
                        ifhopsend=.false.
                        do m=1,hopntouched,1
                            if(orbid1(l,1)==hoptouched(m)) then
                                ifhopsend=.true.   ! have send hopmat to this process
                                exit
                            end if
                        end do
                        if(ifhopsend==.false.) then
                            hopntouched=hopntouched+1
                            hoptouched(hopntouched)=orbid1(l,1)
                            call MPI_ISEND(hoppackbuf,position1,MPI_PACKED,orbid1(l,1),i,MPI_COMM_WORLD,hopsendrequest(hopntouched),ierr)
                        end if
                    end if
                end do
            end if
        end if

        ! recv hopmat
        ! hopping term calculation
        if(myid/=orbid1(i,1) .and. myid/=0) then
            do l=1,nleft+1,1
                if((bondlink(l,i)==1 .or. logic_bondorder==2) .and. myid==orbid1(l,1)) then  
                    call MPI_RECV(hoppackbuf,hoppacksize,MPI_PACKED,orbid1(i,1),i,MPI_COMM_WORLD,status,ierr)
                    position1=0
                    do j=1,2,1
                        call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmatrow(1,j),(4*subM+1),MPI_integer4,MPI_COMM_WORLD,ierr)
                        call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmat(1,j),(hopmatrow(4*subM+1,j)-1),MPI_real8,MPI_COMM_WORLD,ierr)
                        call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopmatcol(1,j),(hopmatrow(4*subM+1,j)-1),MPI_integer4,MPI_COMM_WORLD,ierr)
                    end do
                    exit  ! only recv once
                end if
            end do
        end if
        
        do l=1,nleft+1,1
            if(bondlink(i,l)==1 .or. logic_bondorder==2) then
            if(myid==orbid1(l,1)) then
                operaindex=orbid1(l,2)
                
                !-----------------------------------------------------
                ! the +1 -1 phase added to l' 
                allocate(phase(4*subM),stat=error)
                if(error/=0) stop

                do j=1,4*subM,1
                    phase(j)=(-1.0D0)**(mod(quantabigL(j,1),2))
                end do
                !----------------------------------------------------

                ! construct hopmat
                do j=1,nstate,1
                    do k=1,2,1
                        ! firstly needed to transfer from a(+) to a, so 'T' and 'T' is 'N'
                        ! k=1 a up,k=2 a down
                        call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
                                coeffIF(:,j),coeffIFcolindex(:,j),coeffIFrowindex(:,j), &
                                hopmat(:,k),hopmatcol(:,k),hopmatrow(:,k), &
                                midmat,midmatcol,midmatrow,midnelement,info)
                        call checkinfo(info)

                        ! add phase
                        do p=1,4*subM,1
                        do m=midmatrow(p),midmatrow(p+1)-1,1
                            midmat(m)=midmat(m)*phase(p)
                        end do
                        end do

                        !k<=2 al^+*ar
                        call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
                                operamatbig1(:,operaindex*3-3+k),bigcolindex1(:,operaindex*3-3+k),bigrowindex1(:,operaindex*3-3+k), &
                                midmat,midmatcol,midmatrow,&
                                midmat2,midmatcol2,midmatrow2,midnelement,info)
                        call checkinfo(info)

                        !=========================================================================
                        ! calculate transition density matrix
                        ! includeing the bondorder matrix
                        ! <ex|aR^+*aL|gs>=<gs|aL^+*aR|ex>
                        do istate2=1,nstate,1
                            ! trace(CLR*OLR)
                            call SpMMtrace('T',4*subM, & 
                                    coeffIF(:,istate2),coeffIFcolindex(:,istate2),coeffIFrowindex(:,istate2), &
                                    midmat2,midmatcol2,midmatrow2,itransDM(k,istate2,j))
                        end do
                        !=========================================================================
                    end do
                end do
                call MPI_SEND(itransDM,2*nstate**2,mpi_real8,0,1,MPI_COMM_WORLD,ierr)
                deallocate(phase)

            else if(myid==0) then
                call MPI_RECV(itransDM,2*nstate**2,mpi_real8,orbid1(l,1),1,MPI_COMM_WORLD,status,ierr)
                do istate=1,nstate,1
                do istate2=istate,nstate,1
                    transDM0(i,l,1:2,istate2,istate)=itransDM(1:2,istate,istate2)
                    transDM0(l,i,1:2,istate2,istate)=itransDM(1:2,istate2,istate)
                end do
                end do
            end if
            end if
        end do

        ! confirm that the pppVmat and hopmat can be used again without problem
        if(myid==orbid1(i,1) .and. ifhop==.true.) then
            do j=1,hopntouched,1
                call MPI_WAIT(hopsendrequest(j),status,ierr)
            end do
        end if
    end do

    if(allocated(hopmat)) deallocate(hopmat,hopmatcol,hopmatrow)
    if(allocated(midmat)) deallocate(midmat,midmatcol,midmatrow)
    if(allocated(midmat2)) deallocate(midmat2,midmatcol2,midmatrow2)
    if(allocated(hoppackbuf)) deallocate(hoppackbuf)

return
end subroutine Calc_BOmat_link

!===========================================================================
!===========================================================================
end module BondOrder_mod
