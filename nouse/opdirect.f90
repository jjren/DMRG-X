module ABop

contains

subroutine op(bigdim,smadim,coeff,newcoeff,&
iLrealdim,iRrealdim,isubM,nosymmdim,&
cap_big,cap_bigcol,cap_bigrow,&
cap_Hbig,cap_Hbigcol,cap_Hbigrow,&
cap_quantabigL,cap_quantabigR,cap_goodbasis,cap_goodbasiscol)
! this is the core subroutine to calculate the S*H*S*C or H*C
! the parallel schema follow JCP 12 3174(2004) garnet chan
! if want to save memory, then can write a wrapper, to send one coeff every time

!--------------------------------------------------------
! input bigdim,smadim,coeff
! bigdim is the totaldim 16M*M
! bigdim may be < 16M*M because we use good quantum number and spin symmetry,
! and the dimension  may be half
! if groud state smadim=1
! if gs+ex smadim may be >1
!--------------------------------------------------------
! output newcoeff
! coeff is the input coefficient and in the 1-d arrary format
! new coeff is H cross C result
!---------------------------------------------------------

    use mpi
    use variables
    use symmetry
    use mathlib
    use module_sparse
    use checkmem_mod

    implicit none
    include "mkl_spblas.fi"

    integer,intent(in) :: bigdim,smadim
    real(kind=r8),intent(in) :: coeff(bigdim*smadim)
    real(kind=r8),intent(out) :: newcoeff(bigdim*smadim)
    integer,intent(in) :: iLrealdim,iRrealdim,isubM,nosymmdim
    real(kind=r8),intent(in) :: cap_big(:,:),cap_Hbig(:,:)
    integer(kind=i4),intent(in) :: cap_bigcol(:,:),cap_bigrow(:,:),&
                         cap_Hbigcol(:,:),cap_Hbigrow(:,:),&
                         cap_quantabigL(:,:),cap_quantabigR(:,:),&
                         cap_goodbasis(:,:),cap_goodbasiscol(:)
    ! local
    integer :: operaindex
    
    real(kind=r8),allocatable :: LRcoeffin(:,:),LRcoeffout(:,:),coeffnosymm(:),coeffnosymmreduce(:)
    integer(kind=i4),allocatable :: LRcoeffincol(:,:),LRcoeffinrow(:,:),&
    LRcoeffoutcol(:,:),LRcoeffoutrow(:,:)
    
    real(kind=r8),allocatable :: hopbufmat(:,:),hopmat(:,:),&
        pppVmat(:),pppVbufmat(:),buffmat(:),midmat(:)
    integer(kind=i4),allocatable :: &
        hopmatcol(:,:),hopmatrow(:,:),&
        hopbufcol(:,:),hopbufrow(:,:),&
        pppVmatcol(:),pppVmatrow(:),&
        pppVbufcol(:),pppVbufrow(:),&
        buffmatcol(:),buffmatrow(:),&
        midmatcol(:),midmatrow(:)

    integer(kind=i4) :: pppVnelement,hopnelement,LRoutnelement,&
        nrecv,npppVbuf,npppVmidmat,nhopbuf,nhopmidmat

    character(len=1),allocatable :: hoppackbuf(:),pppVpackbuf(:)
    integer :: position1,pppVpacksize,hoppacksize

    real(kind=r8),allocatable :: phase(:)
    integer :: m,ie,il,ir,istate,ibasis,ispin,ileft,iright,iproc,iside,iirecv
    integer :: info
    
    logical :: ifhop
    
    real(kind=r8) :: tmpratio(4)
    ! MPI flag
    integer :: status(MPI_STATUS_SIZE),pppVrecvrequest,hoprecvrequest
    integer :: ierr
    
!============================================================
! allocate workspace
    ! store nosymmetry coeff
    allocate(coeffnosymm(nosymmdim*smadim))

    ! set the sparse matrix dim
    pppVnelement=CEILING(DBLE(16*isubM*isubM)/pppmatratio)
    hopnelement=CEILING(DBLE(16*isubM*isubM)/hopmatratio)
    LRoutnelement=CEILING(DBLE(16*isubM*isubM)/LRoutratio)
    npppVbuf=CEILING(DBLE(16*isubM*isubM)/pppmatratio)
    npppVmidmat=CEILING(DBLE(16*isubM*isubM)/pppmatratio)
    nhopbuf=CEILING(DBLE(16*isubM*isubM)/hopmatratio)
    nhopmidmat=CEILING(DBLE(16*isubM*isubM)/hopmatratio)

    if(myid/=0) then
        do ileft=1,nleft+1,1
            if(myid==orbid1(ileft,1)) then

                if( .not. allocated(LRcoeffin)) then
                    ! transform the 1-array to 4M*4M form 
                    allocate(LRcoeffin(nosymmdim,smadim))   ! coeff to LR format
                    allocate(LRcoeffincol(nosymmdim,smadim))   
                    allocate(LRcoeffinrow(4*iLrealdim+1,smadim))   
                end if
                
                if( .not. allocated(LRcoeffout)) then
                    allocate(LRcoeffout(LRoutnelement,smadim))  ! newcoeff to LR format
                    allocate(LRcoeffoutcol(LRoutnelement,smadim))  
                    allocate(LRcoeffoutrow(4*iLrealdim+1,smadim)) 
                end if
            end if
        end do

        do ileft=1,nleft+1,1
        do iright=norbs,norbs-nright,-1
            if(bondlink(ileft,iright)==1) then
                if(myid==orbid1(ileft,1) .or. myid==orbid1(iright,1)) then
                    if(.not. allocated(hopmat)) then
                        allocate(hopmat(hopnelement,2)) ! store the hopping matrix
                        allocate(hopmatcol(hopnelement,2)) 
                        allocate(hopmatrow(4*iRrealdim+1,2)) 

                        hoppacksize=(hopnelement*12+4*(4*iRrealdim+1))*2 
                        allocate(hoppackbuf(hoppacksize)) ! packbuf to send hopping matrix
                    end if
                end if
            end if
        end do
        end do
        
        if(logic_PPP==1) then
            allocate(pppVmat(pppVnelement)) ! store the pppV matrix
            allocate(pppVmatcol(pppVnelement)) 
            allocate(pppVmatrow(4*iRrealdim+1)) 

            pppVpacksize=pppVnelement*12+4*(4*iRrealdim+1)
            allocate(pppVpackbuf(pppVpacksize)) ! packbuf to send the pppV matrix
        end if
    else  
        ! 0 process
        allocate(LRcoeffin(nosymmdim,smadim))   ! coeff to LR format
        allocate(LRcoeffincol(nosymmdim,smadim))   
        allocate(LRcoeffinrow(4*iLrealdim+1,smadim))   
        allocate(LRcoeffout(LRoutnelement,smadim))  ! newcoeff to LR format
        allocate(LRcoeffoutcol(LRoutnelement,smadim))  
        allocate(LRcoeffoutrow(4*iLrealdim+1,smadim)) 
    end if

!=================================================================================================
    
    ! unsymmetrize the coeff if needed
    if( myid==0 ) then
        ! if symmetry==.true. then transform the symmetry coeff to the unsymmetry coeff
        if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
            if(bigdim/=nsymmstate) then
                call master_print_message(bigdim,"In symmetry, op bigdim/=nsymmstate wrong!")
                stop
            end if
            do istate=1,smadim,1
                call symmetrizestate(nosymmdim,coeffnosymm(nosymmdim*(istate-1)+1:istate*nosymmdim),&
                    coeff(bigdim*(istate-1)+1:istate*bigdim),'u')
            end do
        else
            if(bigdim/=nosymmdim) then
                call master_print_message(bigdim,"op bigdim/=nosymmdim wrong bigdim:!")
                call master_print_message(nosymmdim,"op bigdim/=nosymmdim wrong! nosymmdim:")
                stop
            end if
            coeffnosymm=coeff
        end if
    end if

    ! send coeffnosymm to L space process
    do iproc=1,nprocs-1,1
        if(myid==0 .or. myid==iproc) then
            do ileft=1,nleft+1,1
                if(iproc==orbid1(ileft,1)) then
                    if(myid==0) then
                        call MPI_SEND(coeffnosymm,nosymmdim*smadim,MPI_real8,iproc,1,MPI_COMM_WORLD,ierr)
                    else if(myid==iproc) then
                        call MPI_RECV(coeffnosymm,nosymmdim*smadim,MPI_real8,0,1,MPI_COMM_WORLD,status,ierr)
                    end if
                    exit
                end if
            end do
        end if
    end do

!-----------------------------------------------------------------------------

! L space process do it
! to transform the 16M*M coeff to 4M*4M(L*R) format ; coeff(16M^2,n) to coeff(4M,4M,n) 
    if(allocated(LRcoeffin)) then
        ! transfer the coeffnosymm to matrix form
        do istate=1,smadim,1
            call coefftosparse(nosymmdim,&
                LRcoeffin(:,istate),LRcoeffincol(:,istate),LRcoeffinrow(:,istate),&
                nosymmdim,coeffnosymm((istate-1)*nosymmdim+1:istate*nosymmdim),&
                iLrealdim,iRrealdim,cap_goodbasis,cap_goodbasiscol(1:4*iRrealdim+1))
        end do
    end if

! L space process and 0 process initializaiton
! L space process have LRcoeffout and send to 0 process at last

    if(allocated(LRcoeffout)) then
        LRcoeffoutrow=1  ! define the LRcoeffout matrix is 0
    end if

!  calculate HL*1 and 1*HR 
    if(myid==0) then 
        allocate(buffmat(LRoutnelement))
        allocate(buffmatcol(LRoutnelement))
        allocate(buffmatrow(4*iLrealdim+1))
        
        do istate=1,smadim,1
            do iside=1,2,1
                if(iside==1) then ! HL*1
                    call mkl_dcsrmultcsr('N',0,8,4*iLrealdim,4*iLrealdim,4*iRrealdim, &
                        cap_Hbig(:,1),cap_Hbigcol(:,1),cap_Hbigrow(:,1), &
                        LRcoeffin(:,istate),LRcoeffincol(:,istate),LRcoeffinrow(:,istate), &
                        buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
                    call checkinfo(info)
                else ! 1*HR
                    call mkl_dcsrmultcsr('N',0,8,4*iLrealdim,4*iRrealdim,4*iRrealdim, &
                        LRcoeffin(:,istate),LRcoeffincol(:,istate),LRcoeffinrow(:,istate), &
                        cap_Hbig(:,2),cap_Hbigcol(:,2),cap_Hbigrow(:,2), &
                        buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
                    call checkinfo(info)
                end if
                ! add LRcoeffout and bufmat
                call SpMatAdd(4*iRrealdim,4*iLrealdim,LRcoeffout(:,istate),LRcoeffoutcol(:,istate),LRcoeffoutrow(:,istate),&
                'N',1.0D0,4*iRrealdim,4*iLrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
            end do
        end do

        deallocate(buffmat,buffmatcol,buffmatrow)
    end if

!------------------------------------------------
! vlr=Hlrl'r'*Cl'r'=sum(opt,l',r')=sum(Lopt,l') parity*Oll'*sum(Ropt,r') Orr'cl'r'
! the parallel schema is that 0 process bcast the coeff matrix to other process
! and 0 process gather the result
    if(myid/=0) then
    if(logic_PPP==1) then
        do ileft=1,nleft+1,1
            
            ! set the pppVmat to 0.0, because we use add algorithm
            ! pppVmat=0.0D0
            pppVmatrow=1
            do iright=norbs,norbs-nright,-1
                if(myid==orbid1(iright,1)) then
                    operaindex=orbid1(iright,2)*3
                    call SpMatAdd(4*iRrealdim,4*iRrealdim,pppVmat,pppVmatcol,pppVmatrow,&
                        'N',pppV(ileft,iright),4*iRrealdim,4*iRrealdim,&
                        cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),pppVnelement)
                end if
            end do

            if(myid==orbid1(ileft,1)) then
                
                ! the number of matrix should recv
                nrecv=0
                do iproc=1,nprocs-1
                    if(iproc/= myid) then
                        do iright=norbs,norbs-nright,-1
                            if(iproc==orbid1(iright,1)) then
                                nrecv=nrecv+1
                                exit
                            end if
                        end do
                    end if
                end do
                
                if(nrecv>0) then
                   ! call MPI_RECV(pppVpackbuf,pppVpacksize,MPI_PACKED,&
                   !      MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                    call MPI_IRECV(pppVpackbuf,pppVpacksize,MPI_PACKED,&
                         MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,pppVrecvrequest,ierr)
                end if
                
                allocate(pppVbufmat(npppVbuf))
                allocate(pppVbufcol(npppVbuf))
                allocate(pppVbufrow(4*iRrealdim+1))

                do iirecv=1,nrecv,1
                    call MPI_WAIT(pppVrecvrequest,status,ierr)

                    position1=0
                    call MPI_UNPACK(pppVpackbuf,pppVpacksize,position1,pppVbufrow(1),&
                        (4*iRrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
                    call MPI_UNPACK(pppVpackbuf,pppVpacksize,position1,pppVbufmat(1),&
                        pppVbufrow(4*iRrealdim+1)-1,MPI_real8,MPI_COMM_WORLD,ierr)
                    call MPI_UNPACK(pppVpackbuf,pppVpacksize,position1,pppVbufcol(1),&
                        pppVbufrow(4*iRrealdim+1)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
                    
                    if(iirecv<nrecv) then
                     !   call MPI_RECV(pppVpackbuf,pppVpacksize,MPI_PACKED,&
                     !        MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
                        call MPI_IRECV(pppVpackbuf,pppVpacksize,MPI_PACKED,&
                            MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,pppVrecvrequest,ierr)
                    end if

                    call SpMatAdd(4*iRrealdim,4*iRrealdim,pppVmat,pppVmatcol,pppVmatrow,&
                        'N',1.0D0,4*iRrealdim,4*iRrealdim,&
                        pppVbufmat,pppVbufcol,pppVbufrow,pppVnelement)
                end do

                deallocate(pppVbufmat,pppVbufcol,pppVbufrow)
                
                allocate(midmat(npppVmidmat))
                allocate(midmatcol(npppVmidmat))
                allocate(midmatrow(4*iLrealdim+1))

                allocate(buffmat(LRoutnelement))
                allocate(buffmatcol(LRoutnelement))
                allocate(buffmatrow(4*iLrealdim+1))
                
                write(*,*) ileft,pppVmatrow(4*iRrealdim+1)-1

                do istate=1,smadim,1
                    ! do sum_{iright \in Rspace} pppV(ileft,iright)*pppVmatbuf_iright_rr*Crr
                    call mkl_dcsrmultcsr('N',0,8,4*iLrealdim,4*iRrealdim,4*iRrealdim, &
                        LRcoeffin(:,istate),LRcoeffincol(:,istate),LRcoeffinrow(:,istate), &
                        pppVmat,pppVmatcol,pppVmatrow, &
                        midmat,midmatcol,midmatrow,npppVmidmat,info)
                    call checkinfo(info)

                    tmpratio(1)=DBLE(midmatrow(4*iLrealdim+1))/DBLE(16*isubM*isubM)
                    !call checkmem_OPmodMat("pppVnelement",tmpratio(1),1)

                    operaindex=orbid1(ileft,2)*3

                    ! buffmat is to save the intermediate matrix
                    call mkl_dcsrmultcsr('N',0,8,4*iLrealdim,4*iLrealdim,4*iRrealdim, &
                        cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex), &
                        midmat,midmatcol,midmatrow, &
                        buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
                    call checkinfo(info)

                    ! add LRcoeffout and buffmat
                    call SpMatAdd(4*iRrealdim,4*iLrealdim,LRcoeffout(:,istate),LRcoeffoutcol(:,istate),LRcoeffoutrow(:,istate),&
                    'N',1.0D0,4*iRrealdim,4*iLrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
                end do

                deallocate(midmat,midmatcol,midmatrow)
                deallocate(buffmat,buffmatcol,buffmatrow)

            else
                do iright=norbs,norbs-nright,-1
                    if(myid==orbid1(iright,1)) then
                        position1=0
                        call MPI_PACK(pppVmatrow(1),(4*iRrealdim+1),MPI_integer4,&
                            pppVpackbuf,pppVpacksize,position1,MPI_COMM_WORLD,ierr)
                        call MPI_PACK(pppVmat(1),pppVmatrow(4*iRrealdim+1)-1,MPI_real8,&
                            pppVpackbuf,pppVpacksize,position1,MPI_COMM_WORLD,ierr)
                        call MPI_PACK(pppVmatcol(1),pppVmatrow(4*iRrealdim+1)-1,MPI_integer4,&
                            pppVpackbuf,pppVpacksize,position1,MPI_COMM_WORLD,ierr)

                        !call MPI_ISEND(pppVpackbuf,position1,MPI_PACKED,orbid1(ileft,1),myid,MPI_COMM_WORLD,pppVsendrequest,ierr)
                        !call MPI_WAIT(pppVsendrequest,status,ierr)
                        call MPI_SEND(pppVpackbuf,position1,MPI_PACKED,orbid1(ileft,1),myid,MPI_COMM_WORLD,ierr)
                        exit
                     end if
                end do
            end if
        end do
    end if
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! hopping term
    if(myid/=0) then
    do ileft=1,nleft+1,1
        
        ! check if need hopping matrix
        ifhop=.false.
        do iright=norbs,norbs-nright,-1
            if(bondlink(ileft,iright)==1) then
                ifhop=.true.
                exit
            end if
        end do

        
        if(ifhop==.true.) then

            ! every node add hopping matrix together
            ! hopmat=0.0D0
            hopmatrow=1
            do iright=norbs,norbs-nright,-1
                if(bondlink(ileft,iright)==1 .and. myid==orbid1(iright,1)) then
                    do ispin=1,2,1
                        operaindex=orbid1(iright,2)*3-3+ispin
                        call SpMatAdd(4*iRrealdim,4*iRrealdim,hopmat(:,ispin),hopmatcol(:,ispin),&
                            hopmatrow(:,ispin),'N',t(ileft,iright),4*iRrealdim,4*iRrealdim,&
                            cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),hopnelement)
                    end do
                end if
            end do
        
            if(myid==orbid1(ileft,1)) then
                ! the number of matrix should recv
                nrecv=0
                do iproc=1,nprocs-1,1
                    if(iproc/=myid) then
                        do iright=norbs,norbs-nright,-1
                            if(bondlink(iright,ileft)==1 &
                                .and. iproc==orbid1(iright,1)) then
                                nrecv=nrecv+1
                                exit
                            end if
                        end do
                    end if
                end do
                
                if(nrecv>0) then
                    call MPI_IRECV(hoppackbuf,hoppacksize,MPI_PACKED,&
                        MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,hoprecvrequest,ierr)
                end if

                allocate(hopbufmat(nhopbuf,2))
                allocate(hopbufcol(nhopbuf,2))
                allocate(hopbufrow(4*iRrealdim+1,2))
                

                do iirecv=1,nrecv,1
                    call MPI_WAIT(hoprecvrequest,status,ierr)

                    position1=0
                    do ispin=1,2,1
                        call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopbufrow(1,ispin),&
                            4*iRrealdim+1,MPI_integer4,MPI_COMM_WORLD,ierr)
                        call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopbufmat(1,ispin),&
                            hopbufrow(4*iRrealdim+1,ispin)-1,MPI_real8,MPI_COMM_WORLD,ierr)
                        call MPI_UNPACK(hoppackbuf,hoppacksize,position1,hopbufcol(1,ispin),&
                            hopbufrow(4*iRrealdim+1,ispin)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
                    end do
                    if(iirecv<nrecv) then
                        call MPI_IRECV(hoppackbuf,hoppacksize,MPI_PACKED,&
                            MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,hoprecvrequest,ierr)
                    end if

                    do ispin=1,2,1
                        call SpMatAdd(4*iRrealdim,4*iRrealdim,hopmat(1,ispin),hopmatcol(1,ispin),hopmatrow(1,ispin),&
                            'N',1.0D0,4*iRrealdim,4*iRrealdim,&
                            hopbufmat(1,ispin),hopbufcol(1,ispin),hopbufrow(1,ispin),hopnelement)
                    end do
                end do

                deallocate(hopbufmat,hopbufcol,hopbufrow)

                allocate(midmat(nhopmidmat))
                allocate(midmatcol(nhopmidmat))
                allocate(midmatrow(4*iLrealdim+1))
                
                allocate(buffmat(LRoutnelement))
                allocate(buffmatcol(LRoutnelement))
                allocate(buffmatrow(4*iLrealdim+1))
                
                ! do sum_{iright \in Rspace} t(ileft,iright)*hopmat_iright_rr*Crr
                !---------------------------------------------------------
                ! the +1 -1 phase added to l' of hopmat
                allocate(phase(4*iLrealdim))

                do il=1,4*iLrealdim,1
                    phase(il)=(-1.0D0)**(mod(cap_quantabigL(il,1),2))
                end do
                
                operaindex=3*orbid1(ileft,2)
                do istate=1,smadim,1
                    do ispin=1,4,1
                        if(ispin<=2) then
                        ! firstly needed to transfer from a(+) to a, so 'T' and 'T' is 'N', (ni-1)^+=(ni-1)
                        ! ispin=1 a up,ispin=2 a down,ispin=3 n,ispin=4 a+ up,ispin=5 a+ down;
                            call mkl_dcsrmultcsr('N',0,8,4*iLrealdim,4*iRrealdim,4*iRrealdim, &
                                    LRcoeffin(:,istate),LRcoeffincol(:,istate),LRcoeffinrow(:,istate), &
                                    hopmat(1,ispin),hopmatcol(1,ispin),hopmatrow(1,ispin), &
                                    midmat,midmatcol,midmatrow,nhopmidmat,info)
                            call checkinfo(info)
                        else
                            call  SpMMtoSp('N','T',4*iLrealdim,4*iRrealdim,4*iRrealdim,4*iRrealdim,4*iLrealdim,&
                                    LRcoeffin(:,istate),LRcoeffincol(:,istate),LRcoeffinrow(:,istate), &
                                    hopmat(1,ispin-2),hopmatcol(1,ispin-2),hopmatrow(1,ispin-2), &
                                    midmat,midmatcol,midmatrow,nhopmidmat)
                        end if
                        tmpratio(1)=DBLE(midmatrow(4*iLrealdim+1))/DBLE(16*isubM*isubM)
                        ! call checkmem_OPmodMat("hopnelement",tmpratio(1:4),4)

                        do il=1,4*iLrealdim,1
                        do ie=midmatrow(il),midmatrow(il+1)-1,1  ! ie=ielement
                            midmat(ie)=midmat(ie)*phase(il)
                        end do
                        end do

                        if(ispin>2) then
                            midmat(1:midmatrow(4*iLrealdim+1)-1)=midmat(1:midmatrow(4*iLrealdim+1)-1)*(-1.0D0)
                        end if
                    
                        !ispin<=2 al^+*ar,ispin>2 al*ar^(+) 
                        if(ispin<=2) then
                            call mkl_dcsrmultcsr('N',0,8,4*iLrealdim,4*iLrealdim,4*iRrealdim, &
                                cap_big(:,operaindex-3+ispin),cap_bigcol(:,operaindex-3+ispin),cap_bigrow(:,operaindex-3+ispin), &
                                midmat,midmatcol,midmatrow, &
                                buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
                            call checkinfo(info)
                        else
                            call mkl_dcsrmultcsr('T',0,8,4*iLrealdim,4*iLrealdim,4*iRrealdim, &
                                cap_big(:,operaindex-5+ispin),cap_bigcol(:,operaindex-5+ispin),cap_bigrow(:,operaindex-5+ispin), &
                                midmat,midmatcol,midmatrow, &
                                buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
                            call checkinfo(info)
                        end if

                        call SpMatAdd(4*iRrealdim,4*iLrealdim,LRcoeffout(:,istate),LRcoeffoutcol(:,istate),LRcoeffoutrow(:,istate),&
                        'N',1.0D0,4*iRrealdim,4*iLrealdim,buffmat,buffmatcol,buffmatrow,LRoutnelement)
                    end do
                end do

                deallocate(midmat,midmatcol,midmatrow)
                deallocate(buffmat,buffmatcol,buffmatrow)
                deallocate(phase)

            else 
                do iright=norbs,norbs-nright,-1
                    if(bondlink(ileft,iright)==1 .and. &
                        myid==orbid1(iright,1)) then
                        
                        position1=0
                        do ispin=1,2,1
                            call MPI_PACK(hopmatrow(1,ispin),(4*iRrealdim+1),MPI_integer4,&
                                hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
                            call MPI_PACK(hopmat(1,ispin),hopmatrow(4*iRrealdim+1,ispin)-1,MPI_real8,&
                                hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
                            call MPI_PACK(hopmatcol(1,ispin),hopmatrow(4*iRrealdim+1,ispin)-1,MPI_integer4,&
                                hoppackbuf,hoppacksize,position1,MPI_COMM_WORLD,ierr)
                        end do
                        !call MPI_ISEND(hoppackbuf,position1,MPI_PACKED,orbid1(ileft,1),myid,MPI_COMM_WORLD,pppVsendrequest,ierr)
                        !call MPI_WAIT(hopsendrequest(j),status,ierr)
                        call MPI_SEND(hoppackbuf,position1,MPI_PACKED,orbid1(ileft,1),myid,MPI_COMM_WORLD,ierr)
                        exit
                    end if
                end do
            end if
        end if
    end do
    end if

    ! every process transfer LRcoeffout to coeffnosymm
    if(allocated(LRcoeffout)) then
        
        do istate=1,smadim,1
            tmpratio(1)=DBLE(LRcoeffoutrow(4*iLrealdim+1,istate))/DBLE(16*isubM*isubM)
           ! call checkmem_OPmodMat("LRoutnelement",tmpratio(1),1)
        end do

        !m=0
        !do istate=1,smadim,1
        !do ir=1,4*iRrealdim,1
        !do il=1,4*iLrealdim,1
        !    if((cap_quantabigL(il,1)+cap_quantabigR(ir,1)==nelecs) .and. &
        !        cap_quantabigL(il,2)+cap_quantabigR(ir,2)==totalSz) then
        !        m=m+1
        !        call SpMatIJ(4*iLrealdim,il,ir,LRcoeffout(:,istate),LRcoeffoutcol(:,istate),LRcoeffoutrow(:,istate),coeffnosymm(m))
        !    end if
        !end do
        !end do
        !end do
        
        m=0
        do istate=1,smadim,1
        do ibasis=1,nosymmdim,1
            m=m+1
            call SpMatIJ(4*iLrealdim,cap_goodbasis(ibasis,1),cap_goodbasis(ibasis,2),LRcoeffout(:,istate),LRcoeffoutcol(:,istate),LRcoeffoutrow(:,istate),coeffnosymm(m))
        end do
        end do
            

        if(m/=nosymmdim*smadim) then
            write(*,*) "========================"
            write(*,*) "m/=nosymmdim*smadim failed!",m,smadim
            write(*,*) "========================"
            stop
        end if
    else
        coeffnosymm=0.0D0   ! other process coeffnosymm does not sum up
    end if
    
    if(allocated(LRcoeffin)) deallocate(LRcoeffin,LRcoeffinrow,LRcoeffincol)
    if(allocated(LRcoeffout)) deallocate(LRcoeffout,LRcoeffoutrow,LRcoeffoutcol)
    if(allocated(pppVmat)) deallocate(pppVmat,pppVmatrow,pppVmatcol)
    if(allocated(hopmat)) deallocate(hopmat,hopmatrow,hopmatcol)
    if(allocated(hoppackbuf)) deallocate(hoppackbuf)
    if(allocated(pppVpackbuf)) deallocate(pppVpackbuf)
    
    if(myid==0) then
        allocate(coeffnosymmreduce(nosymmdim*smadim))
    end if

    call MPI_REDUCE(coeffnosymm,coeffnosymmreduce,nosymmdim*smadim,mpi_real8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    if(myid==0) then
        newcoeff=0.0D0
        if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
            do istate=1,smadim,1
                call symmetrizestate(nosymmdim,coeffnosymmreduce(nosymmdim*(istate-1)+1:istate*nosymmdim),&
                    newcoeff(bigdim*(istate-1)+1:istate*bigdim),'s')
            end do
        else
            call copy(coeffnosymmreduce,newcoeff)
        end if
    end if

    if(allocated(coeffnosymm)) deallocate(coeffnosymm)
    if(allocated(coeffnosymmreduce)) deallocate(coeffnosymmreduce)
    
return

end subroutine op


end module ABop

