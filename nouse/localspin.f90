module LocalSpin_mod
! this module calculate the local spin
! the algorithm follows I.Mayer CPL 478(2009) 323-326

    use module_sparse
    use variables
    use communicate

    implicit none
    
    real(kind=r8),allocatable ::  localspin0(:,:,:)

    contains
!===========================================================================
!===========================================================================

subroutine init_LocalSpinmat(orbindex)
! initiate the on site localspin matrix
    use onesitematrix

    implicit none

    integer :: orbindex
    ! local
    integer :: operaindex

    if(myid==orbid3(orbindex,orbindex,1)) then
        call ConstructOnesiteMatrix(orbindex)
        operaindex=orbid3(orbindex,orbindex,2)
        
        operamatsma3(1:4,operaindex-1)=onesitemat(:,10)
        smacolindex3(1:4,operaindex-1)=osmcolindex(:,10)
        smarowindex3(1:5,operaindex-1)=osmrowindex(:,10)
        smarowindex3(5:subM+1,operaindex-1)=smarowindex3(5,operaindex-1)
        
        operamatsma3(1:4,operaindex)=onesitemat(:,11)
        smacolindex3(1:4,operaindex)=osmcolindex(:,11)
        smarowindex3(1:5,operaindex)=osmrowindex(:,11)
        smarowindex3(5:subM+1,operaindex)=smarowindex3(5,operaindex)
    end if
return
end subroutine init_LocalSpinmat

!===========================================================================
!===========================================================================

subroutine  LocalSpin
! this subroutine calculate the Local spin in the last step
! <SA*SA>,<SA*SB>,SA^+=SA
! <SA*SB>=<SA*SB>^+=<SB^+SA^+>=<SB*SA>
    
    use mpi
    implicit none
    
    integer :: error,ierr
    integer :: i,j,k,tmp
    real(kind=r8) :: totalspin(nstate),halfspin(2,nstate) !,pairspin(norbs,norbs,nstate)
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call master_print_message("enter LocalSpin subroutine")
    
    if(myid==0) then
        allocate(localspin0(norbs,norbs,nstate),stat=error)
        if(error/=0) stop
        localspin0=0.0D0
    end if
    
    call Calc_Localspin_link
    call Calc_Localspin_subspace('L')
    call Calc_Localspin_subspace('R')

    if(myid==0) then
        totalspin=0.0D0
        halfspin=0.0D0
        do k=1,nstate,1
            do i=1,norbs,1
            do j=1,norbs,1
                totalspin(k)=totalspin(k)+localspin0(i,j,k)
                if(i<=(norbs+1)/2 .and. j<=(norbs+1)/2) then
                    halfspin(1,k)=halfspin(1,k)+localspin0(i,j,k)
                else if(i>(norbs+1)/2 .and. j>(norbs+1)/2) then
                    halfspin(2,k)=halfspin(2,k)+localspin0(i,j,k)
                end if
            end do
            end do
            write(*,*) "Total Spin of State",k,"equals",totalspin(k)
            write(*,*) "Half Spin of State",k,"equals",halfspin(:,k)
        end do

        open(unit=170,file="localspin.out",form="unformatted",status="replace")
        write(170) norbs,nstate
        write(170) localspin0
        close(170)
        
        ! pair spin
!       pairspin=0.0D0
!       do k=1,nstate,1
!           do i=1,norbs-1,1
!               pairspin(i,i+1,k)=localspin0(i,i,k)+localspin0(i+1,i+1,k)+localspin0(i,i+1,k)+localspin0(i+1,i,k)
!               write(*,*) i,i+1,k,pairspin(i,i+1,k)
!           end do
!       end do
        call  analysis_cutoffspin 
    end if

return
end subroutine LocalSpin

!===========================================================================
!===========================================================================

subroutine  analysis_cutoffspin 
! this subroutine calculate the spin from 1 to j or j to norbs
! used in polyene big hubbard U ,spin chain
    implicit none
    integer :: k,i,j,istate
    real(kind=r8) :: cutoffspin
    
    open(unit=152,file="cutoffspin.out",status="replace")
    do istate=1,nstate,1
        do i=1,norbs,1
            cutoffspin=0.0D0
            do j=1,i,1
            do k=1,i,1
                cutoffspin=cutoffspin+localspin0(k,j,istate)
            end do
            end do
            write(152,*) i,cutoffspin
        end do
    end do
    close(152)

return
end subroutine analysis_cutoffspin

!===========================================================================
!===========================================================================

subroutine Calc_Localspin_subspace(domain)
! calculate operator in the L/R subspace i,j<=nleft+1,or i,j>=norbs-nright
! <R|<L|CLR OP CL'R'|L'>|R'>
    use exit_mod
    use mpi
    use mathlib
    implicit none
    include "mkl_spblas.fi"
    
    character(len=1) :: domain
    
    ! local
    integer :: orbstart,orbend,operaindex2,operaindex3
    integer :: i,j,istate
    integer :: nmid
    integer :: error,ierr,info
    integer :: status(MPI_STATUS_SIZE) ! MPI flag
    real(kind=r8) :: Spratio
    real(kind=r8)   :: ilocalspin(nstate)
    real(kind=r8),allocatable :: midmat(:),Spaddmat(:)
    integer(kind=i4),allocatable :: midcolindex(:),midrowindex(:),coeffIFrowindexdummy(:,:),&
                 Spaddcolindex(:),Spaddrowindex(:)

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

    if(myid/=0) then
        Spratio=bigratio3/5.0
        Spratio=MAX(Spratio,1.0D0)
        nmid=CEILING(DBLE(16*subM*subM)/Spratio)
        allocate(midmat(nmid),stat=error)
        if(error/=0) stop
        allocate(midcolindex(nmid),stat=error)
        if(error/=0) stop
        allocate(midrowindex(4*subM+1),stat=error)
        if(error/=0) stop
        allocate(coeffIFrowindexdummy(4*subM+1,nstate),stat=error)
        if(error/=0) stop
        
        ! store the sparse matrix add S(i,j) pair
        allocate(Spaddmat(bigdim3),stat=error)
        if(error/=0) stop
        allocate(Spaddcolindex(bigdim3),stat=error)
        if(error/=0) stop
        allocate(Spaddrowindex(4*subM+1),stat=error)
        if(error/=0) stop
    end if

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

    ! <SA*SA>
    ! two operator matrix => no phase
    do i=orbstart,orbend,1
        if(myid==orbid3(i,i,1)) then
            
            operaindex3=orbid3(i,i,2)
            do istate=1,nstate,1
                !<SA*SA>=3/4*(niup-nidown)^2
                call SpMMtoSp('N','N',4*subM,4*subM,4*subM,4*subM,4*subM,&
                    operamatbig3(:,operaindex3),bigcolindex3(:,operaindex3),bigrowindex3(:,operaindex3), &
                    coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
                    midmat,midcolindex,midrowindex,nmid)
                ! trace(CLR*QLR)
                call SpMMtrace('T',4*subM,&
                    coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
                    midmat,midcolindex,midrowindex,ilocalspin(istate))
                ilocalspin(istate)=0.75D0*ilocalspin(istate)
            end do
            call MPI_SEND(ilocalspin,nstate,mpi_real8,0,orbid3(i,i,2),MPI_COMM_WORLD,ierr)
        
        else if(myid==0) then
            call MPI_RECV(ilocalspin,nstate,mpi_real8,orbid3(i,i,1),orbid3(i,i,2),MPI_COMM_WORLD,status,ierr)
            localspin0(i,i,:)=ilocalspin
        end if
    end do

    ! <SA*SB>
    do i=orbstart,orbend,1
    do j=i+1,orbend,1
        if(myid==orbid3(i,j,1)) then

            operaindex3=orbid3(i,j,2)

            do istate=1,nstate,1
                ! ai^+down*aiup*aj^+up*aj^down+1/4*(niup-nidown)(njup-njdown)
                call mkl_dcsradd('N',0,0,4*subM,4*subM, &
                    operamatbig3(:,operaindex3-1),bigcolindex3(:,operaindex3-1),bigrowindex3(:,operaindex3-1), &
                    0.25D0,operamatbig3(:,operaindex3),bigcolindex3(:,operaindex3),bigrowindex3(:,operaindex3), &
                    Spaddmat,Spaddcolindex,Spaddrowindex,bigdim3,info)
                call checkinfo(info)
                
                call SpMMtoSp('N','N',4*subM,4*subM,4*subM,4*subM,4*subM,&
                    Spaddmat,Spaddcolindex,Spaddrowindex,&
                    coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
                    midmat,midcolindex,midrowindex,nmid)
                ! trace(CLR*QLR)
                call SpMMtrace('T',4*subM,&
                    coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindexdummy(:,istate),&
                    midmat,midcolindex,midrowindex,ilocalspin(istate))
            end do
            call MPI_SEND(ilocalspin,nstate,mpi_real8,0,orbid3(i,j,2),MPI_COMM_WORLD,ierr)
        
        else if(myid==0) then
            call MPI_RECV(ilocalspin,nstate,mpi_real8,orbid3(i,j,1),orbid3(i,j,2),MPI_COMM_WORLD,status,ierr)
            localspin0(i,j,:)=ilocalspin
            localspin0(j,i,:)=ilocalspin
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

    if(myid/=0) deallocate(midmat,midcolindex,midrowindex,coeffIFrowindexdummy,Spaddmat,Spaddcolindex,Spaddrowindex)

return

end subroutine Calc_LocalSpin_subspace

!===========================================================================
!===========================================================================

subroutine Calc_Localspin_link
! this subroutine calculate localspin belong to different subspace
! i<=nleft+1,j>=norbs-nright 
    use mathlib
    use mpi
    implicit none
    include "mkl_spblas.fi"

    integer :: midnelement
    real(kind=r8),allocatable :: Spaddmat(:),midmat(:),midmat2(:),Ropmat(:,:)
    integer(kind=i4),allocatable :: &
    Spaddcolindex(:),Spaddrowindex(:),&
    midmatcol(:),midmatrow(:),&
    midmatcol2(:),midmatrow2(:),&
    Ropmatcol(:,:),Ropmatrow(:,:)

    integer :: operaindex2,operaindex3,nonzero
    real(kind=r8) :: ilocalspin(nstate)
    integer :: i,j,k,l,m
    
    character(len=1),allocatable :: packbuf(:)
    integer :: position1,packsize
    integer :: touched(nprocs-1),ntouched
    logical :: ifsend
    
    integer :: error,info
    integer :: status(MPI_STATUS_SIZE),sendrequest(nprocs-1)
    integer :: ierr

    midnelement=CEILING(DBLE(16*subM*subM)/pppmatratio)

    do i=1,nleft+1,1
    do j=norbs,norbs-nright,-1
        
        if(myid==orbid3(i,i,1) .or. myid==orbid3(j,j,1)) then
            if(.not. allocated(packbuf)) then
                packsize=(bigdim2*12+4*(4*subM+1))+(bigdim3*12+4*(4*subM+1))
                allocate(packbuf(packsize),stat=error) ! packbuf to send
                if(error/=0) stop
            end if
        end if

        if(myid==orbid1(i,1)) then
            if(.not. allocated(midmat)) then
                
                allocate(Ropmat(bigdim2,2),stat=error) ! store the intermediate matrix
                if(error/=0) stop
                allocate(Ropmatcol(bigdim2,2),stat=error) 
                if(error/=0) stop
                allocate(Ropmatrow(4*subM+1,2),stat=error) 
                if(error/=0) stop

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

                allocate(Spaddmat(midnelement),stat=error)
                if(error/=0) stop
                allocate(Spaddcolindex(midnelement),stat=error)
                if(error/=0) stop
                allocate(Spaddrowindex(4*subM+1),stat=error)
                if(error/=0) stop
            end if
        end if
    end do
    end do

    do i=norbs,norbs-nright,-1
        if(myid==orbid3(i,i,1)) then
            if(orbid3(i,i,1)/=orbid2(i,i,1)) then
                write(*,*) "orbid3(i,i,1)/=orbid2(i,i,1)",orbid3(i,i,1),orbid2(i,i,1)
                stop
            end if

            ! pack the L space needed matrix
            ! aq^down*aqup,(niup-nidown)
            position1=0
            
            operaindex3=orbid3(i,i,2)-1
            call MPI_PACK(bigrowindex3(1,operaindex3),(4*subM+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            nonzero=bigrowindex3(4*subM+1,operaindex3)-1
            call MPI_PACK(operamatbig3(1,operaindex3),nonzero,MPI_REAL8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            call MPI_PACK(bigcolindex3(1,operaindex3),nonzero,MPI_INTEGER4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)

            operaindex2=orbid2(i,i,2)*2
            call MPI_PACK(bigrowindex2(1,operaindex2),(4*subM+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            nonzero=bigrowindex2(4*subM+1,operaindex2)-1
            call MPI_PACK(operamatbig2(1,operaindex2),nonzero,MPI_REAL8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            call MPI_PACK(bigcolindex2(1,operaindex2),nonzero,MPI_INTEGER4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            
            ! send the packbuf
            touched=0
            ntouched=0
            do l=1,nleft+1,1
                if(orbid3(l,l,1)/=myid) then
                    ifsend=.false.
                    do m=1,ntouched,1
                        if(orbid3(l,l,1)==touched(m)) then
                            ifsend=.true.
                            exit
                        end if
                    end do
                    if(ifsend==.false.) then
                        ntouched=ntouched+1
                        touched(ntouched)=orbid3(l,l,1)
                        call MPI_SEND(packbuf,position1,MPI_PACKED,orbid3(l,l,1),i,MPI_COMM_WORLD,ierr)
                    !   call MPI_ISEND(packbuf,position1,MPI_PACKED,orbid3(l,l,1),i,MPI_COMM_WORLD,sendrequest(ntouched),ierr)
                    end if
                end if
            end do
        end if

        ! recv packbuf
        do l=1,nleft+1,1
            if(myid==orbid3(l,l,1)) then
                if(myid/=orbid3(i,i,1)) then
                    call MPI_RECV(packbuf,packsize,MPI_PACKED,orbid3(i,i,1),i,MPI_COMM_WORLD,status,ierr)
                end if
                position1=0
                do j=1,2,1
                    call MPI_UNPACK(packbuf,packsize,position1,Ropmatrow(1,j),(4*subM+1),MPI_integer4,MPI_COMM_WORLD,ierr)
                    nonzero=Ropmatrow(4*subM+1,j)-1
                    call MPI_UNPACK(packbuf,packsize,position1,Ropmat(1,j),nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
                    call MPI_UNPACK(packbuf,packsize,position1,Ropmatcol(1,j),nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
                end do
                exit  ! only recv once
            end if
        end do

        do l=1,nleft+1,1
            if(myid==orbid3(l,l,1)) then
                if(orbid3(l,l,1)/=orbid2(l,l,1)) then
                    write(*,*) "orbid3(l,l,1)/=orbid2(l,l,1)",orbid3(l,l,1),orbid2(l,l,1)
                    stop
                end if

                do j=1,nstate,1
                    
                ! ai^+down*aiup*aj^+up*ajdown
                    ! firstly needed to transfer R space from aj^+down*ajup to aj^+up*ajdown, so 'T' and 'T' is 'N'
                    ! ORR'*CL'R' no phase
                    call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
                            coeffIF(:,j),coeffIFcolindex(:,j),coeffIFrowindex(:,j), &
                            Ropmat(:,1),Ropmatcol(:,1),Ropmatrow(:,1), &
                            midmat,midmatcol,midmatrow,midnelement,info)
                    call checkinfo(info)
                
                    ! OLL'*(OC)L'R
                    operaindex3=orbid3(l,l,2)-1
                    call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
                            operamatbig3(:,operaindex3),bigcolindex3(:,operaindex3),bigrowindex3(:,operaindex3), &
                            midmat,midmatcol,midmatrow,&
                            Spaddmat,Spaddcolindex,Spaddrowindex,midnelement,info)
                    call checkinfo(info)

                ! (niup-nidown)*(njup-njdown)
                    call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
                            coeffIF(:,j),coeffIFcolindex(:,j),coeffIFrowindex(:,j), &
                            Ropmat(:,2),Ropmatcol(:,2),Ropmatrow(:,2), &
                            midmat,midmatcol,midmatrow,midnelement,info)
                    call checkinfo(info)
                    operaindex2=orbid2(l,l,2)*2
                    call mkl_dcsrmultcsr('N',0,8,4*subM,4*subM,4*subM, &
                            operamatbig2(:,operaindex2),bigcolindex2(:,operaindex2),bigrowindex2(:,operaindex2), &
                            midmat,midmatcol,midmatrow,&
                            midmat2,midmatcol2,midmatrow2,midnelement,info)
                    call checkinfo(info)

                ! add the two conponent
                    call SpMatAdd(4*subM,4*subM,Spaddmat,Spaddcolindex,Spaddrowindex,'N',0.25D0,&
                            4*subM,4*subM,midmat2,midmatcol2,midmatrow2,midnelement)
            
                ! trace(CLR*OLR)
                    call SpMMtrace('T',4*subM, & 
                            coeffIF(:,j),coeffIFcolindex(:,j),coeffIFrowindex(:,j), &
                            Spaddmat,Spaddcolindex,Spaddrowindex,ilocalspin(j))
                end do
                call MPI_SEND(ilocalspin,nstate,mpi_real8,0,orbid3(l,l,2),MPI_COMM_WORLD,ierr)
            else if(myid==0) then
                call MPI_RECV(ilocalspin,nstate,mpi_real8,orbid3(l,l,1),orbid3(l,l,2),MPI_COMM_WORLD,status,ierr)
                localspin0(i,l,:)=ilocalspin(:)
                localspin0(l,i,:)=ilocalspin(:)
            end if
        end do

        ! confirm that the packbuf can be used again without problem
    !   if(myid==orbid3(i,i,1)) then
    !       do j=1,ntouched,1
    !           call MPI_WAIT(sendrequest(j),status,ierr)
    !       end do
    !   end if
    end do

    if(allocated(packbuf)) deallocate(packbuf)
    if(allocated(midmat)) deallocate(midmat,midmatcol,midmatrow)
    if(allocated(midmat2)) deallocate(midmat2,midmatcol2,midmatrow2)
    if(allocated(Spaddmat)) deallocate(Spaddmat,Spaddcolindex,Spaddrowindex)
    if(allocated(Ropmat)) deallocate(Ropmat,Ropmatcol,Ropmatrow)

return
end subroutine Calc_LocalSpin_link

!===========================================================================
!===========================================================================
end module LocalSpin_mod
