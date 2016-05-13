module corrfunc_mod
    use exit_mod
    use kinds_mod
    use communicate
    use variables
    implicit none

contains    
!================================================
!================================================

subroutine corrfunc_driver(direction,lstate,rstate)
    use module_sparse
    use LinkOpExpec_mod
    USE MPI
    implicit none
    character(len=1),intent(in) :: direction
    integer :: lstate,rstate
    
    ! local
    real(kind=r8),allocatable :: ddcorr(:,:),sscorr(:,:),bbcorr(:,:,:),phase(:)
    character(len=20) :: ddfile,ssfile,bbfile
    integer :: orbstart,orbend,orbadd,&
        ileft,iright,iLrealdim,iRrealdim,isubM,&
        rproc,lproc,leadproc,&
        loperaindex,roperaindex
    integer :: i,ispin
    real(kind=r8) :: midratio
    integer :: ierr

    call master_print_message("enter correlation function subroutine")

    if(logic_bondorder==0) then
        call master_print_message(logic_bondorder,"logic_bondorder==0")
        call exit_DMRG(sigAbort,"logic_bondorder==0")
    end if
    
    if(myid==0) then
        allocate(ddcorr(norbs,norbs))
        allocate(sscorr(norbs,norbs))
        allocate(bbcorr(norbs,norbs,2))
    end if
    ddfile=trim(direction)//"ddcorr.out"
    ssfile=trim(direction)//"sscorr.out"
    bbfile=trim(direction)//"bbcorr.out"

    if(direction=="l") then
        orbstart=norbs-nright
        orbend=norbs
        orbadd=nleft+1
    else if(direction=="r") then
        orbstart=1
        orbend=nleft+1
        orbadd=norbs-nright
    end if
    
    midratio=1.0D0

    allocate(phase(4*Lrealdim))
    do i=1,4*Lrealdim,1
        phase(i)=(-1.0D0)**(mod(quantabigL(i,1),2))
    end do

    ! density-density correlation
    do i=orbstart,orbend,1
        
        if(direction=="r") then
            ileft=i
            iright=norbs-nright
        else
            ileft=nleft+1
            iright=i
        end if
        iLrealdim=Lrealdim
        iRrealdim=Rrealdim
        isubM=subM

        lproc=orbid2(ileft,ileft,1)
        rproc=orbid2(iright,iright,1)
        loperaindex=orbid2(ileft,ileft,2)*2-1
        roperaindex=orbid2(iright,iright,2)*2-1
        if(direction=="r") then
            leadproc=lproc
        else
            leadproc=rproc
        end if
        
        if(myid==lproc .or. myid==rproc .or. myid==0) then
            call LinkOpExpec(ddcorr(ileft,iright),'N','N',lproc,rproc,leadproc,&
                    ileft,iright,iLrealdim,iRrealdim,isubM,lstate,rstate,&
                    loperaindex,roperaindex,operamatbig2,bigcolindex2,bigrowindex2,&
                    coeffIF,coeffIfcolindex,coeffIFrowindex,bigratio2,midratio,.false.,phase)
        end if

        ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! spin-spin correlation
        
        lproc=orbid2(ileft,ileft,1)
        rproc=orbid2(iright,iright,1)
        loperaindex=orbid2(ileft,ileft,2)*2
        roperaindex=orbid2(iright,iright,2)*2
        if(direction=="r") then
            leadproc=lproc
        else
            leadproc=rproc
        end if

        if(myid==lproc .or. myid==rproc .or. myid==0) then
            call LinkOpExpec(sscorr(ileft,iright),'N','N',lproc,rproc,leadproc,&
                    ileft,iright,iLrealdim,iRrealdim,isubM,lstate,rstate,&
                    loperaindex,roperaindex,operamatbig2,bigcolindex2,bigrowindex2,&
                    coeffIF,coeffIfcolindex,coeffIFrowindex,bigratio2,midratio,.false.,phase)
        end if
        ! call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! long range bondorder correlation/first order reduced density matrix
        do ispin=1,2,1

            lproc=orbid1(ileft,1)
            rproc=orbid1(iright,1)
            loperaindex=orbid1(ileft,2)*3-3+ispin
            roperaindex=orbid1(iright,2)*3-3+ispin
            if(direction=="r") then
                leadproc=lproc
            else
                leadproc=rproc
            end if

            if(myid==lproc .or. myid==rproc .or. myid==0) then
                call LinkOpExpec(bbcorr(ileft,iright,ispin),'N','T',lproc,rproc,leadproc,&
                        ileft,iright,iLrealdim,iRrealdim,isubM,lstate,rstate,&
                        loperaindex,roperaindex,operamatbig1,bigcolindex1,bigrowindex1,&
                        coeffIF,coeffIfcolindex,coeffIFrowindex,bigratio1,midratio,.true.,phase)
            end if
        !   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do 
    end do
    
    if(myid==0) then
        open(unit=880,file=ddfile,status="replace")
        open(unit=881,file=ssfile,status="replace")
        open(unit=882,file=bbfile,status="replace")
        do i=orbstart,orbend,1
            if(direction=="l") then
                write(880,*) i,ddcorr(orbadd,i)
                write(881,*) i,sscorr(orbadd,i)
                write(882,*) i,bbcorr(orbadd,i,1:2)
            else
                write(880,*) i,ddcorr(i,orbadd)
                write(881,*) i,sscorr(i,orbadd)
                write(882,*) i,bbcorr(i,orbadd,1:2)
            end if
        end do
        close(880)
        close(881)
        close(882)
        deallocate(ddcorr,sscorr,bbcorr)
    end if
    deallocate(phase)

    return
end subroutine corrfunc_driver

!================================================
!================================================
end module corrfunc_mod




