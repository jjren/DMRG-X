Subroutine Enviro_Big(domain)
! this subroutine is to contruct the enviroment R+sigmaR/L+sigmaL space
! all the operamatbig(does not need operamatsma) read from disc

    USE mpi
    USE variables
    use communicate
    use exit_mod
    use module_sparse

    implicit none

    character(len=1) :: domain

    integer :: i,j,k
    integer :: orbnow,orbstart,orbend,orbref,nsuborbs
    integer :: reclength,operaindex,recindex,Hindex
    character(len=50) :: filename,Hfilename,Hcolfilename,Hrowfilename,symmfilename,quantafilename,&
        filenamep,Hfilenamep,Hcolfilenamep,Hrowfilenamep,quantafilenamep,LRdimfilename
    logical :: alive,ifmiddle
    integer :: count1,dim1,nonzero

    integer :: thefile,status(MPI_STATUS_SIZE)  ! MPI flag
    integer(kind=MPI_OFFSET_KIND) ::offset
    integer :: ierr

    ! reclength is the length of direct io
    ! ifort use 4 byte as 1 by default
    call master_print_message("enter in subroutine Enviro_Big")

!   if(myid==1) then
!       write(*,*) "subM=",subM
!       write(*,*) "Lrealdim=",Lrealdim
!       write(*,*) "Rrealdim=",Rrealdim
!   end if
    
    ifmiddle=.false.
    if(domain=='L') then
        orbnow=nleft+1
        write(filename,'(i5.5,a8)') orbnow,'left.tmp'
        orbstart=1
        orbend=nleft+1
        orbref=1
        nsuborbs=nleft+1
        dim1=Lrealdim
    !   if(nleft==(norbs+1)/2-1 .and. mode/='r') then
        if(nleft==(norbs+1)/2-1) then
            ifmiddle=.true.
        end if
    else if (domain=='R') then
        orbnow=norbs-nright
        write(filename,'(i5.5,a9)') orbnow,'right.tmp'
        orbstart=norbs-nright
        orbend=norbs
        orbref=norbs
        nsuborbs=nright+1
        dim1=Rrealdim
    !   if(nright==norbs/2-1 .and. mode/='r') then
        if(nright==norbs/2-1) then
            ifmiddle=.true.
        end if
    else
        call exit_DMRG(sigAbort,"domain/=L .and. domain/='R' failed!")
    end if
    
    ! L/R + sigmaL/R space

    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,thefile,ierr)
    if(myid/=0) then
        do i=orbstart,orbend,1
        if(myid==orbid1(i,1)) then
            operaindex=orbid1(i,2)
            ! get mat-> colindex->rowindex CSR format
            offset=abs(i-orbref)*(bigdim1*12+(4*subM+1)*4)*3
            do j=1,3,1
                if(logic_PPP==0 .and. j==3) exit
                call MPI_FILE_READ_AT(thefile,offset,bigrowindex1(1,3*operaindex-3+j),(4*subM+1),mpi_integer4,status,ierr)
                offset=offset+(4*subM+1)*4
                nonzero=bigrowindex1(4*dim1+1,3*operaindex-3+j)-1
                call MPI_FILE_READ_AT(thefile,offset,operamatbig1(1,3*operaindex-3+j),nonzero,mpi_real8,status,ierr)
                offset=offset+nonzero*8
                call MPI_FILE_READ_AT(thefile,offset,bigcolindex1(1,3*operaindex-3+j),nonzero,mpi_integer4,status,ierr)
                offset=offset+nonzero*4
            end do
        end if
        end do
    ! 0 process
    else if(myid==0) then
        if(domain=='L') then
            Hfilename="0-left.tmp"
            Hcolfilename="0-leftcol.tmp"
            Hrowfilename="0-leftrow.tmp"
            symmfilename="symmlink-left.tmp"
            quantafilename="quantabigL.tmp"
            Hindex=1
        else if(domain=='R') then
            Hfilename="0-right.tmp"
            Hcolfilename="0-rightcol.tmp"
            Hrowfilename="0-rightrow.tmp"
            symmfilename="symmlink-right.tmp"
            quantafilename="quantabigR.tmp"
            Hindex=2
        end if
        recindex=abs(orbnow-orbref)+1

!----------------open a binary file-------------
        ! rowindex
        reclength=subM*4+1
        inquire(file=trim(Hrowfilename),exist=alive)
        if(alive) then
            open(unit=112,file=trim(Hrowfilename),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(Hrowfilename),"doesn't exist!"
            stop
        end if
        ! sigmaR index =norbs is 1; norbs-1 is 2
        read(112,rec=recindex) Hbigrowindex(:,Hindex)
        close(112)

        nonzero=Hbigrowindex(4*dim1+1,Hindex)-1   ! nonzero element numbers in Hbig
        
        reclength=2*Hbigdim
        inquire(file=trim(Hfilename),exist=alive)
        if(alive) then
            open(unit=102,file=trim(Hfilename),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(Hfilename),"doesn't exist!"
            stop
        end if
! norbs is 1; norbs-1 is 2
        read(102,rec=recindex) Hbig(1:nonzero,Hindex)
        close(102)

        ! colindex
        reclength=Hbigdim
        inquire(file=trim(Hcolfilename),exist=alive)
        if(alive) then
            open(unit=110,file=trim(Hcolfilename),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(Hcolfilename),"doesn't exist!"
            stop
        end if
        ! sigmaR index =norbs is 1; norbs-1 is 2
        read(110,rec=recindex) Hbigcolindex(1:nonzero,Hindex)
        close(110)

!------------------------------read the symmlink---
        if(logic_spinreversal/=0) then
            reclength=2*subM
            inquire(file=trim(symmfilename),exist=alive)
            if(alive) then
                open(unit=104,file=trim(symmfilename),access="Direct",form="unformatted",recl=reclength,status="old")
            else
                write(*,*) trim(symmfilename),"doesn't exist!"
                stop
            end if
            read(104,rec=recindex) symmlinkbig(:,1,Hindex)
            close(104)
        end if
!------------------------------read the quantabigR/L(4*subM,2)---
! every process need to read quantabigR/L

        reclength=4*subM*2
        ! quantabigR/L is  integer(kind=4) so without *2
        inquire(file=trim(quantafilename),exist=alive)
        if(alive) then
            open(unit=108,file=trim(quantafilename),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(quantafilename),"doesn't exist!"
            stop
        end if
        if(domain=='L') then
            read(108,rec=recindex) quantabigL
        else if(domain=='R') then
            read(108,rec=recindex) quantabigR
        end if
        close(108)
    end if

    ! broadcast the quantabigR/L to other process
    if(domain=='L') then
        call MPI_BCAST(quantabigL,4*subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
    else if(domain=='R') then
        call MPI_BCAST(quantabigR,4*subM*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
    end if

    call MPI_FILE_CLOSE(thefile,ierr)

!============================================================
    ! read the bondorder matrix
    ! in the standard mode need not read the bomid.tmp
    ! only in the old restart mode
    ! this can be optimized after
    if(logic_bondorder/=0 .and. (nsuborbs==exactsite+1 .or. ifmiddle==.true.)) then
        if(ifmiddle==.true.) then
            write(filename,'(a1,a9)') domain,'bomid.tmp'
        else
            write(filename,'(a1,a6)') domain,'bo.tmp'
        end if
        call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,thefile,ierr)
        if(myid/=0) then
            count1=0
            do i=orbstart,orbend,1
            do j=i,orbend,1
                if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
                    count1=count1+1
                    if(myid==orbid2(i,j,1)) then
                        operaindex=orbid2(i,j,2)
                        offset=(count1-1)*(bigdim2*12+(4*subM+1)*4)*2
                        do k=1,2,1
                            ! read mat-> colindex->rowindex in CSR format
                            call MPI_FILE_READ_AT(thefile,offset,bigrowindex2(1,2*operaindex-2+k),(4*subM+1),mpi_integer4,status,ierr)
                            nonzero=bigrowindex2(4*dim1+1,2*operaindex-2+k)-1
                            offset=offset+(4*subM+1)*4
                            call MPI_FILE_READ_AT(thefile,offset,operamatbig2(1,2*operaindex-2+k),nonzero,mpi_real8,status,ierr)
                            offset=offset+nonzero*8
                            call MPI_FILE_READ_AT(thefile,offset,bigcolindex2(1,2*operaindex-2+k),nonzero,mpi_integer4,status,ierr)
                            offset=offset+nonzero*4
                        end do
                    end if
                end if
            end do
            end do
        end if
        call MPI_FILE_CLOSE(thefile,ierr)
    end if
!============================================================
    ! store local spin matrix(only store the exact site)
    if(logic_localspin==1 .and. (nsuborbs==exactsite+1 .or. ifmiddle==.true.)) then
        if(ifmiddle==.true.) then
            write(filename,'(a1,a16)') domain,'localspinmid.tmp'
        else
            write(filename,'(a1,a13)') domain,'localspin.tmp'
        end if
        call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_RDONLY,MPI_INFO_NULL,thefile,ierr)
        if(myid/=0) then
            count1=0
            do i=orbstart,orbend,1
            do j=i,orbend,1
                count1=count1+2  ! at now, the number of operators
                if(myid==orbid3(i,j,1)) then
                    offset=(count1-2)*(bigdim3*12+(4*subM+1)*4)
                    do k=1,2,1
                        operaindex=orbid3(i,j,2)-2+k
                        ! read rowindex-> mat-> colindex in CSR format
                        call MPI_FILE_READ_AT(thefile,offset,bigrowindex3(1,operaindex),(4*subM+1),mpi_integer4,status,ierr)
                        nonzero=bigrowindex3(4*dim1+1,operaindex)-1
                        offset=offset+(4*subM+1)*4
                        call MPI_FILE_READ_AT(thefile,offset,operamatbig3(1,operaindex),nonzero,mpi_real8,status,ierr)
                        offset=offset+nonzero*8
                        call MPI_FILE_READ_AT(thefile,offset,bigcolindex3(1,operaindex),nonzero,mpi_integer4,status,ierr)
                        offset=offset+nonzero*4
                    end do
                end if
            end do
            end do
        end if
        call MPI_FILE_CLOSE(thefile,ierr)
    end if

!=====================================================
    if(logic_perturbation/=0) then
    if(domain=='L') then
        write(filenamep,'(i5.5,a9)') orbnow,'leftp.tmp'
        LRdimfilename="Lrealdimp.tmp"
    else if (domain=='R') then
        write(filenamep,'(i5.5,a10)') orbnow,'rightp.tmp'
        LRdimfilename="Rrealdimp.tmp"
    end if
    ! read the Lrealdimp,Rrealdimp
    if(myid==0) then
        reclength=1
        inquire(file=LRdimfilename,exist=alive)
        if(alive) then
            open(unit=1001,file=LRdimfilename,access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) LRdimfilename,"doesn't exist!"
            stop
        end if
        recindex=abs(orbnow-orbref)+1
        read(1001,rec=recindex) dim1
        close(1001)
    end if
    ! bcast the L/Rrealdimp
    call MPI_BCAST(dim1,1,MPI_integer4,0,MPI_COMM_WORLD,ierr)
    if(domain=="L") then
        Lrealdimp=dim1
    else
        Rrealdimp=dim1
    end if

    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filenamep),MPI_MODE_RDONLY,MPI_INFO_NULL,thefile,ierr)
    if(myid/=0) then
        do i=orbstart,orbend,1
        if(myid==orbid1(i,1)) then
            operaindex=orbid1(i,2)
            ! get mat-> colindex->rowindex CSR format
            offset=abs(i-orbref)*(bigdim1p*12+(4*subMp+1)*4)*3
            do j=1,3,1
                call MPI_FILE_READ_AT(thefile,offset,bigrowindex1p(1,3*operaindex-3+j),(4*subMp+1),mpi_integer4,status,ierr)
                offset=offset+(4*subMp+1)*4
                nonzero=bigrowindex1p(4*dim1+1,3*operaindex-3+j)-1
            !   write(*,*) "nonzero=",dim1,nonzero,offset
                call MPI_FILE_READ_AT(thefile,offset,operamatbig1p(1,3*operaindex-3+j),nonzero,mpi_real8,status,ierr)
                offset=offset+nonzero*8
                call MPI_FILE_READ_AT(thefile,offset,bigcolindex1p(1,3*operaindex-3+j),nonzero,mpi_integer4,status,ierr)
                offset=offset+nonzero*4
            end do
        end if
        end do
    ! 0 process
    else if(myid==0) then
        if(domain=='L') then
            Hfilenamep="0-leftp.tmp"
            Hcolfilenamep="0-leftcolp.tmp"
            Hrowfilenamep="0-leftrowp.tmp"
            quantafilenamep="quantabigLp.tmp"
            Hindex=1
        else if(domain=='R') then
            Hfilenamep="0-rightp.tmp"
            Hcolfilenamep="0-rightcolp.tmp"
            Hrowfilenamep="0-rightrowp.tmp"
            quantafilenamep="quantabigRp.tmp"
            Hindex=2
        end if
        recindex=abs(orbnow-orbref)+1

!----------------open a binary file-------------
        ! rowindex
        reclength=subMp*4+1
        inquire(file=trim(Hrowfilenamep),exist=alive)
        if(alive) then
            open(unit=112,file=trim(Hrowfilenamep),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(Hrowfilenamep),"doesn't exist!"
            stop
        end if
        ! sigmaR index =norbs is 1; norbs-1 is 2
        read(112,rec=recindex) Hbigrowindexp(:,Hindex)
        close(112)

        nonzero=Hbigrowindexp(4*dim1+1,Hindex)-1   ! nonzero element numbers in Hbig
        
        reclength=2*Hbigdimp
        inquire(file=trim(Hfilenamep),exist=alive)
        if(alive) then
            open(unit=102,file=trim(Hfilenamep),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(Hfilenamep),"doesn't exist!"
            stop
        end if
! norbs is 1; norbs-1 is 2
        read(102,rec=recindex) Hbigp(1:nonzero,Hindex)
        close(102)

        ! colindex
        reclength=Hbigdimp
        inquire(file=trim(Hcolfilenamep),exist=alive)
        if(alive) then
            open(unit=110,file=trim(Hcolfilenamep),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(Hcolfilenamep),"doesn't exist!"
            stop
        end if
        ! sigmaR index =norbs is 1; norbs-1 is 2
        read(110,rec=recindex) Hbigcolindexp(1:nonzero,Hindex)
        close(110)
!------------------------------read the quantabigR/L(4*subM,2)---
! every process need to read quantabigR/L

        reclength=4*subMp*2
        ! quantabigR/L is  integer(kind=4) so without *2
        inquire(file=trim(quantafilenamep),exist=alive)
        if(alive) then
            open(unit=108,file=trim(quantafilenamep),access="Direct",form="unformatted",recl=reclength,status="old")
        else
            write(*,*) trim(quantafilenamep),"doesn't exist!"
            stop
        end if
        if(domain=='L') then
            read(108,rec=recindex) quantabigLp
        else if(domain=='R') then
            read(108,rec=recindex) quantabigRp
        end if
        close(108)
    end if

    ! broadcast the quantabigR/L to other process
    if(domain=='L') then
        call MPI_BCAST(quantabigLp,4*subMp*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
    else if(domain=='R') then
        call MPI_BCAST(quantabigRp,4*subMp*2,MPI_integer4,0,MPI_COMM_WORLD,ierr)
    end if

    call MPI_FILE_CLOSE(thefile,ierr)
    end if

return

end Subroutine Enviro_Big
