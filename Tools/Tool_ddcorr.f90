program  Tool_ddcorr
    implicit none

    integer :: norbs,nstates
    integer :: istate,jstate,iorb,orbindex,orbadd
    real(kind=8) :: ddcorr
    real(kind=8),allocatable :: transDM0(:,:,:,:,:)

    open(unit=398,file="AO-transOpdm.out",form="unformatted",status="old")
    read(398) norbs,nstates
    allocate(transDM0(norbs,norbs,2,nstates,nstates))
    do istate=1,nstates,1
    do jstate=istate,nstates,1
        read(398) transDM0(:,:,:,jstate,istate)
    end do
    end do
    close(398)

    open(unit=11,file="ddcorr.out",status="replace")
    
    open(unit=10,file="rddcorr.out",status="old")
    orbadd=norbs/2+1
    do iorb=1,norbs/2,1
        read(10,*) orbindex,ddcorr
        write(11,*) orbindex,ddcorr-(transDM0(orbadd,orbadd,1,1,1)+transDM0(orbadd,orbadd,2,1,1)) * &
                                (transDM0(orbindex,orbindex,1,1,1)+transDM0(orbindex,orbindex,2,1,1))
    end do
    close(10)
    
    open(unit=10,file="lddcorr.out",status="old")
    orbadd=norbs/2
    do iorb=1,norbs/2,1
        read(10,*) orbindex,ddcorr
        write(11,*) orbindex,ddcorr-(transDM0(orbadd,orbadd,1,1,1)+transDM0(orbadd,orbadd,2,1,1)) * &
                                (transDM0(orbindex,orbindex,1,1,1)+transDM0(orbindex,orbindex,2,1,1))
    end do
    close(10)
    close(11)
    deallocate(transDM0)
end program
