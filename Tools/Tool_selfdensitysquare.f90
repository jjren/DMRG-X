program Tool_selfdensitysquare
    implicit none
    
    integer ::  norbs,nstates,iorb,istate,jstate,index1
    real(kind=8),allocatable :: localmag(:),transDM0(:,:,:,:,:)

    open(unit=398,file="AO-transOpdm.out",form="unformatted",status="old")
    read(398) norbs,nstates
    allocate(transDM0(norbs,norbs,2,nstates,nstates))
    do istate=1,nstates,1
    do jstate=istate,nstates,1
        read(398) transDM0(:,:,:,jstate,istate)
    end do
    end do
    close(398)

    allocate(localmag(norbs))

    open(unit=10,file="Llocalmagmoment.out",status="old")
    open(unit=11,file="Rlocalmagmoment.out",status="old")

    do istate=1,nstates,1
        read(10,*)
        read(11,*)
        do iorb=1,norbs/2,1
            read(10,*) index1,localmag(index1)
            read(11,*) index1,localmag(index1)
        end do
    end do
    close(10)
    close(11)
    
    open(unit=20,file="localdensitysquare.out",status="replace")
    do iorb=1,norbs,1
        write(20,*) iorb,2.0D0*(transDM0(iorb,iorb,1,1,1)+transDM0(iorb,iorb,2,1,1))-localmag(iorb)- &
            (transDM0(iorb,iorb,1,1,1)+transDM0(iorb,iorb,2,1,1))**2
    end do
    close(20)

    deallocate(localmag,transDM0)
        
end program Tool_selfdensitysquare
