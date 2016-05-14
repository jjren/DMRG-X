subroutine Polarization(istate)
    ! point charge model
    ! only useful in 1-d system
    
    implicit none
    real(kind=r8) :: ionicpolar, elecpolar

    allocate(delta(norbs,norbs))

    inquire(file="peierlsdelta.out",exist=alive)
    if(alive) then
        open(unit=779,file="peierlsdelta.out",status="old")
        read(779,*) nterms
        do iterm=1,nterms,1
            read(779,*) iorb,jorb,deltaold(iorb,jorb)
        end do
        close(779)
    else
        write(*,*) "peierlsdelta.out doesn't exist"
        delta=0.0D0
    end if
    
    ionicpolar=0.0D0
    do i=1,norbs,1
        do icoord=1,3,1
            ionicpolar(icoord) = ionicpolar(icoord)+(coord(icoord,i)+delta(i.i+1)-coord(icoord,1))*nuclQ(i)
        end do
    end do

    elecpolar=0.0D0
    allocate(transDM0(norbs,norbs,2,nstate,nstate))
    inquire(file="AO-transOpdm.out",exist=alive)
    if(alive) then
        open(unit=398,file="AO-transOpdm.out",form="unformatted",status="old")
        read(398) norbsdummy,nstatedummy
        if(norbsdummy/=norbs) then
            write(*,*) "norbsdummy/=norbs",norbsdummy,norbs
            stop
        end if
        if(nstatedummy/=nstate) then
            write(*,*) "nstatedummy/=nstate",nstatedummy,nstate
            stop
        end if
        do istate=1,nstate,1
        do jstate=istate,nstate,1
            read(398) transDM0(:,:,:,jstate,istate)
        end do
        end do
        close(398)


    else
        write(*,*) "no AO-transOpdm.out exist!"
    end if
    
    
    return

end subroutine Polarization

