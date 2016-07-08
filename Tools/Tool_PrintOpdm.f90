program Tool_PrintOpdm
    implicit none
    integer :: norbs,nstates,nbonds
    real(kind=8),allocatable :: transDM0(:,:,:,:,:)
    integer,allocatable :: bondindex(:,:)
    integer :: istate,jstate,ibond

    open(unit=398,file="AO-transOpdm.out",form="unformatted",status="old")
    read(398) norbs,nstates
    allocate(transDM0(norbs,norbs,2,nstates,nstates))
    do istate=1,nstates,1
    do jstate=istate,nstates,1
        read(398) transDM0(:,:,:,jstate,istate)
    end do
    end do
    close(398)

    open(unit=10,file="printopdm.inp",status="old")
    open(unit=11,file="printopdm.out",status="replace")
    read(10,*) nbonds
    allocate(bondindex(2,nbonds))
    do ibond=1,nbonds,1
        read(10,*) bondindex(1:2,ibond)
    end do

    do istate=1,nstates,1
        write(11,*) "#######",istate
        do ibond=1,nbonds,1
            write(11,*) ibond,transDM0(bondindex(1,ibond),bondindex(2,ibond),1:2,istate,istate)
        end do
    end do
    
    deallocate(bondindex)
    deallocate(transDM0)
end program Tool_PrintOpdm

