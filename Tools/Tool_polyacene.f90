program Tool_polyacene

    implicit none
    
    
    real(kind=8),parameter :: bonda=1.397D0,bondb=1.397D0
    integer :: nrings,atomindex(6),nbonds,norbs
    real(kind=8) :: t,hubbardU,siteenergy,nuclQ
    integer :: i
    character(len=1) :: symbol
    real(kind=8),allocatable :: coord(:,:)
    
    write(*,*) "input nrings:"
    read(*,*) nrings
    write(*,*) "input hopping integral t:"
    read(*,*) t
    write(*,*) "input hubbardU :"
    read(*,*) hubbardU
    write(*,*) "input siteenergy:"
    read(*,*) siteenergy
    
    nbonds=5*nrings+1
    norbs=4*nrings+2
    
    write(*,*) "norbs=",norbs
    write(*,*) "nbonds=",nbonds
    
    ! creat the bondlink information
    do i=1,6,1
        atomindex(i)=i
    end do
    open(unit=10,file="integral.inp",status="replace")
    
    do i=1,nrings,1
        write(10,*) atomindex(1),atomindex(2),t
        write(10,*) atomindex(1),atomindex(3),t
        write(10,*) atomindex(2),atomindex(4),t
        write(10,*) atomindex(3),atomindex(5),t
        write(10,*) atomindex(4),atomindex(6),t
        atomindex=atomindex+4
    end do

    atomindex=atomindex-4
    write(10,*) atomindex(5),atomindex(6),t

    ! site energy
    do i=1,norbs,1
        write(10,*) siteenergy
    end do

    ! hubbardU
    do i=1,norbs,1
        write(10,*) hubbardU
    end do

    close(10)

    allocate(coord(3,norbs))
    nuclQ=1.0D0
    symbol='C'

    open(unit=12,file="coord.xyz",status="replace")
    write(12,*) norbs
    write(12,*)
    do i=1,norbs,1
        coord(3,i)=0.0D0
        if(Mod(i,4)==1) then
            coord(1,i)=2.0D0*bonda*sqrt(3.0D0)/2.0D0*dble((i-1)/4)
            coord(2,i)=0.0D0
        else if(Mod(i,4)==2) then
            coord(1,i)=2.0D0*bonda*sqrt(3.0D0)/2.0D0*dble((i-2)/4)
            coord(2,i)=-bondb
        else if(Mod(i,4)==3) then
            coord(1,i)=bonda*sqrt(3.0D0)/2.0D0*dble((i-1)/2)
            coord(2,i)=0.5D0*bonda
        else
            coord(1,i)=bonda*sqrt(3.0D0)/2.0D0*dble(i/2-1)
            coord(2,i)=-bondb-0.5D0*bonda
        end if
    
        write(12,'(1A,4F16.12)') symbol,coord(1:3,i),nuclQ
    end do
    close(12)

    deallocate(coord)
end 

