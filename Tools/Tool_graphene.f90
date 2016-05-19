program Tool_polyacene

    implicit none
    
    
    real(kind=8),parameter :: bonda=1.397D0,bondb=1.397D0
    integer :: nrings,atomindex(6),norbs,nwidth,nchainorbs
    real(kind=8) :: t,hubbardU,siteenergy,nuclQ,distance
    integer :: i,j,iorb
    character(len=1) :: symbol
    real(kind=8),allocatable :: coord(:,:,:),coordorder(:,:)
    
    write(*,*) "input nrings:"
    read(*,*) nrings
    write(*,*) "input nwidth:"
    read(*,*) nwidth
    write(*,*) "input hopping integral t:"
    read(*,*) t
    write(*,*) "input hubbardU :"
    read(*,*) hubbardU
    write(*,*) "input siteenergy:"
    read(*,*) siteenergy
    
    nchainorbs=4*nrings+2
    norbs=nchainorbs*nwidth
    
    
    nuclQ=1.0D0
    symbol='C'
    
    open(unit=12,file="coord.xyz",status="replace")
    write(12,*) norbs
    write(12,*)
    
    allocate(coord(3,nchainorbs,nwidth))
    allocate(coordorder(3,norbs))
    coord=0.0D0

    do i=1,nchainorbs,1
        if(Mod(i,4)==1) then
            coord(1,i,1)=2.0D0*bonda*sqrt(3.0D0)/2.0D0*dble((i-1)/4)
            coord(2,i,1)=0.0D0
        else if(Mod(i,4)==2) then
            coord(1,i,1)=2.0D0*bonda*sqrt(3.0D0)/2.0D0*dble((i-2)/4)
            coord(2,i,1)=-bondb
        else if(Mod(i,4)==3) then
            coord(1,i,1)=bonda*sqrt(3.0D0)/2.0D0*dble((i-1)/2)
            coord(2,i,1)=0.5D0*bonda
        else
            coord(1,i,1)=bonda*sqrt(3.0D0)/2.0D0*dble(i/2-1)
            coord(2,i,1)=-bondb-0.5D0*bonda
        end if
    end do

    do i=2,nwidth,1
        coord(1,:,i)=coord(1,:,1)
        coord(2,:,i)=coord(2,:,1)-(i-1)*2*bondb-(i-1)*bonda
    end do

    iorb=0
    do i=1,nchainorbs/2,1
        do j=1,nwidth,1
            write(12,'(1A,4F20.12)') symbol,coord(1:3,i*2-1,j),nuclQ
            write(12,'(1A,4F20.12)') symbol,coord(1:3,i*2,j),nuclQ
            iorb=iorb+2
            coordorder(1:3,iorb-1:iorb)=coord(1:3,i*2-1:i*2,j)
        end do
    end do

    if(iorb/=norbs) then
        write(*,*) iorb,norbs
        stop
    end if
    close(12)
    
    ! creat the bondlink information
    open(unit=10,file="integral.inp",status="replace")
    
    do i=1,norbs,1
    do j=i,norbs,1
        distance=(coordorder(1,i)-coordorder(1,j))**2+(coordorder(2,i)-coordorder(2,j))**2+(coordorder(3,i)-coordorder(3,j))**2
        if(ABS(distance-bonda**2)<=1.0D-4) then
            write(10,*) i,j,t
        end if
    end do
    end do

    ! site energy
    do i=1,norbs,1
        write(10,*) siteenergy
    end do

    ! hubbardU
    do i=1,norbs,1
        write(10,*) hubbardU
    end do

    write(*,*) "norbs=",norbs
    
    close(10)

    deallocate(coord,coordorder)
end 

