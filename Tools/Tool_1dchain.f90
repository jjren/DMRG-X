Program Tool_1dchain
! creat the 1d polyene chain
! the default angle is 120
    implicit none
    real(kind=8) :: bondlength,delta,deltalength,hubbardUD,hubbardUA,t
    integer(kind=4) :: nsite
    real(kind=8),allocatable :: coord(:,:)
    real(kind=8) :: t1,t2,singlelength,doublelength,ts,td,dist1,dist3,siteenergyD,siteenergyA
    integer :: i

    !write(*,*) "input angle:"
    !read(*,*) angle
    write(*,*) "input bondlength:"
    read(*,*) bondlength
    write(*,*) "input deltalength:"
    read(*,*) deltalength
    write(*,*) "input delta:"
    read(*,*) delta
    write(*,*) "nsite:"
    read(*,*) nsite
    write(*,*) "D siteenergy"
    read(*,*) siteenergyD
    write(*,*) "A siteenergy"
    read(*,*) siteenergyA
    write(*,*) "hopping integral:"
    read(*,*) t
    write(*,*) "hubbardU D"
    read(*,*) hubbardUD
    write(*,*) "hubbardU A"
    read(*,*) hubbardUA


    allocate(coord(3,nsite))

    singlelength=bondlength+deltalength
    doublelength=bondlength-deltalength

    dist1=doublelength+singlelength
    dist3=doublelength
    ! coordinate
    coord(1:3,1)=0.0D0
    do i=2,nsite,1
        if(mod(i,2)==1) then
            coord(1,i)=dist1*DBLE((i-1)/2)
            coord(2,i)=0.0D0
            coord(3,i)=0.0D0
        else
            coord(1,i)=coord(1,i-1)+dist3
            coord(2,i)=0.0D0
            coord(3,i)=0.0D0
        end if
    end do

    open(unit=11,file="coord.xyz",status="replace")
    write(11,*) nsite
    write(11,*) 
    do i=1,nsite,1
        write(11,'(1A,4F10.5)') 'C',coord(1:3,i),1.0D0
    end do
    close(11)

    open(unit=13,file="integral.inp",status="replace")
    ts=t*(1.0D0-delta)
    td=t*(1.0D0+delta)
    do i=1,nsite-1,1
        if(mod(i,2)==1) then
            write(13,*) i,i+1,td
        else
            write(13,*) i,i+1,ts
        end if
    end do
    do i=1,nsite,1
        if(mod(i,2)==1) then
            write(13,*) siteenergyD
        else
            write(13,*) siteenergyA
        end if
    end do
    do i=1,nsite,1
        if(mod(i,2)==1) then
            write(13,*) hubbardUD
        else
            write(13,*) hubbardUA
        end if
    end do

    close(13)
    deallocate(coord)
end Program Tool_1dchain
