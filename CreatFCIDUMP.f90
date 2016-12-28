subroutine CreatFCIDUMP
    use variables

    implicit none
    integer,parameter :: isym=1,Ms2=0,izero=0
    integer,allocatable :: orbirrep(:)
    real(kind=8) :: onsite,coreenergy
    integer :: i,j,p

    allocate(orbirrep(norbs))
    orbirrep=1

    open(unit=14,file="FCIDUMP",status="replace")
    
    !! the header part
    write(14,'(1X,''&FCI NORB='',I3,'',NELEC='',I2,'',MS2='',I2,'','')') norbs,nelecs,Ms2
    write(14,'(2X,''ORBSYM='',30(I1,'',''),6(/1X,40(I1,'','')))') (orbirrep(i),i=1,norbs)
    write(14,'(2X,''ISYM='',I1)') isym
    write(14,'(1X,''&END'')')
    
    ! hubbardU and pppV
    do i=1,norbs,1
    do j=1,i,1
        if(i==j) then
            write(14,*) hubbardU(i),i,i,i,i
        else
            write(14,*) pppV(i,j),i,i,j,j
        end if
    end do
    end do

    ! transfer integral and site energy
    do i=1,norbs,1
    do j=1,i,1
        if(bondlink(j,i)==1) then
            write(14,*) t(i,j),i,j,izero,izero
        else if(j==i) then
            onsite=t(i,i)
            do p=1,norbs,1
                if(p/=i) then
                onsite=onsite-nuclQ(p)*pppV(p,i)
                end if
            end do
            write(14,*) onsite,i,j,izero,izero
        !else
        !    write(14,'(E28.20,4I) zero,i,j,izero,izero
        end if
    end do
    end do

    coreenergy=0.0D0
    do i=1,norbs,1
    do j=1,i-1,1
        coreenergy=coreenergy+pppV(j,i)*nuclQ(i)*nuclQ(j)
    end do
    end do

    write(14,*) coreenergy,izero,izero,izero,izero
    
    ! fullci.py format
    open(unit=61,file="MO-active-dipole.out",status="replace")
    write(61,*) norbs
    do i=1,norbs,1
    do j=1,i,1
        if(j==i) then
            write(61,*) coord(1:3,i)
        else
            write(61,*) 0.0D0,0.0D0,0.0D0
        end if
    end do
    end do
    close(61)


    deallocate(orbirrep)
        
end subroutine CreatFCIDUMP
