program Tool_DApolarization
    implicit none
    
    real(kind=8),allocatable :: delta(:,:),density(:,:),ionic(:,:),elec(:,:)
    integer :: nalpha,nhubbard
    integer :: ialpha,ihubbard,idummy
    real :: bondlength

    write(*,*) "nalpha"
    read(*,*) nalpha
    write(*,*) "nhubbard"
    read(*,*) nhubbard
    write(*,*) "bondlength"
    read(*,*) bondlength

    allocate(delta(nhubbard,nalpha))
    allocate(density(nhubbard,nalpha))
    allocate(ionic(nhubbard,nalpha))
    allocate(elec(nhubbard,nalpha))

    open(unit=10,file="DApolarization.out",status="old")
    do ialpha=1,nalpha,1
        read(10,*) 
        read(10,*) 
        do ihubbard=1,nhubbard,1
            read(10,*) idummy,density(ihubbard,ialpha),delta(ihubbard,ialpha)
        end do
        read(10,*) 
        read(10,*) 
    end do
    close(10)

    do ialpha=1,nalpha,1
    do ihubbard=1,nhubbard,1
        ionic(ihubbard,ialpha)=(bondlength-delta(ihubbard,ialpha))*1.0D0 - &
            ((bondlength-delta(ihubbard,1))*1.0D0)
    end do
    end do
        
    do ialpha=1,nalpha,1
    do ihubbard=1,nhubbard,1
        elec(ihubbard,ialpha)=-(bondlength-delta(ihubbard,ialpha))*(2.0D0-density(ihubbard,ialpha))/2.0D0 + &
            ((bondlength-delta(ihubbard,1))*(2.0D0-density(ihubbard,1))/2.0D0)
    end do
    end do

    open(unit=11,file="polarization.out",status="replace")
    do ialpha=1,nalpha,1
        write(11,*) (ialpha-1)*0.1
        write(11,*) "######################"
        do ihubbard=1,nhubbard,1
            write(11,'(1I,3F15.10)') ihubbard-1, ionic(ihubbard,ialpha),elec(ihubbard,ialpha),ionic(ihubbard,ialpha)+elec(ihubbard,ialpha)
        end do
        write(11,*)
        write(11,*)
    end do
    close(11)

    deallocate(delta)
    deallocate(density)
    deallocate(ionic)
    deallocate(elec)


end program Tool_DApolarization
