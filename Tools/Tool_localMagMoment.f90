program Tool_localMagMoment
    implicit none

    integer :: nsites,nstates,istate,isite,idummy
    real(kind=8) :: localmag,totallocalmag
    
    write(*,*) "nsites"
    read(*,*) nsites
    write(*,*) "nstates"
    read(*,*) nstates

    open(unit=10,file="Llocalmagmoment.out",status="old")
    open(unit=11,file="Rlocalmagmoment.out",status="old")
    open(unit=12,file="avelocalmag.out",status="replace")
    
    do istate=1,nstates,1
        read(10,*)
        read(11,*)
        totallocalmag=0.0D0
        do isite=1,nsites/2,1
            read(10,*) idummy,localmag
            totallocalmag=totallocalmag+localmag
            read(11,*) idummy,localmag
            totallocalmag=totallocalmag+localmag
        end do
        totallocalmag=totallocalmag/DBLE(nsites)
        write(12,*) totallocalmag
    end do

            

end program Tool_localMagMoment

