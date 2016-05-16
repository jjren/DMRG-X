program Tool_DAdensitydiff
    implicit none

    integer :: nbonds,nsites,nstates
    integer :: istate,ibond,isite,idummy
    real(kind=8) :: updensity,downdensity,Ddensity,Adensity,diffdensity


    write(*,*) "nbonds"
    read(*,*) nbonds
    write(*,*) "nsites"
    read(*,*) nsites
    write(*,*) "nstates"
    read(*,*) nstates

    open(unit=10,file="bondord.out",status="old")
    open(unit=11,file="DAdensityDiff.out",status="replace")

    do istate=1,nstates,1
        read(10,*) 
        do ibond=1,nbonds,1
            read(10,*) 
        end do
    end do
    

    do istate=1,nstates,1
        read(10,*) 
        Ddensity=0.0D0
        Adensity=0.0D0
        do isite=1,nsites,1
            read(10,*) idummy,idummy,updensity,downdensity
            if(mod(isite,2)==1) then
                Ddensity=Ddensity+updensity+downdensity
            else
                Adensity=Adensity+updensity+downdensity
            end if
        end do
        diffdensity=(Ddensity-Adensity)/nsites*2
        write(11,*) diffdensity 
    end do

    close(10)
    close(11)
end program Tool_DAdensitydiff
