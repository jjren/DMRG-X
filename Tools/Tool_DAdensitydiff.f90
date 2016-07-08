program Tool_DAdensitydiff
    implicit none

    integer :: nbonds,nsites,nstates,ncutoff,nDcount,nAcount
    integer :: istate,ibond,isite,idummy
    real(kind=8) :: updensity,downdensity,Ddensity,Adensity,diffdensity


    write(*,*) "nbonds"
    read(*,*) nbonds
    write(*,*) "nsites"
    read(*,*) nsites
    write(*,*) "nstates"
    read(*,*) nstates
    write(*,*) "ncutoff"
    read(*,*) ncutoff

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
        nDcount=0
        nAcount=0
        do isite=1,nsites,1
            read(10,*) idummy,idummy,updensity,downdensity
            if(isite>ncutoff .and. (nsites-isite)>=ncutoff) then
                if(mod(isite,2)==1) then
                    Ddensity=Ddensity+updensity+downdensity
                    nDcount=nDcount+1
                else
                    Adensity=Adensity+updensity+downdensity
                    nAcount=nAcount+1
                end if
            end if
        end do
        diffdensity=Ddensity/nDcount-Adensity/nAcount
        write(11,*) diffdensity 
    end do

    close(10)
    close(11)
end program Tool_DAdensitydiff
