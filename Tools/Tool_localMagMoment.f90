program Tool_localMagMoment
    implicit none

    integer :: nsites,nstates,istate,isite,idummy,ncutoff,Dncount,Ancount
    real(kind=8) :: localmag,Dtotallocalmag,Atotallocalmag
    
    write(*,*) "nsites"
    read(*,*) nsites
    write(*,*) "nstates"
    read(*,*) nstates
    write(*,*) "ncutoff"
    read(*,*) ncutoff

    open(unit=10,file="Llocalmagmoment.out",status="old")
    open(unit=11,file="Rlocalmagmoment.out",status="old")
    open(unit=12,file="avelocalmag.out",status="replace")
    
    do istate=1,nstates,1
        read(10,*)
        read(11,*)
        Dtotallocalmag=0.0D0
        Atotallocalmag=0.0D0
        Ancount=0
        Dncount=1

        do isite=1,nsites/2,1
            read(10,*) idummy,localmag
            if(idummy>ncutoff .or. nsites-idummy>=ncutoff) then
                if(mod(idummy,2)==1) then
                    Dtotallocalmag=Dtotallocalmag+localmag
                    Dncount=Dncount+1
                else
                    Atotallocalmag=Atotallocalmag+localmag
                    Ancount=Ancount+1
                end if
            end if

            read(11,*) idummy,localmag
            if(idummy>ncutoff .or. nsites-idummy>=ncutoff) then
                if(mod(idummy,2)==1) then
                    Dtotallocalmag=Dtotallocalmag+localmag
                    Dncount=Dncount+1
                else
                    Atotallocalmag=Atotallocalmag+localmag
                    Ancount=Ancount+1
                end if
            end if
        end do
        Dtotallocalmag=Dtotallocalmag/DBLE(Dncount)
        Atotallocalmag=Atotallocalmag/DBLE(Ancount)
        write(12,*) Dtotallocalmag
        write(12,*) Atotallocalmag
        write(12,*) Dtotallocalmag+Atotallocalmag
    end do

    close(10)
    close(11)
    close(12)

end program Tool_localMagMoment

