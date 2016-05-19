program Tool_DAdistortion
    implicit none

    integer :: nbonds,ibond,idummy,ncutoff,ncount
    real(kind=8) :: localdelta,totaldelta

    write(*,*) "nbonds"
    read(*,*) nbonds
    write(*,*) "ncutoff"
    read(*,*) ncutoff

    open(unit=10,file="peierlsdelta.out",status="old")
    read(10,*)
    totaldelta=0.0D0
    ncount=0
    do ibond=1,nbonds,1
        read(10,*) idummy,idummy,localdelta
        if(ibond>ncutoff .and. (nbonds-ibond)>=ncutoff) then
            if(mod(ibond,2)==1) then
                totaldelta=totaldelta+(-1.0D0)*localdelta
            else
                totaldelta=totaldelta+localdelta
            end if
            ncount=ncount+1
        endif
    end do
    totaldelta=totaldelta/DBLE(ncount)

    close(10)

    open(unit=11,file="avedelta.out",status="replace")
    write(11,*) totaldelta
    close(11)
end 

