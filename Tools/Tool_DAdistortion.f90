program Tool_DAdistortion
    implicit none

    integer :: nbonds,ibond,idummy,ncutoff,ncountA,ncountD
    real(kind=8) :: localdelta,totaldelta,totaldeltaA,totaldeltaD

    write(*,*) "nbonds"
    read(*,*) nbonds
    write(*,*) "ncutoff"
    read(*,*) ncutoff

    open(unit=10,file="peierlsdelta.out",status="old")
    read(10,*)
    totaldeltaD=0.0D0
    totaldeltaA=0.0D0
    ncountD=0
    ncountA=0
    do ibond=1,nbonds,1
        read(10,*) idummy,idummy,localdelta
        if(ibond>ncutoff .and. (nbonds-ibond)>=ncutoff) then
            if(mod(ibond,2)==1) then
                totaldeltaD=totaldeltaD+localdelta
                ncountD=ncountD+1
            else
                totaldeltaA=totaldeltaA+localdelta
                ncountA=ncountA+1
            end if
        endif
    end do
    totaldelta=totaldeltaD/DBLE(ncountD)-totaldeltaA/DBLE(ncountA)

    close(10)

    open(unit=11,file="avedelta.out",status="replace")
    write(11,*) totaldelta
    close(11)
end 

