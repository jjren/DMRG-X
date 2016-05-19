program Tool_DAbondorder
    implicit none

    integer :: nbonds,ncutoff,ncount
    integer :: ibond,idummy
    real(kind=8) ::  DBO,ABO,diffBO,upBO,downBO
    
    write(*,*) "nbonds"
    read(*,*) nbonds
    write(*,*) "ncutoff"
    read(*,*) ncutoff

    open(unit=10,file="bondord.out",status="old")
    
    read(10,*) 
    DBO=0.0D0
    ABO=0.0D0
    ncount=0
    do ibond=1,nbonds,1
        read(10,*) idummy,idummy,upBO,downBO
        if(ibond>ncutoff .and. (nbonds-ibond)>=ncutoff) then
            if(mod(ibond,2)==1) then
                DBO=DBO+upBO+downBO
            else
                ABO=ABO+upBO+downBO
            end if
            ncount=ncount+1
        end if
    end do
    diffBO=DBO/(ncount+1)*2-ABO/ncount*2
    
    open(unit=11,file="DABOdiff.out",status="replace")
    write(11,*) diffBO 
    close(11)
end program Tool_DAbondorder
