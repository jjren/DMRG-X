module Peierls_mod
    USE MPI
    use communicate
    use kinds_mod
    use variables

    implicit none
    integer :: peierlsstate,npeierlsloops
    real(kind=r8) :: hopalpha,springK,peierlsconvergethresh
    real(kind=r8),allocatable :: peierlsT(:,:),peierlsD(:,:),peierlsdelta(:,:)
contains
!================================================
!================================================

subroutine Allocate_Peierls
    implicit none
    allocate(peierlsT(norbs,norbs))
    peierlsT=0.0D0
    allocate(peierlsdelta(norbs,norbs))
    peierlsdelta=0.0D0
    if(ifpeierlsD==1) then 
        allocate(peierlsD(norbs,norbs))
        peierlsD=0.0D0
    end if

    return
end subroutine Allocate_Peierls

!================================================
!================================================
 
subroutine Peierls_init(method)
    implicit none
    character(len=*),intent(in) :: method
    logical :: alive
    integer :: iorb,jorb,iterm,nterms

    inquire(file="peierlsdelta.out",exist=alive)
    if(alive) then
        if(myid==0) then
            call Allocate_Peierls
            open(unit=775,file="peierlsdelta.out",status="old")
            read(775,*) nterms
            do iterm=1,nterms,1
                read(775,*) iorb,jorb,peierlsdelta(iorb,jorb)
                peierlsdelta(jorb,iorb)=peierlsdelta(iorb,jorb)
            end do
            close(775)
        end if
        call UpdateIntegral(method)
        
        if(myid==0) call Deallocate_Peierls
    end if

    return
end Subroutine Peierls_init

!================================================
!================================================

subroutine Peierls_driver(ifconverge)
    implicit none
    logical,intent(out) :: ifconverge
    
    if(logic_C2/=0) then
        nstate=nstate*2
    end if

    call Allocate_Peierls
    call PeierlsTD_Read
    call Cal_Peierlsdelta
    
    call PeierlsdeltaCompare(ifconverge)
    
    call Peierls_Write
    call Deallocate_Peierls
    
    ! recover the input nstate
    if(logic_C2/=0) then
        nstate=nstate/2
    end if

    return
end subroutine Peierls_driver

!================================================
!================================================

subroutine PeierlsTD_Read
    implicit none
    !local
    integer :: i,j,istate,k
    integer :: iorb,jorb
    real(kind=r8) :: boup,bodown

    open(unit=777,file="bondord.out",status="old")
    do istate=1,nstate,1
        read(777,*) k
        do i=1,norbs,1
        do j=i,norbs,1
            if(bondlink(i,j)==1) then
                read(777,*) iorb,jorb,boup,bodown
                if(i/=iorb .or. j/=jorb) then
                    write(*,*) "Peierls_read failed!",iorb,jorb,i,j 
                    stop
                end if
                peierlsT(iorb,jorb)=boup+bodown
                peierlsT(jorb,iorb)=peierlsT(iorb,jorb)
            end if
        end do
        end do
        if(k==peierlsstate) exit
    end do
    if(k/=peierlsstate) then
        write(*,*) "nopeierlsstate"
        stop
    end if
    close(777)
    
    if(ifpeierlsD==1) then 
        open(unit=778,file="peierlsD.out",status="old")
        read(778,*) iorb,jorb,peierlsD(iorb,jorb)
        peierlsD(jorb,iorb)=peierlsD(iorb,jorb)
        close(778)
    end if

    return

end subroutine PeierlsTD_Read

!================================================
!================================================

Subroutine Cal_Peierlsdelta
    implicit none
    integer :: iorb,jorb
    real(kind=r8) :: eta
    
    eta=0.0D0
    do iorb=1,norbs,1
    do jorb=iorb+1,norbs,1
        if(bondlink(iorb,jorb)==1) then
            peierlsdelta(iorb,jorb)=-1.0D0*hopalpha*peierlsT(iorb,jorb)*2.0D0
            if(IfpeierlsD==1) then
                peierlsdelta(iorb,jorb)=peierlsdelta(iorb,jorb)+pppw(iorb,jorb)*peierlsD(iorb,jorb)
            end if
            peierlsdelta(jorb,iorb)=peierlsdelta(iorb,jorb)
            eta=-peierlsdelta(iorb,jorb)+eta
        end if
    end do
    end do
    eta=eta/DBLE(nbonds)
    do iorb=1,norbs,1
    do jorb=iorb+1,norbs,1
        peierlsdelta(iorb,jorb)=(peierlsdelta(iorb,jorb)+eta)/springK
        peierlsdelta(jorb,iorb)=peierlsdelta(iorb,jorb)
    end do
    end do

    return
end Subroutine Cal_Peierlsdelta

!================================================
!================================================

Subroutine Peierls_Write
    implicit none
    integer :: iorb,jorb

    open(unit=776,file="peierlsdelta.out",status="replace")
    write(*,*) "peierlsdelta="
    write(776,*) nbonds
    do iorb=1,norbs,1
    do jorb=iorb+1,norbs,1
        if(bondlink(iorb,jorb)==1) then
            write(776,*) iorb,jorb,peierlsdelta(iorb,jorb)
            write(*,*) iorb,jorb,peierlsdelta(iorb,jorb)
        end if
    end do
    end do
    close(776)
    
    return

end Subroutine Peierls_Write

!================================================
!================================================

Subroutine PeierlsdeltaCompare(ifconverge)
    implicit none
    
    logical,intent(out) :: ifconverge
    ! local
    real(kind=r8),allocatable :: deltaold(:,:)
    integer :: nterms,iorb,jorb,iterm
    real(kind=r8) :: diff
    logical :: alive

    allocate(deltaold(norbs,norbs))
    inquire(file="peierlsdelta.out",exist=alive)
    if(alive) then
        ifconverge=.true.
        open(unit=779,file="peierlsdelta.out",status="old")
        read(779,*) nterms
        do iterm=1,nterms,1
            read(779,*) iorb,jorb,deltaold(iorb,jorb)
            diff=deltaold(iorb,jorb)-peierlsdelta(iorb,jorb)
            if(abs(diff) > peierlsconvergethresh) then
                ifconverge=.false.
                exit
            end if
        end do
        deallocate(deltaold)
        close(779)
    else
        ifconverge=.false.
    end if

    return
end subroutine PeierlsdeltaCompare
!================================================
!================================================

Subroutine UpdateIntegral(method)
    use  PPP_term_mod
    implicit none
    character(len=*),intent(in) :: method

    ! local
    integer :: iorb,jorb
    ! MPI
    integer :: ierr

    if(myid==0) then
        call ReadIntegral
        call PPP_term
        do iorb=1,norbs,1
        do jorb=iorb+1,norbs,1
            if(bondlink(iorb,jorb)==1) then
                t(iorb,jorb)=t(iorb,jorb)+hopalpha*peierlsdelta(iorb,jorb)
                t(jorb,iorb)=t(iorb,jorb)
                if(IfpeierlsD==1) then
                    pppV(iorb,jorb)=pppV(iorb,jorb)-pppw(iorb,jorb)*peierlsdelta(iorb,jorb)
                    pppV(jorb,iorb)=pppV(iorb,jorb)
                end if
            end if
        end do
        end do
    end if
    
    if(method=="DMRG") then
        call MPI_BCAST(t,norbs*norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        if(IfpeierlsD==1) then
            call MPI_BCAST(pppV,norbs*norbs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        end if
    end if

    if(myid==0) then
        write(*,*) "new t="
        write(*,*) t
        if(IfpeierlsD==1) then
            write(*,*) "new pppV="
            write(*,*) pppV
        end if
    end if
    return
end Subroutine UpdateIntegral

!================================================
!================================================

subroutine Deallocate_Peierls
    implicit none
    deallocate(peierlsT)
    deallocate(peierlsdelta)
    if(ifpeierlsD==1) then 
        deallocate(peierlsD)
    end if

    return
end subroutine Deallocate_Peierls

!================================================
!================================================

end module Peierls_mod
