Module PPP_term_mod
!this module can include many potential terms like ohno potential
    use variables, only : hubbardU,coord,pppV,norbs,logic_PPP,PPPpot,&
        logic_Peierls,IfpeierlsD,pppw
    use communicate
    use kinds_mod
    
    implicit none

    integer(kind=i4),allocatable :: pppVlink(:,:)
    
    contains

!============================================================
!============================================================

subroutine PPP_term
   
    implicit none
    if(logic_PPP==1) then
        if(PPPpot=="OK") then
            call Ohno_Potential
        else if(PPPpot=="EHP") then
            call Extented_Hubbard_Peierls_Potential
        end if
    else if(logic_PPP==0) then
        ! hubbard model
        pppV=0.0D0
        pppVlink=0
    end if
            
return
end subroutine PPP_term

!============================================================
!============================================================

subroutine Extented_Hubbard_Peierls_Potential
    ! only the nearest neighbor potential
    use variables,only : bondlink
    implicit none
    
    real(kind=r8) :: EHPpot
    integer :: i,j

    open(unit=21,file="EHP.inp",status="old")
    read(21,*) EHPpot
    close(21)

    do i=1,norbs,1
    do j=1,i-1,1
        if(bondlink(i,j)==1) then
            pppV(i,j)=EHPpot
            pppV(j,i)=EHPpot
            pppVlink(i,j)=1
            pppVlink(j,i)=1
        end if
    end do
    end do

    return
end subroutine Extented_Hubbard_Peierls_Potential 

!============================================================
!============================================================

!ohno_potential begins
subroutine Ohno_Potential

    implicit none
    ! local
    real(kind=r8) :: distance,U,midtmp  ! distance is atom-atom distance**2
    integer :: i,j,k
    
    call master_print_message("enter in PPP_term subroutine!")
    
    do j=1,norbs,1
        do i=j+1,norbs,1
            distance=0.0D0
            do k=1,3,1
                distance=(coord(k,i)-coord(k,j))*(coord(k,i)-coord(k,j))+distance
            end do
            U=(hubbardU(i)+hubbardU(j))/2.0D0
            midtmp=1.0D0+U*U*distance/14.397D0/14.397D0
            pppV(i,j)=U/sqrt(midtmp)
            pppV(j,i)=pppV(i,j)
            if(logic_Peierls==1 .and. ifpeierlsD==1) then
                pppw(i,j)=(sqrt(distance)*(U**3)/14.397D0/14.397D0)/(midtmp**1.5D0)
                pppw(j,i)=pppw(i,j)
            end if
            pppVlink(i,j)=1
            pppVlink(j,i)=1
        end do
    end do
    ! 14.397=e^2/4*pai*epsion0/e/angstrom
    return
end subroutine

!=========================================================  
end Module PPP_term_mod
 
