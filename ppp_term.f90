Module PPP_term_mod
!this module can include many potential terms like ohno potential
    use variables, only : hubbardU,coord,pppV,norbs,logic_PPP,PPPpot
    use communicate
    use kinds_mod
    
    implicit none
    
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

    pppV=0.0D0
    do i=1,norbs,1
    do j=1,i-1,1
        if(bondlink(i,j)==1) then
            pppV(i,j)=EHPpot
            pppV(j,i)=EHPpot
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
    real(kind=r8) :: distance,U  ! distance is atom-atom distance**2
    integer :: i,j,k
    
    call master_print_message("enter in PPP_term subroutine!")
    
    do j=1,norbs,1
        do i=j+1,norbs,1
            distance=0.0D0
            do k=1,3,1
                distance=(coord(k,i)-coord(k,j))*(coord(k,i)-coord(k,j))+distance
            end do
            U=(hubbardU(i)+hubbardU(j))/2.0D0
            pppV(i,j)=U/sqrt(1+U*U*distance/14.397D0/14.397D0)
            pppV(j,i)=pppV(i,j)
        end do
    end do
    ! 14.397=e^2/4*pai*epsion0/e/angstrom
    return
end subroutine

!=========================================================  
end Module PPP_term_mod
 
