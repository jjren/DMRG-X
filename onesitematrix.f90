Module OnesiteMatrix
! this subroutine is construct the onesite operator matrix
! in sparse CSR format
    
    use variables
    use communicate
    use kinds_mod

    implicit none
    
    integer,parameter :: noperators=11

    real(kind=r8) :: onesitemat(4,noperators)            ! one site matrix in 4*4 basis 
    integer(kind=i4) :: osmrowindex(5,noperators),osmcolindex(4,noperators)
    ! onesitemat(:,:,x)
    ! x=1 means a(+)up
    ! x=2 means a(+)down
    ! x=3 means occpuation operator-nulcearQ
    ! x=4 means a up
    ! x=5 means a down
    ! x=6 means small Hamiltonian
    ! x=7 means niup+nidown operator
    ! x=8 means niup-nidown operator 
    ! x=9 means ai^+up*aidown
    ! x=10 means ai^+down*aiup
    ! x=11 means (niup-nidown)^2

    contains
!====================================================
!====================================================

subroutine ConstructOnesiteMatrix(orbindex)

    implicit none

    integer :: orbindex
    integer :: i
    ! orbindex is the orbital index
    ! this is the new site index maybe nleft+1 or nright-1
    real(kind=r8) :: zero

    call master_print_message("enter subroutine onesitematrix")
    
    onesitemat=0.0D0
    osmrowindex=0
    osmcolindex=0
    
    ! x=1 means a(+)up
    onesitemat(1,1)=1.0D0
    onesitemat(2,1)=1.0D0
    osmcolindex(1,1)=1
    osmcolindex(2,1)=3
    osmrowindex(1,1)=1
    osmrowindex(2,1)=1
    osmrowindex(3,1)=2
    osmrowindex(4,1)=2
    osmrowindex(5,1)=3

    ! x=2 means a(+)down
    onesitemat(1,2)=1.0D0
    onesitemat(2,2)=-1.0D0
    osmcolindex(1,2)=1
    osmcolindex(2,2)=2
    osmrowindex(1,2)=1
    osmrowindex(2,2)=1
    osmrowindex(3,2)=1
    osmrowindex(4,2)=2
    osmrowindex(5,2)=3

    ! x=3 means occpuation operator-nulcearQ
    ! save memory
    zero=1.0D-2
    if(abs(nuclQ(orbindex)-2.0D0) < zero) then
        onesitemat(1,3)=-2.0D0
        onesitemat(2,3)=-1.0D0
        onesitemat(3,3)=-1.0D0
        osmcolindex(1,3)=1
        osmcolindex(2,3)=2
        osmcolindex(3,3)=3
        osmrowindex(1,3)=1
        osmrowindex(2,3)=2
        osmrowindex(3,3)=3
        osmrowindex(4,3)=4
        osmrowindex(5,3)=4
    else if( abs(nuclQ(orbindex)-1.0D0) < zero ) then
        onesitemat(1,3)=-1.0D0
        onesitemat(2,3)=1.0D0
        osmcolindex(1,3)=1
        osmcolindex(2,3)=4
        osmrowindex(1,3)=1
        osmrowindex(2,3)=2
        osmrowindex(3,3)=2
        osmrowindex(4,3)=2
        osmrowindex(5,3)=3
    else
        onesitemat(1,3)=-nuclQ(orbindex)
        onesitemat(2,3)=1.0D0-nuclQ(orbindex)
        onesitemat(3,3)=1.0D0-nuclQ(orbindex)
        onesitemat(4,3)=2.0D0-nuclQ(orbindex)
        do i=1,4,1
            osmcolindex(i,3)=i
            osmrowindex(i,3)=i
        end do
        osmrowindex(5,3)=5
    end if

    ! x=4 means a up
    onesitemat(1,4)=1.0D0
    onesitemat(2,4)=1.0D0
    osmcolindex(1,4)=2
    osmcolindex(2,4)=4
    osmrowindex(1,4)=1
    osmrowindex(2,4)=2
    osmrowindex(3,4)=2
    osmrowindex(4,4)=3
    osmrowindex(5,4)=3

    ! x=5 means a down
    onesitemat(1,5)=1.0D0
    onesitemat(2,5)=-1.0D0
    osmcolindex(1,5)=3
    osmcolindex(2,5)=4
    osmrowindex(1,5)=1
    osmrowindex(2,5)=2
    osmrowindex(3,5)=3
    osmrowindex(4,5)=3
    osmrowindex(5,5)=3

    ! x=6 means small Hamiltonian
    onesitemat(1,6)=t(orbindex,orbindex)
    onesitemat(2,6)=t(orbindex,orbindex)
    onesitemat(3,6)=2.0D0*t(orbindex,orbindex)+HubbardU(orbindex)
    osmcolindex(1,6)=2
    osmcolindex(2,6)=3
    osmcolindex(3,6)=4
    osmrowindex(1,6)=1
    osmrowindex(2,6)=1
    osmrowindex(3,6)=2
    osmrowindex(4,6)=3
    osmrowindex(5,6)=4

    ! x=7 means niup+nidown operator  ; used in bond order calculation
    onesitemat(1,7)=1.0D0
    onesitemat(2,7)=1.0D0
    onesitemat(3,7)=2.0D0
    osmcolindex(1,7)=2
    osmcolindex(2,7)=3
    osmcolindex(3,7)=4
    osmrowindex(1,7)=1
    osmrowindex(2,7)=1
    osmrowindex(3,7)=2
    osmrowindex(4,7)=3
    osmrowindex(5,7)=4

    ! x=8 means niup-nidown operator  ; used in bond order calculation
    onesitemat(1,8)=1.0D0
    onesitemat(2,8)=-1.0D0
    osmcolindex(1,8)=2
    osmcolindex(2,8)=3
    osmrowindex(1,8)=1
    osmrowindex(2,8)=1
    osmrowindex(3,8)=2
    osmrowindex(4,8)=3
    osmrowindex(5,8)=3

    ! x=9 means ai^+up*aidown
    onesitemat(1,9)=1.0D0
    osmcolindex(1,9)=3
    osmrowindex(1,9)=1
    osmrowindex(2,9)=1
    osmrowindex(3,9)=2
    osmrowindex(4,9)=2
    osmrowindex(5,9)=2

    ! x=10 means ai^+down*aiup
    onesitemat(1,10)=1.0D0
    osmcolindex(1,10)=2
    osmrowindex(1,10)=1
    osmrowindex(2,10)=1
    osmrowindex(3,10)=1
    osmrowindex(4,10)=2
    osmrowindex(5,10)=2

    ! x=11 means (niup-nidown)^2
    onesitemat(1,11)=1.0D0
    onesitemat(2,11)=1.0D0
    osmcolindex(1,11)=2
    osmcolindex(2,11)=3
    osmrowindex(1,11)=1
    osmrowindex(2,11)=1
    osmrowindex(3,11)=2
    osmrowindex(4,11)=3
    osmrowindex(5,11)=3

return

end subroutine ConstructOnesiteMatrix

!====================================================
!====================================================
end module OnesiteMatrix
