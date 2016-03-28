Module CoeffTrans
    use kinds_mod
    use communicate
    implicit none
    save
    private

    public :: coefftosparse,Dense16LRtoNgood

contains
!========================================================================
!========================================================================

subroutine coefftosparse(maxnelement,coeffmat,coeffmatcol,coeffmatrow,&
    nosymmdim,coeffnosymm,iLrealdim,iRrealdim,cap_goodbasis,cap_goodbasiscol)
! this subroutine is to transfer the fortran column major coeff mat
! like coeffnosymm to CSR format
! the problem in that 
! CSR format and fortran's column major format not corresponds

! column major -> CSC format -> CSR format
    use variables
    use mathlib
    implicit none
    
    integer,intent(in) :: iLrealdim,iRrealdim,nosymmdim,maxnelement
    real(kind=r8),intent(in) :: coeffnosymm(nosymmdim)
    integer(kind=i4),intent(in) :: cap_goodbasiscol(4*iRrealdim+1),cap_goodbasis(cap_goodbasiscol(4*iRrealdim+1)-1,2)
    real(kind=r8),intent(out) :: coeffmat(maxnelement)
    integer(kind=i4),intent(out) :: coeffmatrow(4*iLrealdim+1),coeffmatcol(maxnelement)

    ! local
    integer :: i

    if(maxnelement<nosymmdim) then
        write(*,*) "==============================="
        write(*,*) "coefftosparse maxnelement<nosymmdim",maxnelement,nosymmdim
        write(*,*) "==============================="
        stop
    end if
    if(nosymmdim/=cap_goodbasiscol(4*iRrealdim+1)-1) then
        write(*,*) "nosymmdim/=cap_goodbasiscol(4*iRrealdim+1)-1",nosymmdim,cap_goodbasiscol(4*iRrealdim+1)-1
        stop
    end if

    ! in the CSC form
    call copy(coeffnosymm(1:nosymmdim),coeffmat(1:nosymmdim))
    
    ! copy the every column non zero nelements
    coeffmatcol(1:nosymmdim)=cap_goodbasis(1:nosymmdim,1)
    !forall(i=1:nosymmdim)
    !    coeffmatcol(i)=cap_goodbasis(i,1)
    !end forall

    ! CSC transfer to CSR form
    call CSCtoCSR('CR',4*iLrealdim,4*iRrealdim,coeffmat,coeffmatcol,cap_goodbasiscol,coeffmatrow)

end subroutine coefftosparse

!========================================================================
!========================================================================

subroutine Dense16LRtoNgood(iLrealdim,iRrealdim,nbasis,cap_goodbasis,LRcoeff,coeffnosymm)
    implicit none
    integer(kind=i4),intent(in) :: iLrealdim,iRrealdim,nbasis,cap_goodbasis(:,:)
    real(kind=r8),intent(in) :: LRcoeff(:)
    real(kind=r8),intent(out) :: coeffnosymm(:)
    !local
    integer :: ibasis

    do ibasis=1,nbasis,1
        coeffnosymm(ibasis)=LRcoeff((cap_goodbasis(ibasis,2)-1)*4*iLrealdim+cap_goodbasis(ibasis,1))
    end do

    return
end subroutine Dense16LRtoNgood

!========================================================================
!========================================================================
end Module CoeffTrans
