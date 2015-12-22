subroutine coefftosparse(maxnelement,coeffmat,coeffmatcol,coeffmatrow,&
	num,coeffnosymm,iLrealdim,iRrealdim,cap_quantabigL,cap_quantabigR)
! this subroutine is to transfer the fortran column major coeff mat
! like coeffnosymm to CSR format
! the problem in that 
! CSR format and fortran's column major format not corresponds

! column major -> CSC format -> CSR format
	use variables
	use mathlib
	implicit none
	
	integer,intent(in) :: iLrealdim,iRrealdim,cap_quantabigL(4*iLrealdim,2),cap_quantabigR(4*iRrealdim,2)
	integer :: num,maxnelement
	real(kind=r8) :: coeffnosymm(num),coeffmat(maxnelement)
	integer(kind=i4) :: coeffmatrow(4*iLrealdim+1),coeffmatcol(maxnelement)

	! local
	integer(kind=i4),allocatable :: coeffrowdummy(:)
	integer :: error
	integer :: m,i,j

	if(maxnelement<num) then
		write(*,*) "==============================="
		write(*,*) "coefftosparse maxnelement<num",maxnelement,num
		write(*,*) "==============================="
		stop
	end if


	allocate(coeffrowdummy(4*iRrealdim+1),stat=error)
	if(error/=0) stop

	! in the CSC form
	m=0
	coeffrowdummy(1)=1

	do i=1,4*iRrealdim,1
	do j=1,4*iLrealdim,1
		if((cap_quantabigL(j,1)+cap_quantabigR(i,1)==nelecs) .and. &
			cap_quantabigL(j,2)+cap_quantabigR(i,2)==totalSz) then
			m=m+1
			coeffmat(m)=coeffnosymm(m)
			coeffmatcol(m)=j
		end if
	end do
	coeffrowdummy(i+1)=m+1
	end do

	! CSC transfer to CSR form
	call CSCtoCSR('CR',4*iLrealdim,4*iRrealdim,coeffmat,coeffmatcol,coeffrowdummy,coeffmatrow)

	deallocate(coeffrowdummy)

end subroutine coefftosparse
