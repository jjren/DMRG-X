subroutine coefftosparse(nrows,ncols,maxnelement,coeffmat,coeffmatcol,coeffmatrow,&
	num,coeffnosymm)
! this subroutine is to transfer the fortran column major coeff mat
! like coeffnosymm( only ngoodstates number ) to CSR format
! the problem in that 
! CSR format and fortran's column major format not corresponds

! column major -> CSC format -> CSR format
	use variables
	use mathlib
	implicit none
	
	integer :: num,nrows,ncols,maxnelement
	real(kind=r8) :: coeffnosymm(num),coeffmat(maxnelement)
	integer(kind=i4) :: coeffmatrow(nrows+1),coeffmatcol(maxnelement)

	! local
	integer(kind=i4),allocatable :: coeffrowdummy(:)
	integer :: error
	integer :: m,i,j

	if(nrows/=4*Lrealdim .or. ncols/=4*Rrealdim) then
		write(*,*) "==============================================="
		write(*,*) "coefftosparse subroutine nrows/=4*Lrealdim &
		ncols/=4*Rrealdim",nrows,ncols,4*Lrealdim,4*Rrealdim
		write(*,*) "==============================================="
		stop
	end if

	if(num/=ngoodstates) then
		write(*,*) "==============================="
		write(*,*) "coefftosparse num/=ngoodstates",num,ngoodstates
		write(*,*) "==============================="
		stop
	end if

	if(maxnelement<num) then
		write(*,*) "==============================="
		write(*,*) "coefftosparse maxnelement<num",maxnelement,num
		write(*,*) "==============================="
		stop
	end if


	allocate(coeffrowdummy(4*Rrealdim+1),stat=error)
	if(error/=0) stop

	! in the CSC form
	m=0
	coeffrowdummy(1)=1

	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
			quantabigL(j,2)+quantabigR(i,2)==totalSz) then
			m=m+1
			coeffmat(m)=coeffnosymm(m)
			coeffmatcol(m)=j
		end if
	end do
	coeffrowdummy(i+1)=m+1
	end do

	! CSC transfer to CSR form
	call CSCtoCSR('CR',4*Lrealdim,4*Rrealdim,coeffmat,coeffmatcol,coeffrowdummy,coeffmatrow)

	deallocate(coeffrowdummy)

end subroutine coefftosparse
