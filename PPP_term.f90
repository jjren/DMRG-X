	Module PPP_term
!this module can include many potential terms like ohno potential
	Use Variables
	USE MPI
	implicit none
  
    
	contains

!ohno_potential begins------------------------
	subroutine ohno_potential

	implicit none
	real(kind=8) :: distance,U
	integer :: i,j,k
	
	if(myid==0) then 
		write(*,*) "enter in PPP_term subroutine!"
	end if
	
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
!	pppV=0.0D0
! 14.397=e^2/4*pai*epsion0/e/angstrom
	return
	end subroutine
! ohno_potential ends------------------------
  
	end Module
 
