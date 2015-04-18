Module PPP_term
!this module can include many potential terms like ohno potential
	use variables, only : hubbardU,coord,pppV,norbs
	use communicate
	use kinds_mod
	
	implicit none
    
	contains

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
!	pppV=0.0D0
! 14.397=e^2/4*pai*epsion0/e/angstrom
	return
end subroutine

!=========================================================  
end Module PPP_term
 
