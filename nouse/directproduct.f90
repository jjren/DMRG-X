Subroutine directproduct(a,dima,b,dimb,c,phase)
! this subroutine is to do direct product of two matrix and there will
! be a phase
! for example :<sigmaL|<L|O|L>|sigmaL>
! a is the L matrix, b is the sigmaL matrix, phase is the antisymmetry
! phase

USE variables
USE mpi

implicit none

integer(kind=4) :: dima,dimb
real(kind=8) ::a(dima,dima),b(dimb,dimb),c(dima*dimb,dima*dimb)
integer(kind=4),optional :: phase(dima*dimb,dima*dimb)
integer :: i,j



c=0.0D0

do i=1,dimb,1
	do j=1,dimb,1
		c((j-1)*dima+1:j*dima,(i-1)*dima+1:i*dima)=a*b(j,i)
	end do
end do

if(present(phase)) then
do i=1,dima*dimb,1
	do j=1,dima*dimb,1
		c(j,i)=c(j,i)*DBLE(phase(j,i))
	end do
end do
end if


return
end Subroutine directproduct
