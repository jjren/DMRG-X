Module densesparse
	USE mpi
	USE variables
	
	implicit none

	contains
	(dense,m,n,nonzero,columnindex,rowindex)

	use mpi

	implicit none

	integer :: bigdim
	real(kind=8) coeff(bigdim)
