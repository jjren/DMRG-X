subroutine sortcoeffLR(coeff,bigdim,smadim,nonzerocoeffdummy)
! nonzerocoeffdummy 

	USE mpi
	USE variables

	implicit none
	include "mkl_spblas.h"

	integer :: bigdim,smadim
	real(kind=8) :: coeff(bigdim,smadim)
	real(kind=8),allocatable :: LRcoeff(:,:,:)
	integer :: i,j,k,l
	real(kind=8) :: nonzerocoeffdummy(:,:)
	integer,allocatable:: columnindexdummy(:,:)
	integer :: job(8),info,rowindex(4*Lrealdim+1,smadim),upperdim,count1

	allocate(LRcoeff(4*Lrealdim,4*Rrealdim,smadim),stat=error)
	if(error/=0) stop
	allocate(nonzerocoeffdummy(16*Lrealdim*Rrealdim,smadim),stat=error)
	if(error/=0) stop
	allocate(columnindexdummy(16*Lrealdim*Rrealdim,smadim),stat=error)
	if(error/=0) stop

	
	job(1)=0
	job(2)=1
	job(3)=1
	job(4)=2
	job(5)=16*Lrealdim*Rrealdim
	job(6)=1

	
	upperdim=0

	do i=1,smadim,1
		do j=1,4*Rrealdim,1
			do k=1,4*Lrealdim,1
! convert the vecter to 4*Lrealdim*4*Rrealdim matrix
				LRcoeff(k,j,i)=coeff((j-1)*Lrealdim+k)
			end do
		end do
! convert the dense matrix to CSR format(3 array variation)
		call mkl_ddnscsr(job,4*Lrealdim,4*Rrealdim,LRcoeff(:,:,i),4*Lrealdim,nonzerocoeffdummy(:,i),columnindexdummy(:,i),rowindex,info)
		if(info/=0) then
			write(*,*) "coeff matrix dense to sparse(CSR format) failed !"
		end if

	end do

		
	allocate(nonzerocoeff(upperdim,smadim),stat=error)
	if(error/=0) stop
	allocate(columnindex(upperdim,smadim),stat=error)
	if(error/=0) stop

	do i=1,smadim,1
		nonzerocoeff(:,i)=nonzerocoeffdummy(1:upperdim,i)
		columnindex(:,i)=columnindexdummy(1:upperdim,i)
	end do



	
 






	








	deallocate(LRcoeff)




