subroutine op(bigdim,smadim,coeff,newcoeff)
! this is the core subroutine to calculate the S*H*S*C or H*C
! the parallel schema follow JCP 12 3174(2004) garnet chan

	use mpi
	use variables

	implicit none
	include "mkl_spblas.h"

	integer :: error,i,j,k
	integer :: bigdim,smadim
	! bigdim is the totaldim 16M*M,smadim is the davidson small block dimension
	! bigdim may be < 16M*M because we use spin reversal symmetry, and the dimension is smaller may be half
	! if groud state smadim=1
	! if gs+ex smadim may be >1
	real(kind=8) :: coeff(bigdim,smadim),newcoeff(bigdim,smadim)
	! coeff is the input coefficient
	! new coeff is H cross C result
	integer :: operaindex
	real(kind=8),allocatable ::nonzerocoeff(:,:),nonzerocoeffdummy(:,:),LRcoeff(:,:)
	integer(kind=4),allocatable :: coeffcolumnindex(:,:),coeffcolumnindexdummy(:,:)
	integer(kind=4) :: coeffrowindex(4*Lrealdim+1,smadim)
	integer :: job(8)
	integer(kind=4) :: nonzeronum
	real(kind=8),allocatable :: Rcomponentmat(:,:,:)


	if(error/=0) stop
	job(1)=0
	job(2)=1
	job(3)=1
	job(4)=2
	job(5)=16*Lrealdim*Rrealdim
	job(6)=1

	if( myid==0 ) then
		write(*,*) "enter H*C subroutine"
		
		if(bigdim/=Lrealdim*Rrealdim*16) then
			write(*,*) "-----------------------------------"
			write(*,*) "total bigdim/= 16*Lrealdim*Rrealdim"
			write(*,*) "-----------------------------------"
			stop
		end if

		if(logic_spinreversal/=0) then
			write(*,*) "spinreversal symmetry. first do S*C, then do H*(SC), at last S(+)*(HSC)"
			! call S*C subroutine
		else
			write(*,*) "do H*C"
		end do
!-----------------------------------------------------------------------------
!to sort the 16M*M coeff to 4M*4M(L*R) format coeff(16M^2,n) to coeff(4M,4M,n) 
!and convert it to CRS sparse format
		allocate(nonzerocoeffdummy(16*Lrealdim*Rrealdim,smadim),stat=error)
		if(error/=0) stop
		allocate(coeffcolumnindexdummy(16*Lrealdim*Rrealdim,smadim),stat=error)
		if(error/=0) stop
		allocate(LRcoeff(4*Lrealdim,4*Rrealdim,smadim),stat=error)
		if(error/=0) stop
		
		nonzeronum=0
		do i=1,smadim,1
			do j=1,4*Rrealdim,1
				do k=1,4*Lrealdim,1
! convert the vecter to 4*Lrealdim*4*Rrealdim matrix
					LRcoeff(k,j,i)=coeff((j-1)*Lrealdim+k)
				end do
			end do
! convert the dense matrix to CSR format(3 array variation)
			call mkl_ddnscsr(job,4*Lrealdim,4*Rrealdim,LRcoeff(:,:,i),4*Lrealdim,nonzerocoeffdummy(:,i),coeffcolumnindexdummy(:,i),coeffrowindex(:,i),info)
			if(info/=0) then
				write(*,*) "coeff matrix dense to sparse(CSR format) failed !"
				stop
			end if
			if(nonzeronum < coeffrowindex(4*Lrealdim+1,i)) then
				nonzeronum=coeffrowindex(4*Lrealdim+1,i)
			end if
		end do
		
		allocate(nonzerocoeff(nonzeronum,smadim),stat=error)
		if(error/=0) stop
		allocate(coeffcolumnindex(nonzeronum,smadim),stat=error)
		if(error/=0) stop

		nonzerocoeff(:,:)=nonzerocoeffdummy(1:nonzeronum,:)
		coeffcolumnindex(:,:)=coeffcolumnindexdummy(1:nonzeronum,:)

		deallocate(nonzerocoeffdummy)
		deallocate(coeffcolumnindexdummy)
	end if

		call MPI_bcast(nonzeronum,1,MPI_integer4,0,MPI_COMM_WORLD,ierr)
	
	if(myid/=0) then
		allocate(nonzerocoeff(nonzeronum,smadim),stat=error)
		if(error/=0) stop
		allocate(coeffcolumnindex(nonzeronum,smadim),stat=error)
		if(error/=0) stop
	end if
		call MPI_bcast(nonzerofcoeff,nonzeronum*smadim,MPI_real8,0,MPI_COMM_WORLD,ierr)
		call MPI_bcast(coeffcolumnindex,nonzeronum*smadim,MPI_integer4,0,MPI_COMM_WORLD,ierr)
		call MPI_bcast(coeffrowindex,(4*Lrealdim+1)*smadim,MPI_integer4,0,MPI_COMM_WORLD,ierr)
!------------------------------------------------
! vlr=Hlrl'r'*Cl'r'=sum(opt,l',r')=sum(Lopt,l') parity*Oll'*sum(Ropt,r') Orr'cl'r'
! the parallel shema is that 0 process bcast the coeff matrix to other process
! and 0 process gather the R space component matrix(sum of Ropt)
! and 0 process bcast the R space componet matrix
! and 0 process gather the result

	do i=norbs,norbs-nright,-1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			


			
			
			

			
			
	if(myid==0) then



end subroutine op





