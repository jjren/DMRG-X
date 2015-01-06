Subroutine InitialGuess(guesscoeff)
! when nstate=1 then we can use the last stored matrix to
! contruct the Initial Guess
! only used in the finit MPS 
	use mpi
	use variables
	USE BLAS95
	USE F95_PRECISION

	implicit none
	real(kind=8) :: guesscoeff(16*Lrealdim*Rrealdim)
	real(kind=8),allocatable :: leftu(:,:),rightv(:,:),singularvalue(:)&
	,LRcoeff(:,:)
	logical :: alive
	integer :: reclength
	integer :: error,i

	if(myid==0) then
		write(*,*) "enter InitialGuess subroutine"
		! two site dmrg
		if((nright+nleft+2)/=norbs) then
			write(*,*) "-----------------------------------"
			write(*,*) "two site dmrg nright+nleft+2/=norbs"
			write(*,*) "-----------------------------------"
		end if

		allocate(leftu(4*Lrealdim,subM),stat=error)
		if(error/=0) stop
		allocate(rightv(subM,4*Rrealdim),stat=error)
		if(error/=0) stop
		allocate(singularvalue(subM),stat=error)
		if(error/=0) stop

		reclength=2*subM*subM

		inquire(file="wavefunction.tmp",exist=alive)
		if(alive) then
			open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=105,file="wavefunction.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		open(unit=106,file="singularvalue.tmp",status="old")
		
		do i=1,4,1
		read(105,rec=4*nleft+i) leftu((i-1)*Lrealdim+1:i*Lrealdim,1:subM)
		end do
		do i=1,4,1
		read(105,rec=4*(norbs-nright-1)+i) rightv(1:subM,i:4*Rrealdim:4)
		end do
		read(106,*) singularvalue(1:subM) 
		singularvalue=sqrt(singularvalue)

		do i=1,subM,1
			rightv(i,:)=rightv(i,:)*singularvalue(i)
		end do

		! recombine the two site sigmaL sigmaR coefficient
		
		allocate(LRcoeff(4*Lrealdim,4*Rrealdim),stat=error)
		if(error/=0) stop
		call gemm(leftu,rightv,LRcoeff,'N','N',1.0D0,0.0D0)
		do i=1,4*Rrealdim,1
			guesscoeff((i-1)*4*Lrealdim+1:i*4*Lrealdim)=LRcoeff(:,i)
		end do


		close(105)
		close(106)


		deallocate(leftu)
		deallocate(rightv)
		deallocate(singularvalue)
		deallocate(LRcoeff)
	end if
return
end subroutine
