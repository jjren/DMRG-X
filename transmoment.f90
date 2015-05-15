subroutine transmoment
! this subroutine is to calculate the transition moment between gs and ex
! only used in the state average method
!=========================================================
! the transition operator is (sum(rii*ni))
!(Rxi,Ryi,Rzi) <R|<L|ni|L'>|R'>*CLR*CL'R' 
! if i in L space
! =<L|ni|L'>CLRCL'R
!=========================================================

	use mpi
	use variables
	use communicate
	use module_sparse
	use exit_mod
	use blas95
	use f95_precision


	implicit none
	include "mkl_spblas.fi"

	integer :: operaindex
	integer :: i,j,k,error,ierr
	integer :: status(MPI_STATUS_SIZE)
	real(kind=8),allocatable :: midmat(:,:)
	real(kind=8),allocatable :: imoment(:),zeromoment(:,:)
	! imoment(1) is a mediate number
	! imoment(2:nstate) is <faiex|ni|faigs>
	real(kind=8) :: moment(4,nstate)  ! store the last transition moment
	character(len=1) :: matdescra(6)
	logical :: iftransposed
	
	call master_print_message("enter transmoment subroutine")

	if(nstate==1) then
		call exit_DMRG(sigAbort,"nstate==1 wrong!")
	end if
	if(Lrealdim/=subM .or. Rrealdim/=subM) then
		call exit_DMRG(sigAbort,"Lrealdim/=subM or Rrealdim/=subM")
	end if
	
	matdescra(1)='G'
	matdescra(2)='L'
	matdescra(3)='N'
	matdescra(4)='F'

	if(myid/=0) then
		allocate(coeffIF(4*subM,4*subM,nstate),stat=error)  ! other process need coeffIF
		if(error/=0) stop
		allocate(midmat(4*subM,4*subM),stat=error)
		if(error/=0) stop
		allocate(imoment(nstate),stat=error)
		if(error/=0) stop
	else
		allocate(zeromoment(nstate,norbs),stat=error)
		if(error/=0) stop
		moment=0.0D0
	end if
	
	! 0 process broadcast the coeffIF
	call MPI_BCAST(coeffIF,16*subM*subM*nstate,MPI_real8,0,MPI_COMM_WORLD,ierr)
	iftransposed=.false.

	! operator ni loop
	do i=1,norbs,1
		if(myid==orbid1(i,1)) then
			operaindex=orbid1(i,2)
		!	it does not need to remove the nuclearQ
		!	because the gs and ex state is orthogonal(ni-Qi)
			if(i<=nleft+1) then
				call mkl_dcsrmm('N',4*subM,4*subM,4*subM,1.0D0,matdescra, &
					operamatbig1(:,3*operaindex),bigcolindex1(:,3*operaindex),bigrowindex1(1:4*subM,3*operaindex), &
					bigrowindex1(2:4*subM+1,3*operaindex),coeffIF(:,:,1),4*subM,0.0D0,midmat,4*subM)
				
				! the trace of coeffIF*transpose(midmat)
				imoment=0.0D0
				do j=2,nstate,1
					do k=1,4*subM,1
						imoment(1)=dot(coeffIf(k,:,j),midmat(k,:))
						imoment(j)=imoment(j)+imoment(1)
					end do
				end do
			else
				! in R space transpose the coeffIF(:,:,1)
				if(iftransposed==.false.) then
					coeffIF(:,:,1)=transpose(coeffIF(:,:,1))
					iftransposed=.true.
				end if

				call mkl_dcsrmm('N',4*subM,4*subM,4*subM,1.0D0,matdescra, &
					operamatbig1(:,3*operaindex),bigcolindex1(:,3*operaindex),bigrowindex1(1:4*subM,3*operaindex), &
					bigrowindex1(2:4*subM+1,3*operaindex),coeffIF(:,:,1),4*subM,0.0D0,midmat,4*subM)
				
				imoment=0.0D0
				do j=2,nstate,1
					do k=1,4*subM,1
						imoment(1)=dot(coeffIf(k,:,j),midmat(:,k))
						imoment(j)=imoment(j)+imoment(1)
					end do
				end do
			end if
			call MPI_SEND(imoment,nstate,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		
		else if(myid==0) then
			call MPI_RECV(zeromoment(:,i),nstate,mpi_real8,orbid1(i,1),i,MPI_COMM_WORLD,status,ierr)
			! Rix/y/z *<ex|ni|gs>
			do j=2,nstate,1
				do k=1,3,1
					moment(k,j)=moment(k,j)+zeromoment(j,i)*coord(k,i)
				end do
			end do
		end if
	end do
	
	if(myid==0) then
		write(*,*) "transition moment"
		do i=2,nstate,1
			moment(4,i)=sqrt(moment(1,i)**2+moment(2,i)**2+moment(3,i)**2)
			write(*,'(4F10.5)') moment(:,i)
		end do
	end if


end subroutine
