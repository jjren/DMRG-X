Module transmoment_mod
	use variables
	use communicate
	use module_sparse

	implicit none

	real(kind=r8),allocatable :: moment(:,:)  ! store the last transition moment
	
	contains

!=======================================================================
!=======================================================================

subroutine transmoment
! this subroutine is to calculate the transition moment between gs and ex
! only used in the state average method
!=========================================================
! the transition operator is (sum(rii*ni))
!(Rxi,Ryi,Rzi) <R|<L|ni|L'>|R'>*CLR*CL'R' 
! if i in L space
! =<L|ni|L'>CLRCL'R
!=========================================================

	use exit_mod
	implicit none
	integer :: i

	call master_print_message("enter transmoment subroutine")
	
	if(nstate==1) then
		call exit_DMRG(sigAbort,"nstate==1 wrong!")
	end if
	if(Lrealdim/=subM .or. Rrealdim/=subM) then
		call exit_DMRG(sigAbort,"Lrealdim/=subM or Rrealdim/=subM")
	end if

	if(myid==0) then
		allocate(moment(4,nstate))
		moment=0.0D0
	end if

	call transmoment_subspace('L')
	call transmoment_subspace('R')

	if(myid==0) then
		write(*,*) "transition moment"
		do i=2,nstate,1
			moment(4,i)=sqrt(moment(1,i)**2+moment(2,i)**2+moment(3,i)**2)
			write(*,'(4F10.5)') moment(:,i)
		end do
	end if

return
end subroutine transmoment

!=======================================================================
!=======================================================================

subroutine transmoment_subspace(domain)
! calculate the left and right space transmoment contribution

	use mpi
	use blas95
	use f95_precision
	use exit_mod
	use module_sparse
	use mathlib

	implicit none
	include "mkl_spblas.fi"
	
	character(len=1) :: domain
	! local
	integer :: operaindex,orbstart,orbend
	integer :: status(MPI_STATUS_SIZE)
	real(kind=r8) :: localratio
	integer :: nmid
	real(kind=r8),allocatable :: midmat(:)
	integer(kind=i4),allocatable :: midcolindex(:),midrowindex(:)
	real(kind=r8),allocatable :: imoment(:)
	character(len=1) :: trans,transB
	! imoment(1) is a mediate number
	! imoment(2:nstate) is <faiex|ni|faigs>
	integer :: i,j,k,error,ierr,info
	

	! set the parameters
	if(domain=='L') then
		orbstart=1
		orbend=nleft+1
	else if(domain=='R') then
		orbstart=norbs-nright
		orbend=norbs
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
	end if
	
	
	if(myid/=0) then
		localratio=5.0
		nmid=CEILING(DBLE(16*subM*subM)/localratio)
		allocate(midmat(nmid),stat=error)
		if(error/=0) stop
		allocate(midcolindex(nmid),stat=error)
		if(error/=0) stop
		allocate(midrowindex(4*subM+1),stat=error)
		if(error/=0) stop
	end if

	! store the intermediate transmoment
	allocate(imoment(nstate),stat=error)
	if(error/=0) stop
	
	! operator ni loop
	do i=orbstart,orbend,1
		if(myid==orbid1(i,1)) then
			operaindex=orbid1(i,2)
		!	it does not need to remove the nuclearQ
		!	because the gs and ex state is orthogonal(ni-Qi)
			
			if(domain=='R') then
				transB='T'
			else 
				transB='N'
			end if

			call SpMMtoSp('N',transB,4*subM,4*subM,4*subM,4*subM,4*subM, &
					operamatbig1(:,3*operaindex),bigcolindex1(:,3*operaindex),bigrowindex1(:,3*operaindex), &
					coeffIF(:,1),coeffIFcolindex(:,1),coeffIFrowindex(:,1), &
					midmat,midcolindex,midrowindex,nmid)
			
			! the trace of coeffIF*transpose(midmat)
			imoment=0.0D0
			do j=2,nstate,1
				if(domain=='L') then
					trans='T'
				else if(domain=='R') then
					trans='N'
				end if
				call SpMMtrace(trans,4*subM,coeffIF(:,j),coeffIFcolindex(:,j),coeffIFrowindex(:,j),&
						midmat,midcolindex,midrowindex,imoment(j))
			end do
			
			call MPI_SEND(imoment,nstate,mpi_real8,0,i,MPI_COMM_WORLD,ierr)
		
		else if(myid==0) then

			call MPI_RECV(imoment,nstate,mpi_real8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
			! Rix/y/z *<ex|ni|gs>
			do j=2,nstate,1
				do k=1,3,1
					moment(k,j)=moment(k,j)+imoment(j)*coord(k,status(MPI_TAG))
				end do
			end do
		end if
	end do
	
	if(myid/=0) then
		deallocate(midmat,midcolindex,midrowindex)
	end if
return
end subroutine transmoment_subspace

!=======================================================================
!=======================================================================

end Module transmoment_mod
