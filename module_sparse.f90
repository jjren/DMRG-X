module module_sparse
! this module contains the sparse format matrix in 3 array CSR format
! operamatbig and Hbig
! the core workarray of this program

	use kinds_mod
	use variables
	
	implicit none
	private
	save

	public :: AllocateArray

	! sparse form in 3 array CSR format
	real(kind=r8),allocatable,public :: &
	operamatbig1(:,:) , &         ! sparse form 1 electron operamatbig
	operamatsma1(:,:) , &         ! sparse form 1 electron operamatsma
	operamatbig2(:,:) , &         ! sparse form 2 electron operamatbig
	operamatsma2(:,:) , &         ! sparse form 2 electron operamatsma
	Hbig(:,:)         , &         ! Hbig in sparse form
	Hsma(:,:)                     ! Hsma in sparse form
	
	integer(kind=i4),allocatable,public :: &
	bigrowindex1(:,:) , &         ! 1 electron operamatbig rowindex
	bigcolindex1(:,:) , &         ! 1 electron oepramatbig columnindex
	smarowindex1(:,:) , &         ! 1 electron operamatsma rowindex
	smacolindex1(:,:) , &         ! 1 electron operamatsma columnindex
	bigrowindex2(:,:) , &         ! 2 electron operamatbig rowindex
	bigcolindex2(:,:) , &         ! 2 electron oepramatbig columnindex
	smarowindex2(:,:) , &         ! 2 electron operamatsma rowindex
	smacolindex2(:,:) , &         ! 2 electron operamatsma columnindex
	Hbigcolindex(:,:) , &         ! Hbig colindex
	Hbigrowindex(:,:) , &         ! Hbig rowindex
	Hsmacolindex(:,:) , &         ! Hsma colindex
	Hsmarowindex(:,:)             ! Hsma rowindex

	integer(kind=i4),public :: bigdim1,smadim1,bigdim2,smadim2,Hbigdim,Hsmadim  ! in sparse form operamatbig/operamatsma,Hbig/Hsma dim
	
	! parameter
	real(kind=r4),parameter,public :: pppmatratio=2.0,&
	hopmatratio=2.0,LRoutratio=2.0,UVmatratio=2.0
	
	real(kind=r4),parameter :: bigratio1=5.0,smaratio1=5.0,bigratio2=5.0,smaratio2=5.0,Hbigratio=5.0,Hsmaratio=5.0  ! sparse radio
	logical,parameter,public :: sparseform=.true.

	contains

!=========================================================================================================
!=========================================================================================================

subroutine AllocateArray(operanum1,operanum2)
	
	use communicate
	implicit none
	
	! store the number of operators on every process
	integer :: operanum1(nprocs-1),operanum2(nprocs-1)
	
	! local
	integer :: error

! set the sparse mat dim
	bigdim1=NINT(DBLE(16*subM*subM)/bigratio1,i4)
	bigdim2=NINT(DBLE(16*subM*subM)/bigratio2,i4)
	smadim1=NINT(DBLE(subM*subM)/smaratio1,i4)
	smadim2=NINT(DBLE(subM*subM)/smaratio2,i4)
	Hbigdim=NINT(DBLE(16*subM*subM)/Hbigratio,i4)
	Hsmadim=NINT(DBLE(subM*subM)/Hsmaratio,i4)

! allocate memory 
	if(myid/=0) then
		allocate(operamatbig1(bigdim1,3*operanum1(myid)),stat=error)
		if(error/=0) stop
		allocate(bigcolindex1(bigdim1,3*operanum1(myid)),stat=error)
		if(error/=0) stop
		allocate(bigrowindex1(4*subM+1,3*operanum1(myid)),stat=error)
		if(error/=0) stop

		allocate(operamatsma1(smadim1,3*operanum1(myid)),stat=error)
		if(error/=0) stop
		allocate(smacolindex1(smadim1,3*operanum1(myid)),stat=error)
		if(error/=0) stop
		allocate(smarowindex1(subM+1,3*operanum1(myid)),stat=error)
		if(error/=0) stop

		allocate(operamatbig2(bigdim2,2*operanum2(myid)),stat=error)
		if(error/=0) stop
		allocate(bigcolindex2(bigdim2,2*operanum2(myid)),stat=error)
		if(error/=0) stop
		allocate(bigrowindex2(4*subM+1,2*operanum2(myid)),stat=error)
		if(error/=0) stop

		allocate(operamatsma2(smadim2,2*operanum2(myid)),stat=error)
		if(error/=0) stop
		allocate(smacolindex2(smadim2,2*operanum2(myid)),stat=error)
		if(error/=0) stop
		allocate(smarowindex2(subM+1,2*operanum2(myid)),stat=error)
		if(error/=0) stop
		
		bigrowindex1=1   ! set the matrix to be 0
		smarowindex1=1
		bigrowindex2=1
		smarowindex2=1
	else
	! 2 means the R space ;1 means the L space
	
		allocate(Hbig(Hbigdim,2),stat=error)
		if(error/=0) stop
		allocate(Hbigcolindex(Hbigdim,2),stat=error)
		if(error/=0) stop
		allocate(Hbigrowindex(4*subM+1,2),stat=error)
		if(error/=0) stop

		allocate(Hsma(Hsmadim,2),stat=error)
		if(error/=0) stop
		allocate(Hsmacolindex(Hsmadim,2),stat=error)
		if(error/=0) stop
		allocate(Hsmarowindex(subM+1,2),stat=error)
		if(error/=0) stop
		
		Hbigrowindex=1
		Hsmarowindex=1
	end if

return

end subroutine AllocateArray

!=========================================================================================================
!=========================================================================================================

end module module_sparse
