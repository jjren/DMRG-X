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
	Hsma(:,:)         , &         ! Hsma in sparse form
	coeffIF(:,:)      , &         ! coeffIF is the inital and final wavefunction coefficient 
	realcoeffIF(:,:)              ! coeffIF is the inital and final wavefunction coefficient 
	
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
	Hsmarowindex(:,:) , &         ! Hsma rowindex
	coeffIFcolindex(:,:) ,&       ! coeffIF colindex
	coeffIFrowindex(:,:) ,&       ! coeffIF rowindex
	realcoeffIFcolindex(:,:) ,&       ! realcoeffIF colindex
	realcoeffIFrowindex(:,:)          ! realcoeffIF rowindex

	integer(kind=i4),public :: bigdim1,smadim1,bigdim2,smadim2,Hbigdim,Hsmadim,coeffIFdim  ! in sparse form operamatbig/operamatsma,Hbig/Hsma dim
	
	! sparse parameter
	real(kind=r8),public :: pppmatratio,hopmatratio,LRoutratio,UVmatratio,coeffIFratio
	real(kind=r8),public :: bigratio1,smaratio1,bigratio2,smaratio2,Hbigratio,Hsmaratio  ! sparse radio

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
	
	call sparse_default

! set the sparse mat dim
	bigdim1=CEILING(DBLE(16*subM*subM)/bigratio1)
	bigdim2=CEILING(DBLE(16*subM*subM)/bigratio2)
	smadim1=CEILING(DBLE(subM*subM)/smaratio1)
	smadim2=CEILING(DBLE(subM*subM)/smaratio2)
	Hbigdim=CEILING(DBLE(16*subM*subM)/Hbigratio)
	Hsmadim=CEILING(DBLE(subM*subM)/Hsmaratio)
	coeffIFdim=CEILING(DBLE(16*subM*subM)/coeffIFratio)

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
		bigrowindex1=1   ! set the matrix to be 0
		smarowindex1=1
		
		if(logic_bondorder==1) then
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
			bigrowindex2=1
			smarowindex2=1
		end if
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
		if(exscheme/=4) then
			allocate(coeffIF(coeffIFdim,nstate),stat=error)
			if(error/=0) stop
			allocate(coeffIFcolindex(coeffIFdim,nstate),stat=error)
			if(error/=0) stop
			allocate(coeffIFrowindex(4*subM+1,nstate),stat=error)
			if(error/=0) stop
		else
			if(highestStateIndex<nstate) then
				write(*,*) "highestStateIndex<nstate error!" 
				stop
			end if
			allocate(coeffIF(coeffIFdim,highestStateIndex),stat=error)
			if(error/=0) stop
			allocate(coeffIFcolindex(coeffIFdim,highestStateIndex),stat=error)
			if(error/=0) stop
			allocate(coeffIFrowindex(4*subM+1,highestStateIndex),stat=error)
			if(error/=0) stop
			allocate(RealcoeffIF(coeffIFdim,highestStateIndex),stat=error)
			if(error/=0) stop
			allocate(RealcoeffIFcolindex(coeffIFdim,highestStateIndex),stat=error)
			if(error/=0) stop
			allocate(RealcoeffIFrowindex(4*subM+1,highestStateIndex),stat=error)
			if(error/=0) stop
		end if
		
		Hbigrowindex=1
		Hsmarowindex=1
		coeffIFrowindex=1
	end if

return

end subroutine AllocateArray

!=========================================================================================================
!=========================================================================================================

subroutine sparse_default
! set the default ratio according to the subM
	use communicate
	implicit none

	if(subM<100) then
		bigratio1=1.0
		smaratio1=1.0
		bigratio2=1.0
		smaratio2=1.0
		Hbigratio=1.0
		Hsmaratio=1.0
		pppmatratio=1.0
		hopmatratio=1.0
		LRoutratio=1.0
		UVmatratio=1.0
		coeffIFratio=1.0
	else if(abs(subM-120)<20) then
		bigratio1=32.0
		smaratio1=8.0
		bigratio2=32.0
		smaratio2=8.0
		Hbigratio=12.0
		Hsmaratio=8.0
		pppmatratio=10.0
		hopmatratio=12.0
		LRoutratio=10.0
		UVmatratio=10.0
		coeffIFratio=10.0
	else if (abs(subM-256)<50) then
		bigratio1=20.0
		smaratio1=5.0
		bigratio2=20.0
		smaratio2=5.0
		Hbigratio=10.0
		Hsmaratio=5.0
		pppmatratio=10.0
		hopmatratio=12.0
		LRoutratio=10.0
		UVmatratio=10.0
		coeffIFratio=10.0
	else if (abs(subM-512)<50) then
		bigratio1=25.0
		smaratio1=5.0
		bigratio2=25.0
		smaratio2=5.0
		Hbigratio=10.0
		Hsmaratio=5.0
		pppmatratio=10.0
		hopmatratio=12.0
		LRoutratio=10.0
		UVmatratio=10.0
		coeffIFratio=10.0
	else if (abs(subM-1024)<100) then
		bigratio1=25.0
		smaratio1=5.0
		bigratio2=25.0
		smaratio2=5.0
		Hbigratio=10.0
		Hsmaratio=5.0
		pppmatratio=10.0
		hopmatratio=12.0
		LRoutratio=10.0
		UVmatratio=10.0
		coeffIFratio=10.0
	else if (abs(subM-2048)<200) then
		bigratio1=30.0
		smaratio1=8.0
		bigratio2=30.0
		smaratio2=8.0
		Hbigratio=12.0
		Hsmaratio=8.0
		pppmatratio=12.0
		hopmatratio=12.0
		LRoutratio=10.0
		UVmatratio=10.0
		coeffIFratio=10.0
	else
		bigratio1=10.0
		smaratio1=5.0
		bigratio2=10.0
		smaratio2=5.0
		Hbigratio=10.0
		Hsmaratio=2.0
		pppmatratio=5.0
		hopmatratio=5.0
		LRoutratio=5.0
		UVmatratio=5.0
		coeffIFratio=5.0
	end if

	if(myid==0) then
		write(*,*) "bigratio1=",    bigratio1
		write(*,*) "smaratio1=",    smaratio1
		write(*,*) "bigratio2=",    bigratio2
		write(*,*) "smaratio2=",    smaratio2
		write(*,*) "Hbigratio=",    Hbigratio
		write(*,*) "Hsmaratio=",    Hsmaratio
		write(*,*) "pppmatratio=",  pppmatratio
		write(*,*) "hopmatratio=",  hopmatratio
		write(*,*) "LRoutratio=" ,  LRoutratio
		write(*,*) "UVmatratio=" ,  UVmatratio
		write(*,*) "coeffIFratio=", coeffIFratio
	end if
return

end subroutine sparse_default

!=========================================================================================================
!=========================================================================================================

end module module_sparse
