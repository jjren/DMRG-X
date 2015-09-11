Subroutine Davidson_Wrapper(direction,lim,ilow,ihigh,iselec,niv,mblock,&
				   crite,critc,critr,ortho,maxiter)
! this is the davidson wrapper to call DVDSON written by Andreas
! aim to allocate memory and set variables in the global array
! mainly allocate memory on 0 process

	USE mpi
	USE variables
	USE InitialGuess
	USE symmetry
	use communicate
	use module_sparse

	implicit none

	character(len=1) :: direction
	integer :: lim,ilow,ihigh,niv,mblock,maxiter
	integer :: iselec(lim)
	real(kind=8) :: crite,critc,critr,ortho

	! local
	! davidson parameter
	integer :: dimN,nloops,nmv,ierror,smadim,IWRSZ
	logical :: hiend
	external op

	real(kind=r8),allocatable :: HDIAG(:),DavidWORK(:)
	real(kind=r8),allocatable :: dummycoeff(:),dummynewcoeff(:) ! have no use in fact
	real(kind=r8),allocatable :: nosymmout(:)
	integer :: reclength
	integer :: error,i,j,k,m
	
	integer :: ierr ! MPI flag

	call master_print_message("enter in davidson diagonalization subroutine")

	! check how many states fullfill good quantum number
	! every process do it
	! ngoodstates is the number of good quantum number states
	
	! allocate the symmetry work array
	if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
		call SymmAllocateArray
	end if

	ngoodstates=0
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
			ngoodstates=ngoodstates+1
			! construct the symmlinkgood 
			if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
				call CreatSymmlinkgood(ngoodstates,j,i)
			end if
		end if
	end do
		! construct the symmlinkcol
		! the nonzero LRcoeff element of every column
		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
			symmlinkcol(i+1)=ngoodstates+1
		end if
	end do

	dimN=ngoodstates
	
	if((logic_spinreversal/=0 .or. &
		(logic_C2/=0 .and. nleft==nright))) then
		! construct the symmetry matrix S in sparse format
		call SymmetryMat
		dimN=nsymmstate
	end if
	call master_print_message(ngoodstates,"ngoodstates:")
	call master_print_message(dimN,"total Hamiltonian dimension:")

!---------------------------------------------------
! can do direct diagonalization
!	call fullmat
!-----------------------------------------------------

! allocate the davidson workarray needed by DVDSON
	IWRSZ=lim*(2*dimN+lim+9)+lim*(lim+1)/2+nstate+100
	if(myid==0) then
		allocate(HDIAG(dimN),stat=error)
		if(error/=0) stop
		allocate(DavidWORK(IWRSZ),stat=error)
		if(error/=0) stop
	end if
	
! Get the diagonal element of hamiltonian
	if(logic_spinreversal/=0 .or. &
		(logic_C2/=0 .and. nleft==nright)) then
		call SymmHDiag(HDIAG)
	else 
		call GetHDiag(HDIAG)
	end if

! Get the Initialcoeff Guess
	if(myid==0) then
		call InitialStarter(direction,dimN,niv,DavidWORK)
	end if

!--------------------------------------------------------------------
! The core part of davidson diagnolization
	if(myid==0) then
		call DVDSON(op,dimN,lim,HDIAG,ilow,ihigh,iselec &
		    ,niv,mblock,crite,critc,critr,ortho,maxiter,DavidWORK,&
		    IWRSZ,hiend,nloops,nmv,ierror)
		    smadim=0
		call MPI_BCAST(smadim,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	end if

	if(myid/=0) then
		do while(.true.)
			call MPI_BCAST(smadim,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
			if(smadim>0) then
				allocate(dummycoeff(smadim),stat=error)
				if(error/=0) stop
				allocate(dummynewcoeff(smadim),stat=error)
				if(error/=0) stop
				call op(1,smadim,dummycoeff,dummynewcoeff)
				deallocate(dummycoeff)
				deallocate(dummynewcoeff)
			else
				exit
			end if
		end do
	end if

!-----------------------------------------------------------------------------
	
	if(myid==0) then
		if(hiend/=.false.) then
			call master_print_message("didn't get the lowest state")
			stop
		end if
		
		! transfer the symmetry state to the non-symmetry state S*fai
		allocate(nosymmout(ngoodstates*nstate),stat=error)
		if(error/=0) stop

		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
			do i=1,nstate,1
				call SymmetrizeState(ngoodstates,&
					nosymmout((i-1)*ngoodstates+1:i*ngoodstates),Davidwork((i-1)*nsymmstate+1:i*nsymmstate),'u')
			end do
		else
			nosymmout=DavidWORK(1:ngoodstates*nstate)
		end if
		!  no need do GramSchmit; the result fullfill orthnormal
		!  call GramSchmit(niv,ngoodstates,nosymmout,norm)
		!  write(*,*) "davidson nosymmout norm=",norm
		 
		! transfer the nosymmout to coeffIF
		do k=1,IHIGH,1
			call coefftosparse(4*Lrealdim,4*Rrealdim,&
				coeffIFdim,coeffIF(:,k),coeffIFcolindex(:,k),coeffIFrowindex(:,k),&
				ngoodstates,nosymmout((k-1)*ngoodstates+1:k*ngoodstates))
			coeffIFrowindex(4*Lrealdim+1:4*subM+1,k)=coeffIFrowindex(4*Lrealdim+1,k)
		end do

	!	write the coeffIF in two partical density matrix calculation

	!	reclength=nstate*coeffIFdim*2
	!	open(unit=109,file="coeffIF.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
	!	write(109,rec=1) coeffIF
	!	close(109)
	!	reclength=nstate*coeffIFdim
	!	open(unit=109,file="coeffIFcol.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
	!	write(109,rec=1) coeffIFcolindex
	!	close(109)
	!	reclength=(4*subM+1)*nstate
	!	open(unit=109,file="coeffIFrow.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
	!	write(109,rec=1) coeffIFrowindex
	!	close(109)

!=================================================================================
! write the final out

		write(*,*) "low state energy"
		do i=1,ihigh,1
			write(*,*) nleft+1,norbs-nright,i,"th energy=",DavidWORK(IHIGH*dimN+i)
			write(*,*) "energy converge:",DavidWORK(IHIGH*dimN+IHIGH+i)
			write(*,*) "residual norm:",DavidWORK(IHIGH*dimN+2*IHIGH+i)
		end do
		write(*,*) "NLOOPS=",nloops
		Write(*,*) "IERROR=",ierror
		write(*,*) "NMV=",nmv

		if(ierror/=0) then
			call master_print_message("caution! IERROR/=0")
			if(ierror/=2048) then
				call master_print_message("failed!")
				stop
			end if
		end if

!==================================================================================
		
! update the sweepenergy
! use the middle site as the sweepenergy
		if(nleft==(norbs+1)/2-1) then
			do i=1,ihigh,1
				sweepenergy(isweep,i)=DavidWORK(ihigh*dimN+i)
				dmrgenergy(i)=sweepenergy(isweep,i)
			end do
		end if
		
		deallocate(HDIAG)
		deallocate(DavidWORK)
		deallocate(nosymmout)
	end if

! deallocate symmetry workarray
	if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
		call DestroySymm
	end if

return
end subroutine Davidson_Wrapper
