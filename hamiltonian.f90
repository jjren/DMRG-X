subroutine Hamiltonian(direction)
! this is the core subroutine
! contruct the hamiltonian on the fly and muliply it with coefficient of
! the wavefunction and transfter it to the 0 process

	use mpi
	use communicate
	use variables
	use symmetry
	use InitialGuess
	use stateOverlap

	implicit none
    
	! davidson input
	integer :: &
		lim   ,   &       ! the expanding small space in davidson iteration
		ilow  ,   &       ! index of lowest eigenpair
		ihigh ,   &       ! index of the highest eigenpair
		niv   ,   &       ! number of initial vector
		mblock,   &       ! number of vector to be targeted in each iteration
		maxiter           ! maxiter iterations
	real(kind=r8) :: &
		crite ,    &      ! convergence eigenvalue
		critc ,    &      ! convergence coefficiency
		critr ,    &      ! convergence residual vector norm
		ortho             ! orthogonality threshold
	integer,allocatable :: iselec(:)        ! selected eigenpair
	
	! davidson output 
	integer :: dimN,nloops,nmv,IERROR,smadim,IWRSZ,NUME
	logical :: hiend
    
	real(kind=r8),allocatable :: HDIAG(:),DavidWORK(:)
	real(kind=r8),allocatable :: dummycoeff(:),dummynewcoeff(:) ! have no use in fact
	real(kind=r8),allocatable :: nosymmout(:)
	integer             :: reclength
	integer             :: i,j,k,m
	integer             :: ierr             ! MPI flag
	integer             :: error
	character(len=1)    :: direction        ! i,l,r direction l is L space to R space sweep
	logical             :: davidFinished    ! there might be multiple davidson to be done
    
	external op
	
	call master_print_message("enter hamiltonian subroutine")

	! initialize workspace of davidson diagonalization
	! including calculate dimension, process symmetry,
	! allocate DavidWORK and get diagonal elements of H
	call initDavidWORK        

!	if(logic_fullmat == .true. ) then 
!		call fullmat   ! direct diagonalization
!	end if
                                
	davidFinished = .false.    
	do while(davidFinished == .false.)
		call getDavidParameters         ! initialize parameters
		call Davidson                   ! core part of davidson diagnolization
		if(myid == 0) then
			call processDavidOut        
			call printResults
		end if
		call sync               ! synchronize davidFinished, targetStateIndex, targetStateFlag
	end do
	
	call clean      ! deallocate arrays

contains
!================================================================
!================================================================
    
subroutine initDavidWORK
	
	lim = 20+nstate
	allocate(iselec(lim),stat=error)
	if(error/=0) stop    
	iselec = -1

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

! allocate the davidson workarray needed by DVDSON
	!lim = 20    ! only tmp value to calculate biggest spaced might be used
	IWRSZ=lim*(2*dimN+lim+9)+lim*(lim+1)/2+nstate+100
	if(myid==0) then
		allocate(HDIAG(dimN),stat=error)
		if(error/=0) stop
		allocate(DavidWORK(IWRSZ),stat=error)
		if(error/=0) stop
		DavidWORK = 0.0
	end if
  
! Get the diagonal element of hamiltonian
	if(logic_spinreversal/=0 .or. &
		(logic_C2/=0 .and. nleft==nright)) then
		call SymmHDiag(HDIAG)
	else 
		call GetHDiag(HDIAG)
	end if
end subroutine initDavidWORK
    
!================================================================
!================================================================

subroutine getDavidParameters
	! convergence threshholds and iteration parameters
	call getCrit
	maxiter = 400
	
	! state indices to be targeted in davidson diagonalization
	if(exscheme == 4 .and. startedStateSpecific) then
		if(targetStateFlag == 'none') then
			call master_print_message(formerStateIndex,'formerStateIndex = ')
			targetStateIndex = formerStateIndex
			targetStateFlag = 'trysame'   ! first try the same state index
		end if
		if(targetStateFlag == 'trysame' .or. targetStateFlag == 'tryhigher' .or. targetStateFlag == 'reachedmax') then
			call master_print_message(targetStateIndex,"!!!!Targetting state index")
			ilow = 0
			ihigh = 0
			NUME = targetStateIndex
			niv = NUME
			mblock = 1
			iselec(1) = targetStateIndex
			!write(*,*) "parameters",dimN,lim,ilow,ihigh,iselec(1),iselec(2), &
			!niv,mblock,crite,critc,critr,ortho,maxiter,&
			!DavidWORK,IWRSZ
		else if(targetStateFlag == 'trylower')  then
			call master_print_message(formerStateIndex,"!!!!Targetting states lower than former index")
			if(formerStateIndex/=targetStateIndex) then
				call master_print_message("formerStateIndex/=targetStateIndex error")
				stop
			end if
			ilow = 1
			ihigh = formerStateIndex
			NUME = ihigh
			niv = NUME
			mblock = ihigh
		else
			if(myid==0) then
				write(*,*) "unexpected targetStateFlag in subroutine getDavidParameters", targetStateFlag
			end if
			stop
		end if
	else
		call master_print_message(nstate,"Calculating states .LE. nstate:")
		ilow = 1
		ihigh = nstate
		mblock = nstate
		NUME = nstate
		niv = NUME
	end if
end subroutine getDavidParameters
    
!================================================================
!================================================================
    
subroutine Davidson
! The core part of davidson diagnolization
! Get the Initialcoeff Guess
	if(myid==0) then
		call InitialStarter(direction,dimN,niv,DavidWORK)
	end if
	
	if(myid==0) then
		!if(targetStateFlag=='tryhigher') then
		!    write(*,*) "davidWORK", DavidWORK
		!end if
		write(*,*) "Begin Running Davidson Diagonalization"
		call DVDSON(op,dimN,lim,HDIAG,ilow,ihigh,iselec &
		    ,niv,mblock,crite,critc,critr,ortho,maxiter,DavidWORK,&
		    IWRSZ,hiend,nloops,nmv,IERROR)
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
    
end subroutine Davidson

!================================================================
!================================================================

subroutine processDavidOut
	if(hiend/=.false.) then
		call master_print_message("didn't get the lowest state")
		stop
	end if
	! transfer the symmetry state to the non-symmetry state S*fai
	if(allocated(nosymmout)) then
		deallocate(nosymmout)
	end if
	allocate(nosymmout(ngoodstates*NUME),stat=error)
	if(error/=0) stop
	if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
		do i=1,NUME,1
			call SymmetrizeState(ngoodstates,nosymmout((i-1)*ngoodstates+1:i*ngoodstates),&
				Davidwork((i-1)*nsymmstate+1:i*nsymmstate),'u')
		end do
	else
		nosymmout=DavidWORK(1:ngoodstates*NUME)
	end if
	
	! calculate state overlap with state from former step (unsymmetrized state)
	if(exscheme == 4 .and. startedStateSpecific) then 
		call getStateOverlap(nosymmout, NUME, direction, IERROR)
	end if

	! transfer the nosymmout to coeffIF
	do k=1,NUME,1
		call coefftosparse(4*Lrealdim,4*Rrealdim,&
			coeffIFdim,coeffIF(:,k),coeffIFcolindex(:,k),coeffIFrowindex(:,k),&
			ngoodstates,nosymmout((k-1)*ngoodstates+1:k*ngoodstates))
		coeffIFrowindex(4*Lrealdim+1:4*subM+1,k)=coeffIFrowindex(4*Lrealdim+1,k)
	end do

!	write the coeffIF in two partical density matrix calculation

!	if(exscheme/=4) then
!		reclength=nstate*32*subM*subM
!	else
!		reclength=highestStateIndex*32*subM*subM
!	end if        
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
!
end subroutine processDavidOut

!================================================================
!================================================================

subroutine printResults
	write(*,*) "==============DAVIDSON RESULTS============="
	write(*,*) "position:",nleft+1,norbs-nright
	write(*,*) "direction = ",direction
	write(*,*) "NLOOPS=",nloops
	Write(*,*) "IERROR=",IERROR
	write(*,*) "NMV=",nmv
	write(*,*) "crite=",crite
	if(exscheme == 4 .and. startedStateSpecific) then  ! state specific calculation, computing overlaps
		
		if(targetStateFlag == 'getsame' .or. targetStateFlag == 'ngetsame' .or. &
			targetStateFlag == 'gethigher' .or. targetStateFlag == 'ngethigher') then
			write(*,*) "Energy of the", targetStateIndex, "state is", &
				DavidWORK(NUME*dimN+targetStateIndex)
			write(*,*) "energy converge:",DavidWORK(NUME*dimN+NUME+targetStateIndex)
			write(*,*) "residual norm:",DavidWORK(NUME*dimN+2*NUME+targetStateIndex)
			select case(targetStateFlag)
				case('getsame')
					write(*,*) "GET STATE! keeps SAME", targetStateIndex
					write(*,*) "with overlap =", stateOverlapValue(targetStateIndex)
					targetStateFlag = 'none'
					formerStateIndex = targetStateIndex
					formerStateEnergy = DavidWORK(NUME*dimN+targetStateIndex)
					davidFinished = .true.
				case('ngetsame')
					write(*,*) "STATE INDEX MAY NOT BE THE SAME",  targetStateIndex
					write(*,*) "with overlap =", stateOverlapValue(targetStateIndex)
					if(targetStateIndex /= 1) then
						write(*,*) "try lower states"
						targetStateFlag = 'trylower'
					else
						write(*,*) "try higher states"
						targetStateFlag = 'tryhigher'
						targetStateIndex = targetStateIndex + 1
					end if                    
				case('gethigher')
					write(*,*) "GET STATE! from HIGHER", targetStateIndex,&
						"from former index", formerStateIndex
					write(*,*) "with overlap =", stateOverlapValue(targetStateIndex)
					targetStateFlag = 'none'
					formerStateIndex = targetStateIndex
					formerStateEnergy = DavidWORK(NUME*dimN+targetStateIndex)
					davidFinished = .true.
				case('ngethigher')
					write(*,*) "STATE INDEX MAY NOT BE",  targetStateIndex
					write(*,*) "with overlap =", stateOverlapValue(targetStateIndex)
					write(*,*) "try higher states"
					targetStateFlag = 'tryhigher'
					targetStateIndex = targetStateIndex + 1    ! try one higher state
			end select
		else if(targetStateFlag == 'getlower' .or. targetStateFlag == 'ngetlower') then
			write(*,*) "Energy of lowest", NUME, "states:"
			do i=1,NUME,1
				write(*,*) i,"th energy=",DavidWORK(NUME*dimN+i)
				write(*,*) "energy converge:",DavidWORK(NUME*dimN+NUME+i)
				write(*,*) "residual norm:",DavidWORK(NUME*dimN+2*NUME+i)
			end do
			do i=1,NUME,1
				write(*,*) "overlap between davidson solution", i, "and initial guess is"
				write(*,*) stateOverlapValue(i)
			end do
			select case(targetStateFlag)
			case('getlower')
				write(*,*) "GET STATE! from LOWER", maxOverlapIndex,&
				           "from former index", formerStateIndex
				write(*,*) "with overlap =", maxOverlapValue
				targetStateFlag = 'none'
				formerStateIndex = maxOverlapIndex
				formerStateEnergy = DavidWORK(NUME*dimN+maxOverlapIndex)
				davidFinished = .true.
			case('ngetlower')
				write(*,*) "STATE INDEX MAY NOT BE BELOW",  targetStateIndex
				write(*,*) "with overlap =", maxOverlapValue
				write(*,*) "try higher states"
				targetStateFlag = 'tryhigher'
				targetStateIndex = targetStateIndex + 1    ! try one higher state
			end select
		else if(targetStateFlag == 'reachedmax') then
			write(*,*) "overlap of state",  highestStateIndex,&
			           "from last step is only", stateOverlapValue(highestStateIndex)
			write(*,*) "already tried all states with index less than", highestStateIndex
			write(*,*) "can't trace former state with overlap higher than", overlapThresh
			write(*,*) "highest overlap is", maxOverlapValue, "from state", maxOverlapIndex
			!write(*,*) "try to retrieve wavefunction from last step......"
			!call retrieveFormerState(direction,davidWORK,IWRSZ,NUME,dimN)
			if(maxOverlapValue <= 0.1) then
				write(*,*) "maxOverlapValue <= 0.1, we have lost the correct state" 
				stop
			else
				write(*,*) "using state ", maxOverlapIndex
				targetStateFlag = 'none'
				formerStateIndex = maxOverlapIndex
				davidFinished = .true.
			end if

		else if(targetStateFlag == 'stoptrying') then
			write(*,*) "When targetting high index state, NLOOPS>MAXITER. Stop trying..."
			write(*,*) "already tried all states with index less than", targetStateIndex - 1
			write(*,*) "can't trace former state with overlap higher than", overlapThresh
			write(*,*) "highest overlap is", maxOverlapValue, "from state", maxOverlapIndex
			!write(*,*) "try to retrieve wavefunction from last step......"
			!call retrieveFormerState(direction,davidWORK,IWRSZ,NUME,dimN)
			if(maxOverlapValue <= 0.1) then
				write(*,*) "maxOverlapValue <= 0.1, we have lost the correct state" 
				stop
			else
				write(*,*) "using state ", maxOverlapIndex
				targetStateFlag = 'none'
				formerStateIndex = maxOverlapIndex
				davidFinished = .true.
			end if
		else
			write(*,*) "unexpected targetStateFlag"
			stop
		end if

		if(targetStateIndex >= highestStateIndex) then  
			targetStateFlag = 'reachedmax'
		end if

	else   ! usual davidson (ground state or state average)
		write(*,*) "Energy of lowest", NUME, "states"
		do i=1,NUME,1
			write(*,*) i,"th energy=",DavidWORK(NUME*dimN+i)
			write(*,*) "energy converge:",DavidWORK(NUME*dimN+NUME+i)
			write(*,*) "residual norm:",DavidWORK(NUME*dimN+2*NUME+i)
		end do            
		davidFinished = .true.
	end if
	
	if(IERROR/=0 .and. (exscheme/=4 .or. targetStateFlag/='none')) then
		call master_print_message("failed! IERROR/=0")
		stop
	end if

! update the sweepenergy
	if(nleft==(norbs+1)/2-1) then
		if(exscheme == 4 .and. startedStateSpecific) then
			stateSpecificSweepEnergy(isweep) = DavidWORK(NUME*dimN+targetStateIndex)
			storedStateIndex(isweep) = formerStateIndex
		else
			do i=1,NUME,1
				sweepenergy(isweep,i)=DavidWORK(NUME*dimN+i)
			end do
		end if
	end if

end subroutine printResults

!================================================================
!================================================================
    
subroutine sync
	call MPI_bcast(formerStateIndex,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(targetStateIndex,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(targetStateFlag,10,MPI_character,0,MPI_COMM_WORLD,ierr)
	call MPI_bcast(davidFinished,1,MPI_logical,0,MPI_COMM_WORLD,ierr)
	DavidWORK(ngoodstates*NUME+1:)=0.0
end subroutine sync

!================================================================
!================================================================

subroutine clean
	if(myid == 0) then
		if(exscheme==4) then
			do i=1,highestStateIndex,1
				call copy(coeffIF(:,i),realcoeffIF(:,i))
			end do
			RealcoeffIFrowindex=coeffIFrowindex
			RealcoeffIFcolindex=coeffIFcolindex
		end if
		deallocate(HDIAG)
		deallocate(DavidWORK)
		deallocate(nosymmout)
	end if
	if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
		call DestroySymm   ! deallocate symmetry workarray
	end if
	deallocate(iselec)
end subroutine clean
    
!================================================================
!================================================================

subroutine getCrit
	
	implicit none
	
	real(kind=r8)     ::   originCrite,finalCrite,infiniteCrite
	integer           ::   pastSweep

	! default value
	critc = 1.0D-9
	critr = 1.0D-9
	ortho = 1.0D-8
	infiniteCrite = 1.0D-5
	originCrite = 1.0D-7
	finalCrite  = 1.0D-9
	
	if(direction == 'i') then  ! infinite DMRG
		crite = infiniteCrite
	else                       ! finite sweeps
		pastSweep = isweep - 1
		if(exscheme==4 .and. startedStateSpecific) then
			pastSweep = pastSweep - sweeps
		end if
		crite = originCrite * (0.2**pastSweep)
		if(crite <= 0.2*finalCrite) then
			reachedEnergyThresh = .true.
		end if
		if(crite <= finalCrite) then
			crite = finalCrite
		end if
	end if

	if(exscheme/=4) then
		! in the initial few sweeps
		if(1.0D-5*(1.0D-1)**isweep>1.1D-10) then
			crite=1.0D-5*(1.0D-1)**isweep
			critc=1.0D-5*(1.0D-1)**isweep
			critr=1.0D-5*(1.0D-1)**isweep
			ortho=1.0D-4*(1.0D-1)**isweep
		else
			crite=1.0D-9
			critc=1.0D-9
			critr=1.0D-9
			ortho=1.0D-8
		end if

		! in the last 3 sweeps
		if(isweep==MaxSweeps .or. isweep==MaxSweeps-1 .or. &
		isweep==MaxSweeps-2) then
			crite=1.0D-10
			critc=1.0D-10
			critr=1.0D-10
			ortho=1.0D-9
		end if
	end if
return

end subroutine getCrit

!================================================================
!================================================================

end subroutine Hamiltonian
