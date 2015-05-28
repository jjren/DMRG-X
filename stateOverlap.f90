module stateOverlap
    use variables
    use Symmetry
    use InitialGuess
    use blas95
    use lapack95
    use kinds_mod
    
    implicit none
    
    integer         :: realTargetStateIndex           ! user specified index which to be traced
    integer(kind=4),allocatable &
                    :: storedStateIndex(:)            ! after every sweep is over the state index to be traced is stored
    real(kind=r8),allocatable &
                    :: stateSpecificSweepEnergy(:)    ! sweepenergy for state specific DMRG sweep
    real(kind=r8)   :: formerStateEnergy
    integer(kind=4) :: maxStateSpecificSweeps = 100   ! maximum number of max overlap sweeps
    integer(kind=4) :: maxStateSpecificSteps
    integer(kind=4) :: highestStateIndex              ! highest state considered when tracing excited state
    real(kind=8)    :: overlapThresh                  ! If two states have a bigger overlap than this value,
                                                      ! they are considered equal
    real(kind=8)    :: singularvalueThresh = 1.0D-3   ! determine whether two singular values are equal
    logical         :: printCorrectR   = .false.      ! whether print detail of correctR subroutine
    
    real(kind=r8),allocatable &
                    :: stateOverlapValue(:)
    real(kind=r8)   :: maxOverlapValue             ! working variable, to store maximum overlap value within a step
    integer         :: maxOverlapIndex             ! working variable, to store maximum overlap index within a step
    
    real(kind=8),allocatable &
                    :: smallOverlapValue(:)        ! When the overlap from last step is small(usually when smaller than 0.9),
    integer,allocatable &                          ! state specific DMRG calculation based on max-overlap algorithm will meet
                    :: smallOverlapiSweep(:)&      ! problems, e.g. traced an unexpected state. These variables stores the 
                      ,smallOverlapPosition1(:)&   ! information when small overlap value occurs during the state-specific
                      ,smallOverlapPosition2(:)    ! finit DMRG process
    character(len=1),allocatable &
                    :: smallOverlapDirection(:)
    real(kind=8)    :: smallOverlapThresh = 0.9 
    integer         :: smallOverlapCounter      
    
    real(kind=8)    :: alpha1                      ! when enter in subroutine retrieveFormerState, 
    real(kind=8)    :: alpha2                      ! coeffIF = alpha*coeffIF of this step + (1-alpha)*coeffIF of last step
                                                   ! alpha1 used when "reachedmax", alpha2 used when "stoptrying"
    contains
    subroutine initStateSpecific
        ! initialize state-specific DMRG 
        implicit none
        integer  :: error
        
        if(myid==0) then
            write(*,*) "**************************"
		    write(*,*) "enter in max overlap sweep"
            write(*,*) "**************************"
        end if
        
        formerStateIndex = realTargetStateIndex
        reachedEnergyThresh = .false.
        
        if(myid==0) then
            allocate(stateSpecificSweepEnergy(sweeps:sweeps+maxStateSpecificSweeps),stat=error)
            if(error/=0) stop
            stateSpecificSweepEnergy = 0.0
            write(*,*) "sweeps",sweeps
            write(*,*) "realTargetStateIndex",realTargetStateIndex
            stateSpecificSweepEnergy(sweeps) = sweepenergy(sweeps,realTargetStateIndex)
            allocate(storedStateIndex(sweeps+1:sweeps+maxStateSpecificSweeps),stat=error)
            if(error/=0) stop
            storedStateIndex=0
        
            maxStateSpecificSteps = stepPerSweep * maxStateSpecificSweeps
            allocate(smallOverlapValue(maxStateSpecificSteps),stat=error)
            if(error/=0) stop
            allocate(smallOverlapiSweep(maxStateSpecificSteps),stat=error)
            if(error/=0) stop
            allocate(smallOverlapPosition1(maxStateSpecificSteps),stat=error)
            if(error/=0) stop
            allocate(smallOverlapPosition2(maxStateSpecificSteps),stat=error)
            if(error/=0) stop        
            allocate(smallOverlapDirection(maxStateSpecificSteps),stat=error)
            if(error/=0) stop   
            smallOverlapCounter = 0
        end if
        
        call Renormalization(nleft+1,norbs-nright,'l')      !renormalize according to specific state
        
    end subroutine initStateSpecific

    subroutine getStateOverlap(nosymmout, ngoodstates, NUME, direction, ierror)
        implicit none
        integer           ::   NUME, ngoodstates, ierror
        character(len=1)  ::   direction
        character(len=10) ::   tmpTargetStateFlag
        real(kind=r8)     ::   nosymmout(ngoodstates*NUME)
        real(kind=r8)     ::   initial_vector(ngoodstates)
        integer           ::   m, n, i, j, k, error
        
        write(*,*) "enter subroutine getStateOverlap"
        
        !tmpTargetStateFlag = targetStateFlag
        !targetStateFlag = 'noretrieve'
        call SingleInitialFinite(initial_vector,ngoodstates,direction)  ! the initial vector got here is the wavefunction of last step
        !targetStateFlag = tmpTargetStateFlag
        
        ! special case 1: reaching highestStateIndex
        if(targetStateFlag == 'reachedmax') then
            stateOverlapValue(highestStateIndex) = dot(nosymmout((1+ngoodstates*(highestStateIndex-1)) : ngoodstates*highestStateIndex) &
                                                 , initial_vector(1:ngoodstates))
            stateOverlapValue(highestStateIndex) = abs(stateOverlapValue(highestStateIndex))
            if(stateOverlapValue(highestStateIndex)>=overlapThresh) then
                targetStateFlag = 'gethigher'
                maxOverlapValue = stateOverlapValue(highestStateIndex)
                maxOverlapIndex = highestStateIndex
            else
                maxOverlapValue = 0.0
                maxOverlapIndex = 0  
                do i=1,highestStateIndex,1
                    if (stateOverlapValue(i) > maxOverlapValue) then
                        maxOverlapValue = stateOverlapValue(i)
                        maxOverlapIndex = i
                    end if
                end do
            end if
        end if
        
        ! special case 2: NLOOPS>MAXITER when targetting high index state
        if(targetStateFlag=='tryhigher' .and. ierror==2048) then   
            targetStateFlag = 'stoptrying'
            maxOverlapValue = 0.0
            maxOverlapIndex = 0  
            do i=1,targetStateIndex - 1,1
                if (stateOverlapValue(i) > maxOverlapValue) then
                    maxOverlapValue = stateOverlapValue(i)
                    maxOverlapIndex = i
                end if
            end do
        end if
        
        ! commen cases
        select case(targetStateFlag)
        case('trysame')
            stateOverlapValue(targetStateIndex) = dot(nosymmout((1+ngoodstates*(targetStateIndex-1)) : ngoodstates*targetStateIndex) &
                                                 , initial_vector(1:ngoodstates))
            stateOverlapValue(targetStateIndex) = abs(stateOverlapValue(targetStateIndex))
            if(stateOverlapValue(targetStateIndex)>=overlapThresh) then
                targetStateFlag = 'getsame'
                maxOverlapValue = stateOverlapValue(targetStateIndex)
                maxOverlapIndex = targetStateIndex
            else
                targetStateFlag = 'ngetsame'
            end if
        case('trylower')
            maxOverlapValue = 0.0
            maxOverlapIndex = 0  
            do i=1,targetStateIndex,1
                stateOverlapValue(i) = dot(nosymmout((1+ngoodstates*(i-1)) : ngoodstates*i) &
                                      , initial_vector(1:ngoodstates))
                stateOverlapValue(i) = abs(stateOverlapValue(i))
                if (stateOverlapValue(i) > maxOverlapValue) then
                    maxOverlapValue = stateOverlapValue(i)
                    maxOverlapIndex = i
                end if
            end do
            if(maxOverlapValue>=overlapThresh) then
                targetStateFlag = 'getlower'
            else
                targetStateFlag = 'ngetlower'
            end if
        case('tryhigher')
            stateOverlapValue(targetStateIndex) = dot(nosymmout((1+ngoodstates*(targetStateIndex-1)) : ngoodstates*targetStateIndex) &
                                                 , initial_vector(1:ngoodstates))
            stateOverlapValue(targetStateIndex) = abs(stateOverlapValue(targetStateIndex))
            if(stateOverlapValue(targetStateIndex) >= overlapThresh) then
                targetStateFlag = 'gethigher'
                maxOverlapValue = stateOverlapValue(targetStateIndex)
                maxOverlapIndex = targetStateIndex
            else
                targetStateFlag = 'ngethigher'
            end if
!        case default
!            write(*,*) "unexpected targetStateFlag in subroutine getStateOverlap"
        end select
        
        ! if the overlap is small, record information about it
        if(targetStateFlag=='getsame' .or. targetStateFlag=='gethigher' .or. targetStateFlag=='getlower' .or. &
           targetStateFlag=='reachedmax' .or. targetStateFlag=='stoptrying')  then
            if(maxOverlapValue < smallOverlapThresh) then
                smallOverlapCounter = smallOverlapCounter + 1                  
                smallOverlapValue(smallOverlapCounter) = maxOverlapValue
                smallOverlapiSweep(smallOverlapCounter) = isweep
                smallOverlapPosition1(smallOverlapCounter) = nleft+1
                smallOverlapPosition2(smallOverlapCounter) = norbs-nright
                smallOverlapDirection(smallOverlapCounter) = direction
            end if
        end if

    end subroutine getStateOverlap
    
    subroutine correctR(singularvalue,leftu,rightv)
        implicit none
        real(kind=8)    ::   singularvalue(subM)
        real(kind=8)    ::   leftu(4*Lrealdim,subM),rightv(subM,4*Rrealdim)
        real(kind=8)    ::   matbuffer(subM,4*Rrealdim)        ! store the result of leftu(T)*coeffIF
        real(kind=8)    ::   rightvBuffer(subM,4*Rrealdim)     ! store the corrected rightv
        real(kind=8)    ::   S(subM,subM),S0(subM),absS(subM,subM)
        integer         ::   i,j,iGet,i1,i2,j1,j2
        integer         ::   SnonzeroNum,S0nonzeroNum,nonzeroNum,numDisorder,numNegative
        integer         ::   oldIndexToNew(subM),newIndexToOld(subM)
        logical         ::   usedIndex(subM)
        logical         ::   isDiagPositive

        write(*,*) "Enter in subroutine correctR"
        write(*,*) "S matrix is defined as (U+)*psi*V"
        write(*,*) "only consider singular values that are >=", singularvalueThresh
        
        do i=1, subM, 1
            S0(i) = sqrt(singularvalue(i))   !S0 stores the (true) singular value
        end do
        
        ! S stores the singular value matrix calculted by leftu^(T) * coeffIF * rightv
        call gemm(leftu,coeffIF(1:4*Lrealdim,1:4*Rrealdim,formerStateIndex),matbuffer,'T','N',1.0D0,0.0D0)
        call gemm(matbuffer,rightv,S,'N','T',1.0D0,0.0D0) 
        
        if(printCorrectR == .true.) then 
            write(*,*) "S matrix"
                do i=1, subM, 1
                    write(*,'(32F10.6)') S(i,:)
                end do
            write(*,*) "S0 vector"
            write(*,*) S0
        end if
        
        S0nonzeroNum = 0
        do i=1, subM, 1
            if(S0(i)>=singularvalueThresh) then
                S0nonzeroNum = S0nonzeroNum + 1
            end if
        end do
        write(*,*) "S0 has ", S0nonzeroNum, "nonzero elements"
        SnonzeroNum = 0
        do i=1, subM, 1
           do j=1, subM, 1
               if(abs(S(i,j))>singularvalueThresh) then
                   SnonzeroNum = SnonzeroNum + 1
               end if
           end do
        end do
        write(*,*) "S matrix has ", SnonzeroNum, "nonzero elements"
        nonzeroNum = min(S0nonzeroNum,SnonzeroNum)

        if(S0nonzeroNum/=SnonzeroNum) then
            write(*,*) "the difference may be due to rightv row with small singular value"
            write(*,*) "don't have corresponding leftu column (discarded)"
        end if

        oldIndexToNew = 0
        newIndexToOld = 0
        numDisorder = 0
        numNegative = 0
        usedIndex = .false.
        ! compare S with S0 to correct rightv
        do i=1, nonzeroNum, 1        ! only consider singular values that are bigger than singularvalueThresh
            do j=1, subM, 1
                if(abs(abs(S(i,j))-S0(i)) < 0.5*singularvalueThresh) then
                    if(i/=j) then
                        numDisorder = numDisorder + 1
                    end if
                    usedIndex(j) = .true.
                    if(S(i,j)>=0 ) then  ! ith row of correct rightv = jth row of uncorrected rightv
                        oldIndexToNew(j) = i
                        newIndexToOld(i) = j
                    else    ! find unmatched phase
                        oldIndexToNew(j) = -i
                        newIndexToOld(i) = -j
                        numNegative = numNegative +1
                    end if
                    exit
                end if
            end do
        end do
        do i=nonzeroNum+1,subM,1
            iGet = 0
            do j=1,subM,1
                if(usedIndex(j)==.false.) then
                    oldIndexToNew(j) = i
                    newIndexToOld(i) = j
                    usedIndex(j) = .true.
                    iGet = 1
                    exit
                end if
            end do
            if(iGet==0) then
                write(*,*) "iGet=0, cannot find any row of rightv+ unused"
                stop
            end if
        end do
        do j=1,subM,1
            if(usedIndex(j)==.false.) then
                write(*,*) "there are unused index"
                stop
            end if
        end do
        do j=1,subM,1
            if(newIndexToOld(j)==0) then
                write(*,*) "newIndexToOld=0 error"
                stop
            else if(oldIndexToNew(j)==0) then
                write(*,*) "oldIndexToNew=0 error"
                stop
            end if
        end do

        write(*,*) "Found",numDisorder,"disordered V+ rows"
        write(*,*) "Found",numNegative,"mismatched phase between column of U and row of V+"
        
        if(printCorrectR) then
            if(logic_spinreversal/=0) then
                write(*,*) "symmlinksma(:,1,2):"
                write(*,*) symmlinksma(:,1,2)
            end if
            write(*,*) "newIndexToOld:"
            write(*,*) newIndexToOld
            write(*,*) "quantasmaL"
            write(*,*) quantasmaL
            write(*,*) "quantasmaR"
            write(*,*) quantasmaR
        end if
        
        rightvBuffer = 0.0      
        do i=1, subM, 1
            if(newIndexToOld(i)>0) then
                rightvBuffer(i,:) = rightv(newIndexToOld(i),:)
            else if(newIndexToOld(i)<0) then
                rightvBuffer(i,:) = - rightv(-newIndexToOld(i),:)
            else
                write(*,*) "newIndexToOld=0 error"
                stop
            end if
        end do
        rightv= rightvBuffer     ! write corrected rightv
        
        isDiagPositive = .true.  ! test the new rightv
        call gemm(matbuffer,rightv,S,'N','T',1.0D0,0.0D0) 
        do i=1, subM, 1
            do j=1, subM, 1
                if(i==j) then
                    if(S(i,j)<-singularvalueThresh) then
                        isDiagPositive = .false.
                    end if
                else
                    if(abs(S(i,j))>singularvalueThresh) then
                        isDiagPositive = .false.
                    end if
                end if
            end do
        end do
            
        if(printCorrectR == .true.) then 
            write(*,*) "S matrix"
            do i=1, subM, 1
                write(*,'(32F10.6)') S(i,:)
            end do
        end if
            
        if(isDiagPositive == .true.) then
            write(*,*) "After correction, majority part of S matrix become positive and diagonal"
        else
            write(*,*) "After correction, majority part of S matrix is not positive and diagonal"
            stop
        end if

        if(logic_spinreversal/=0) then    ! correct spin matrix
            if(printCorrectR == .true.) then 
                write(*,*) "symmlinksma"
                write(*,*) symmlinksma(:,1,2)
            end if
            do i1=1,subM,1
                do i2=i1,subM,1
                    j1 = abs(newIndexToOld(i1))
                    j2 = abs(newIndexToOld(i2))
                    if(abs(symmlinksma(j1,1,2))==j2) then  ! if i1 and i2 is pair
                        symmlinksma(i1,1,2) = i2 * sign(1,newIndexToOld(i1)) &
                                                 * sign(1,symmlinksma(j1,1,2)) &
                                                 * sign(1,oldIndexToNew(j2))
                        symmlinksma(i2,1,2) = i1 * sign(1,newIndexToOld(i2)) &
                                                 * sign(1,symmlinksma(j2,1,2)) &
                                                 * sign(1,oldIndexToNew(j1))
                    end if
                end do
            end do
            if(printCorrectR == .true.) then 
                write(*,*) "after correction symmlinksma"
                write(*,*) symmlinksma(:,1,2)   
            end if
        end if
    end subroutine correctR
    
    subroutine checkStateSpecificResults
        ! check the results of state specific DMRG
        integer :: i
        if(storedStateIndex(isweep)/=realTargetStateIndex) then
            write(*,*) "#############################################################"
            write(*,*) "Warning: state index in the last step /= user specified index"
            write(*,*) "#############################################################"
        end if
        if(smallOverlapCounter/=0) then
            write(*,*) "#############################################################"
            write(*,*) "Warning: there are steps with overlap <", smallOverlapThresh
            write(*,*) "#############################################################"
            do i=1,smallOverlapCounter,1
                write(*,*) "At sweep",smallOverlapiSweep(i)
                write(*,*) "At position  ", smallOverlapPosition1(i),smallOverlapPosition2(i)
                write(*,*) "When direction is  ", smallOverlapDirection(i)
                write(*,*) "OVERLAP =", smallOverlapValue(i)
                write(*,*) "=========================================================="
            end do
        end if
    end subroutine checkStateSpecificResults
    
    subroutine cleanStateSpecificVariables
        if(myid==0) then
            deallocate(stateSpecificSweepEnergy)
            deallocate(storedStateIndex)
            deallocate(smallOverlapValue)
            deallocate(smallOverlapiSweep)
            deallocate(smallOverlapPosition1)
            deallocate(smallOverlapPosition2)
            deallocate(smallOverlapDirection)
        end if
    end subroutine cleanStateSpecificVariables
    
    subroutine retrieveFormerState(direction,davidWORK,IWRSZ,NUME,dimN)
        implicit none
        real(kind=8), allocatable   :: tmpCoeffIF(:,:)
        real(kind=8), allocatable   :: nosymmguess(:)
        character(len=1)  :: direction 
        integer           :: IWRSZ,NUME,dimN
        real(kind=r8)     :: davidWORK(IWRSZ)
        integer           :: i,j,m,error

        write(*,*) "enter in subroutine retrieveFormerState"
        write(*,*) "coeffIF = alpha*coeffIF of this step + (1-alpha)*coeffIF of last step"

        allocate(nosymmguess(ngoodstates),stat=error)
        if(error/=0) stop
        call SingleInitialFinite(nosymmguess,ngoodstates,direction)
        
        allocate(tmpCoeffIF(4*subM,4*subM),stat=error)
        if(error/=0) stop
        tmpCoeffIF=0.0D0
        m=1
	    do i=1,4*Rrealdim,1
	        do j=1,4*Lrealdim,1
		        if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
			        quantabigL(j,2)+quantabigR(i,2)==totalSz) then
			        tmpCoeffIF(j,i) = nosymmguess(m)!=LRcoeff(j,i)
			        m=m+1
		        end if
	        end do
        end do
        if(m-1/=ngoodstates) then
            write(*,*) "m-1/=ngoodstates error in subroutine retrieveFormerState"
            stop
        end if
        
        if(targetStateFlag=='reachedmax') then
            write(*,*) "alpha = alpha1 =", alpha1 
            coeffIF(:,:,formerStateIndex) = alpha1 * coeffIF(:,:,formerStateIndex) &
                                            + (1 - alpha1) * tmpCoeffIF(:,:)
        else if(targetStateFlag=='stoptrying') then
            write(*,*) "alpha = alpha2 =", alpha2
            coeffIF(:,:,formerStateIndex) = alpha2 * coeffIF(:,:,formerStateIndex) &
                                            + (1 - alpha2) * tmpCoeffIF(:,:)
        else
            write(*,*) "unexpected targetStateFlag in subroutine retrieveFormerState"
        end if        
        
        targetStateIndex = formerStateIndex
        davidWORK(NUME*dimN+targetStateIndex) = formerStateEnergy
        write(*,*) "retrieved last state wavefunction with energy=", formerStateEnergy
        
        deallocate(tmpCoeffIF)
        deallocate(nosymmguess)
        
    end subroutine retrieveFormerState
        
end module stateOverlap