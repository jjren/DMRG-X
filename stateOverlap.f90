module stateOverlap
    use variables
    use Symmetry
    use InitialGuess
    use blas95
    use lapack95
    
    real(kind=8),allocatable  &
                    :: stateOverlapValue(:)
    integer(kind=4) :: maxOverlapSweeps = 100              ! maximum iteration of max overlap sweeps
    integer(kind=4) :: highestStateIndex = 6               ! highest state considered when tracing excited state
    real(kind=8)    :: overlapThresh = 0.6                 ! if two states have a bigger overlap than this value,
                                                           ! they are considered equal
    real(kind=8)    :: singularvalueThresh   = 1.0D-3      ! determine whether two singular values are equal
    logical         :: printSMat             = .false.     ! whether print S matrix
    
    real(kind=r8)::   maxOverlapValue
    integer      ::   maxOverlapIndex
    
    contains
    subroutine getStateOverlap(nosymmout, ngoodstates, NUME, direction)
        integer      ::   NUME, ngoodstates
        character(len=1)  ::   direction
        real(kind=r8)     ::   nosymmout(ngoodstates*NUME)
        real(kind=r8)     ::   initial_vector(ngoodstates)
        integer      ::   m, n, i, j, k, error
        
        call SingleInitialFinite(initial_vector,ngoodstates,direction)  ! the initial vector got here is the wavefunction of last step
        
        select case(targetStateFlag)
        case('trysame')
            stateOverlapValue(targetStateIndex) = dot(nosymmout((1+ngoodstates*(targetStateIndex-1)) : ngoodstates*targetStateIndex) &
                                                 , initial_vector(1:ngoodstates))
            stateOverlapValue(targetStateIndex) = abs(stateOverlapValue(targetStateIndex))
            if(stateOverlapValue(targetStateIndex)>=overlapThresh) then
                targetStateFlag = 'getsame'
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
            else
                targetStateFlag = 'ngethigher'
            end if
            if(targetStateIndex >= highestStateIndex) then
                targetStateFlag = 'reachedmax'
                maxOverlapValue = 0.0
                maxOverlapIndex = 0  
                do i=1,highestStateIndex,1
                    if (stateOverlapValue(i) > maxOverlapValue) then
                        maxOverlapValue = stateOverlapValue(i)
                        maxOverlapIndex = i
                    end if
                end do
            end if
        case default
            write(*,*) "unexpected targetStateFlag in subroutine getStateOverlap"
        end select
        write(*,*) "subroutine stateOverlap: targetStateFlag=", targetStateFlag

    end subroutine getStateOverlap
    
    subroutine correctR(singularvalue,leftu,rightv)
        real(kind=8)    ::   singularvalue(subM)
        real(kind=8)    ::   leftu(4*Lrealdim,subM),rightv(subM,4*Rrealdim)
        real(kind=8)    ::   matbuffer(subM,4*Rrealdim)        ! store the result of leftu(T)*coeffIF
        real(kind=8)    ::   rightvBuffer(subM,4*Rrealdim)     ! store the corrected rightv
        real(kind=8)    ::   S(subM,subM),S0(subM),absS(subM,subM)
        integer         ::   i,j,iFound,iGet
        integer         ::   SnonzeroNum,S0nonzeroNum,nonzeroNum,numDisorder
        integer         ::   correctIndex(subM)
        logical         ::   usedIndex(subM)
        logical         ::   isDiagPositive
        
        write(*,*) "Enter in subroutine correctR"
        write(*,*) "S matrix is defined as (U+)*psi*V"
        
        do i=1, subM, 1
            S0(i) = sqrt(singularvalue(i))   !S0 stores the (true) singular value
        end do
        
        ! S stores the singular value matrix calculted by leftu^(T) * coeffIF * rightv
        call gemm(leftu,coeffIF(1:4*Lrealdim,1:4*Rrealdim,formerStateIndex),matbuffer,'T','N',1.0D0,0.0D0)
        call gemm(matbuffer,rightv,S,'N','T',1.0D0,0.0D0) 
        
        if(printSMat == .true.) then 
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

        correctIndex = 0
        numDisorder = 0
        usedIndex = .false.
        ! compare S with S0 to correct rightv
        do i=1, nonzeroNum, 1
            iFound = 0
            do j=1, subM, 1
                if(abs(abs(S(i,j))-S0(i))<singularvalueThresh) then
                    iFound = iFound + 1
                    if(i/=j) then
                        numDisorder = numDisorder + 1
                    end if
                    usedIndex(j) = .true.
                    if(S(i,j)>=0 ) then  
                        correctIndex(i) = j    ! ith row of correct rightv = jth row of uncorrected rightv
                    else    ! find unmatched phase
                        correctIndex(i) = -j
                    end if
                end if
            end do
            if(iFound==0) then
                write(*,*) "iFound=0"
                stop
            end if
        end do
        write(*,*) "Found",numDisorder,"disorder"
        
        rightvBuffer = 0.0        
        iGet = 0
        do i=1, nonzeroNum, 1
            j = correctIndex(i)
            if(j > 0) then
                rightvBuffer(i,:) = rightv(j,:)
            else if(j < 0) then
                j = -j
                rightvBuffer(i,:) = -rightv(j,:)
            end if
        end do
        do i=nonzeroNum+1,subM,1
            do j=1,subM,1
                if(usedIndex(j)==.false.) then
                    rightvBuffer(i,:) = rightv(j,:)
                    usedIndex(j) = .true.
                    iGet = 1
                    exit
                end if
            end do
            if(iGet==0) then
                write(*,*) "iGet==0 error"
                stop
            end if
        end do
        do j=1,subM,1
            if(usedIndex(j)==.false.) then
                write(*,*) "there are unused index"
                stop
            end if
        end do
        
        rightv= rightvBuffer   ! write corrected rightv
        
        isDiagPositive = .true. ! test the subrountine
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
            
        if(printSMat == .true.) then 
            write(*,*) "S matrix"
            do i=1, subM, 1
                write(*,'(32F10.6)') S(i,:)
            end do
        end if
            
        if(isDiagPositive == .true.) then
            write(*,*) "After correction, S matrix become positive and diagonal"
        else
            write(*,*) "After correction, S matrix is not positive and diagonal"
            stop
        end if

    end subroutine correctR
        
end module stateOverlap