module maxOverlap
    use variables
    use Symmetry
    use InitialGuess
    use blas95
    use lapack95
    
    contains
    subroutine getStateOverlap(Davidwork,dimN, IHIGH, direction)
        integer      ::   IHIGH, dimN
        character(len=1)  ::   direction
        real(kind=r8)     ::   Davidwork(dimN*IHIGH)
        real(kind=r8)     ::   initial_vector(dimN)
        real(kind=r8)::   workingOverlap, maxOverlap
        integer      ::   workingIndex, maxOverlapIndex
        integer      ::   m, n, i, j, k, error
        
        write(*,*) "Enter in subroutine getStateOverlap"
        write(*,*) "direction = ", direction
        write(*,*) "former target state index = ", targetStateIndex
        
        workingOverlap = 0.0
        maxOverlap = 0.0
        
        call InitialStarter(direction,dimN,1,initial_vector)
        
        select case(targetStateFlag)
        case('uncertain')
            do workingIndex=1,IHIGH,1
                workingOverlap = dot(Davidwork((1+dimN*(workingIndex-1)) : dimN*workingIndex) &
                                  , initial_vector(1:dimN))
                workingOverlap = abs(workingOverlap)
                write(*,*) "overlap between davidson solution", workingIndex, "and initial vector is"
                write(*,*) workingOverlap
                if (workingOverlap > maxOverlap) then
                    maxOverlap = workingOverlap
                    maxOverlapIndex = workingIndex
                end if
            end do      
            if(maxOverlap < 0.9) then
                write(*,*) "Caution! Max overlap < 0.9"
            end if
            if(targetStateIndex .NE. maxOverlapIndex) then
                write(*,*) "targetStateIndex .NE. maxOverlapIndex "
            end if
            targetStateIndex =  maxOverlapIndex 
            write(*,*) "new targetStateIndex = ", targetStateIndex
            targetStateFlag = 'finished'
        case('trysame')
            workingOverlap = dot(Davidwork((1+dimN*(targetStateIndex-1)) : dimN*targetStateIndex) &
                                  , initial_vector(1:dimN))
            workingOverlap = abs(workingOverlap)
            if(workingOverlap>=0.9) then
                write(*,*) "target state index keeps the same, with overlap =", workingOverlap
                targetStateFlag = 'getsame'
            else 
                write(*,*) "target state index may changed, with overlap=", workingOverlap
                write(*,*) "going back to davidon diagonalization"
                targetStateFlag = 'uncertain'
            end if
        case default
            write(*,*) "targetStateFlag case default error"
            stop
        end select

    end subroutine getStateOverlap
    
    subroutine correctR(singularvalue,leftu,rightv)
        real(kind=8)    ::   singularvalue(subM)
        real(kind=8)    ::   leftu(4*Lrealdim,subM),rightv(subM,4*Rrealdim)
        real(kind=8)    ::   matbuffer(subM,4*Rrealdim)        ! store the result of leftu(T)*coeffIF
        real(kind=8)    ::   rightvBuffer(subM,4*Rrealdim)     ! store the corrected rightv
        real(kind=8)    ::   S(subM,subM),S0(subM),absS(subM,subM)
        integer         ::   i,j,iFound
        integer         ::   SnonzeroNum,S0nonzeroNum,numDisorder
        integer         ::   correctIndex(subM)
        logical         ::   isDiagPositive
        
        write(*,*) "Enter in subroutine correctR"
        write(*,*) "S matrix is defined as (U+)*psi*V"
        
        do i=1, subM, 1
            S0(i) = sqrt(singularvalue(i))   !S0 stores the (true) singular value
        end do

        ! S stores the singular value matrix calculted by leftu^(T) * coeffIF * rightv
        call gemm(leftu,coeffIF(1:4*Lrealdim,1:4*Rrealdim,targetStateIndex),matbuffer,'T','N',1.0D0,0.0D0)
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
            if(S0(i)>=singularvalueThreshA) then
                S0nonzeroNum = S0nonzeroNum + 1
            end if
        end do
        
        SnonzeroNum = 0
        do i=1, subM, 1
           do j=1, subM, 1
               if(abs(S(i,j))>singularvalueThresh) then
                   SnonzeroNum = SnonzeroNum + 1
               end if
           end do
        end do        
        
        write(*,*) "S0 has ", S0nonzeroNum, "nonzero elements"
        write(*,*) "S matrix has ", SnonzeroNum, "nonzero elements"
        
        if(S0nonzeroNum/=SnonzeroNum) then
            write(*,*) "the difference may be due to rightv row with small singular value"
            write(*,*) "don't have corresponding leftu column (discarded)"
        end if
        
        correctIndex = 0
        numDisorder = 0
        ! compare S with S0 to correct rightv
        do i=1, subM, 1
            iFound = 0
            do j=1, subM, 1
                if(abs(abs(S(i,j))-S0(i))<singularvalueThresh) then
                    iFound = iFound + 1
                    if(i/=j) then
                        numDisorder = numDisorder + 1
                    end if
                    if(S(i,j)>=0 ) then  
                        correctIndex(i) = j    ! ith row of correct rightv = jth row of uncorrected rightv
                    else    ! find unmatched phase
                        correctIndex(i) = -j
                    end if
                end if
            end do
            if(iFound>1) then
                write(*,*) "iFound>1"
                stop
            end if
        end do
        write(*,*) "Found",numDisorder,"disorder"
        
        rightvBuffer = 0.0        
        do i=1, subM, 1
            j = correctIndex(i)
            if(j > 0) then
                rightvBuffer(i,:) = rightv(j,:)
            else if(j < 0) then
                j = -j
                rightvBuffer(i,:) = -rightv(j,:)
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
        
end module