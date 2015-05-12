module max_overlap
    use variables
    use Symmetry
    use InitialGuess
    use blas95
    use lapack95
    
    contains
    subroutine getMaxOverlapStateIndex(Davidwork,dimN, IHIGH, direction)
        integer      ::   IHIGH, dimN
        character(len=1)  ::   direction
        real(kind=r8)     ::   Davidwork(dimN*IHIGH)
        real(kind=r8)     ::   initial_vector(dimN)
        real(kind=r8)::   working_overlap, max_overlap
        integer      ::   working_index, max_overlap_index
        integer      ::   m, n, i, j, k, error
        
        write(*,*) "Enter in subroutine getMaxOverlapStateIndex"
        write(*,*) "former targetted index = ", targettedStateIndex
        
        working_overlap = 0.0
        max_overlap = 0.0
        
        
        call InitialStarter(direction,dimN,1,initial_vector)
        
        write(*,*) "direction = ", direction
        write(*,*) "startedMaxOverlap", startedMaxOverlap
        do working_index=1,IHIGH,1
            working_overlap = dot(Davidwork((1+dimN*(working_index-1)) : dimN*working_index) &
                              , initial_vector(1:dimN))
            working_overlap = abs(working_overlap)
            write(*,*) "overlap between davidson solution", working_index, "and initial vector is"
            write(*,*) working_overlap
            if (working_overlap > max_overlap) then
                max_overlap = working_overlap
                max_overlap_index = working_index
            end if
        end do      
        
        if(targettedStateIndex .NE. max_overlap_index) then
            write(*,*) "targettedStateIndex .NE. max_overlap_index ", targettedStateIndex, max_overlap_index
        end if
        
        targettedStateIndex =  max_overlap_index 
        
        write(*,*) "new targettedStateIndex = ", targettedStateIndex

    end subroutine getMaxOverlapStateIndex
    
    subroutine correctR(singularvalue,leftu,rightv)
        real(kind=8)    ::   singularvalue(subM)
        real(kind=8)    ::   leftu(4*Lrealdim,subM),rightv(subM,4*Rrealdim)
        real(kind=8)    ::   matbuffer(subM,4*Rrealdim)        ! store the result of leftu(T)*coeffIF
        real(kind=8)    ::   rightvBuffer(subM,4*Rrealdim)     ! store the corrected rightv
        real(kind=8)    ::   S(subM,subM),S0(subM),absS(subM,subM)
        integer         ::   i,j,iFound
        integer         ::   SnonzeroNum, S0nonzeroNum
        integer         ::   correctIndex(subM)
        
        write(*,*) "Enter in subroutine correctR"
        
        do i=1, subM, 1
            S0(i) = sqrt(singularvalue(i))   !S0 stores the (true) singular value
        end do

        ! S stores the singular value matrix calculted by leftu^(T) * coeffIF * rightv
        call gemm(leftu,coeffIF(1:4*Lrealdim,1:4*Rrealdim,targettedStateIndex),matbuffer,'T','N',1.0D0,0.0D0)
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
        ! compare S with S0 to correct rightv
        do i=1, subM, 1
            iFound = 0
            do j=1, subM, 1
                if(abs(abs(S(i,j))-S0(i))<singularvalueThresh) then
                    iFound = iFound + 1
                    if(i/=j) then
                        write(*,*) "Found disorder"
                    end if
                    if(S(i,j)>=0 ) then  
                        correctIndex(i) = j    ! ith row of correct rightv = jth row of uncorrected rightv
                    else    ! find unmatched phase
                        correctIndex(i) = -j
                    end if
                end if
            end do
            if(iFound/=1) then
                write(*,*) "iFound/=1"
                stop
            end if
        end do
        
        rightvBuffer = 0.0        
        do i=1, subM, 1
            j = correctIndex(i)
            if(j > 0) then
                rightvBuffer(i,:) = rightv(j,:)
            else if(j < 0) then
                j = -j
                rightvBuffer(i,:) = -rightv(j,:)
            else
                write(*,*) "correctIndex = 0 error"
                stop
            end if
        end do
        
        rightv= rightvBuffer
        
        if(printSMat == .true.) then ! test the subrountine
            call gemm(matbuffer,rightv,S,'N','T',1.0D0,0.0D0) 
            write(*,*) "After correction, S matrix"
            do i=1, subM, 1
                write(*,'(32F10.6)') S(i,:)
            end do
        end if
        
    end subroutine correctR
        
end module