subroutine Hamiltonian(direction)
! this is the core subroutine
! contruct the hamiltonian on the fly and muliply it with coefficient of
! the wavefunction and transfter it to the 0 process

	USE variables
	use communicate
    use stateOverlap

	implicit none

	integer :: &
		lim   ,  &       ! the expanding small space in davidson iteration
		ilow  ,  &       ! index of lowest eigenpair
		ihigh ,  &       ! index of the highest eigenpair
		niv   ,  &       ! number of initial vector
		mblock,  &       ! number of vector to be targeted in each iteration
		maxiter          ! maxiter iterations
	real(kind=8) :: &
		crite ,   &      ! convergence eigenvalue
		critc ,   &      ! convergence coefficiency
		critr ,   &      ! convergence residual vector norm
		ortho            ! orthogonality threshold
	integer,allocatable :: iselec(:)   ! selected eigenpair not used
	integer :: error
	character(len=1) :: direction      
	! i,l,r direction l is L space to R space sweep
	
	call master_print_message("enter hamiltonian subroutine")
    
	! default value
	ilow=1
	ihigh=nstate
	lim=nstate+20    ! this number 20 can be changed consider the efficiency
	niv=nstate
	mblock=nstate
	maxiter=400
    allocate(iselec(lim),stat=error)
	if(error/=0) stop
    iselec=-1
    
    critc=1.0D-9
    critr=1.0D-9
    ortho=1.0D-8
    call getCrite(crite,direction)  !He Ma
    
    if(exscheme == 4 .and. startedMaxOverlap .and. targetStateFlag == 'finished') then
        if(myid==0) then
            write(*,*) "In davidson diagonalization we first only target state", targetStateIndex
        end if
        targetStateFlag = 'trysame'
        iselec(1) = targetStateIndex
        iselec(2) = -1
        ilow = 0
        mblock = 1
        call Davidson_Wrapper(direction,lim,ilow,ihigh,iselec,niv,mblock,crite,critc,critr,ortho,maxiter) ! ilow=0, iselec will be used
        if(targetStateFlag == "finished") then
            deallocate(iselec)
            return
        else    !restore the parameters and going back to davidson(targetting nstate states)
            iselec=-1
            mblock=nstate
            ilow=1
            ihigh=nstate
        end if
    end if

	call Davidson_Wrapper(direction,lim,ilow,ihigh,iselec,niv,mblock,crite,critc,critr,ortho,maxiter)

	deallocate(iselec)
return

end subroutine Hamiltonian


subroutine getCrite(crite,direction)   ! He Ma
!determine the crite(convergence criteria) according to isweep
    use variables
    implicit none
    
    real(kind=8)     ::   crite,originCrite,finalCrite,infiniteCrite
    character(len=1) ::   direction
    integer          ::   pastSweep
    
    infiniteCrite = 1.0D-5
    originCrite = 1.0D-7
    finalCrite  = 1.0D-9
    
!    if(energyThresh >= originCrite) then   ! if energyThresh is too large, no need to refine crite
!        crite = originCrite
!       return
!    end if

    
    if(direction == 'i') then  ! infinite DMRG
        crite = infiniteCrite
    else                       ! finite sweeps
        pastSweep = isweep - 1
        if(exscheme==4 .and. startedMaxOverlap) then
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

end subroutine getCrite


