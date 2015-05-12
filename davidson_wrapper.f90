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
    use maxOverlap

	implicit none

	character(len=1) :: direction
	integer :: lim,ilow,ihigh,niv,mblock,maxiter
	integer :: iselec(lim)
	real(kind=8) :: crite,critc,critr,ortho

	! local
	! davidson parameter
	integer :: dimN,nloops,nmv,ierror,smadim,IWRSZ,NUME
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
    
    if(exscheme == 4 .and. startedMaxOverlap .and. targetStateFlag == 'trysame') then
        NUME = targetStateIndex
    else
        NUME = nstate
    end if
    niv = NUME

	ngoodstates=0
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
			ngoodstates=ngoodstates+1
		end if
	end do
	end do
	dimN=ngoodstates
	
	if((logic_spinreversal/=0 .or. &
		(logic_C2/=0 .and. nleft==nright))) then
		! construct the symmlinkgood 
		call CreatSymmlinkgood
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
		allocate(nosymmout(ngoodstates*NUME),stat=error)
		if(error/=0) stop

		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
			do i=1,NUME
				call SymmetrizeState(ngoodstates,&
					nosymmout((i-1)*ngoodstates+1:i*ngoodstates),Davidwork((i-1)*nsymmstate+1:i*nsymmstate),'u')
			end do
		else
			nosymmout=DavidWORK(1:ngoodstates*NUME)
		end if
		!  no need do GramSchmit; the result fullfill orthnormal
		!  call GramSchmit(NUME,ngoodstates,nosymmout,norm)
		!  write(*,*) "davidson nosymmout norm=",norm
		 
!=================================================================================
! He Ma
        ! To determine in the davidson solutions which state has the maximum overlap with the previous-step excited state. 
        ! write global variable targettedStateIndex
        if(exscheme == 4 .and. startedMaxOverlap) then 
            call getStateOverlap(Davidwork, dimN, ihigh, direction)
        end if
!=================================================================================
        
		coeffIF=0.0D0
! the DavidWORK only contains the ngoodstates coeff other nongoodstates should be set to 0
        m=1
        do k=1,NUME,1
        do i=1,4*Rrealdim,1
        do j=1,4*Lrealdim,1
            if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
                quantabigL(j,2)+quantabigR(i,2)==totalSz) then
                coeffIF(j,i,k)=nosymmout(m)
                m=m+1
            end if
        end do
        end do
        end do
! write the coeffIF
		reclength=nstate*32*subM*subM
		open(unit=109,file="coeffIF.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
		write(109,rec=1) coeffIF
		close(109)

!=================================================================================
! write the final out

        if(exscheme == 4 .and. startedMaxOverlap .and. targetStateFlag == 'trysame') then
            write(*,*) "Energy of the", targetStateIndex, "state is", &
                        DavidWORK(NUME*dimN+targetStateIndex)
            write(*,*) "energy converge:",DavidWORK(NUME*dimN+NUME+targetStateIndex)
            write(*,*) "residual norm:",DavidWORK(NUME*dimN+2*NUME+targetStateIndex)
        else
		    write(*,*) "Energy of lowest", NUME, "states"
		    do i=1,NUME,1
			    write(*,*) nleft+1,norbs-nright,i,"th energy=",DavidWORK(NUME*dimN+i)
			    write(*,*) "energy converge:",DavidWORK(NUME*dimN+NUME+i)
			    write(*,*) "residual norm:",DavidWORK(NUME*dimN+2*NUME+i)
            end do
        end if
        
		write(*,*) "NLOOPS=",nloops
		Write(*,*) "IERROR=",ierror
		write(*,*) "NMV=",nmv
        write(*,*) "crite=",crite

		if(ierror/=0) then
			call master_print_message("failed! IERROR/=0")
			stop
		end if

!==================================================================================
		
! update the sweepenergy
! use the middle site as the sweepenergy

        if(exscheme == 4 .and. startedMaxOverlap) then
            sweepenergy(isweep,:) = DavidWORK(NUME*dimN+targetStateIndex)
        else
		    if(nleft==(norbs+1)/2-1) then
			    do i=1,NUME,1
				    sweepenergy(isweep,i)=DavidWORK(NUME*dimN+i)
			    end do
            end if
        end if
		
		deallocate(HDIAG)
		deallocate(DavidWORK)
		deallocate(nosymmout)
	end if

! deallocate symmetry workarray
	if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
		call DestorySymm
	end if

return
end subroutine Davidson_Wrapper
