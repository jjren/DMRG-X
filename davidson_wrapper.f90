Subroutine davidson_wrapper(direction,LIM,ILOW,IHIGH,ISELEC,NIV,MBLOCK,&
                           CRITE,CRITC,CRITR,ORTHO,MAXITER)
! this is the davidson wrapper to call DVDSON written by Andreas
! aim to allocate memory and set variables in the global array
! mainly allocate memory on 0 process

	USE mpi
	USE variables
	USE InitialGuess

	implicit none

	integer :: error,i,j,k,m
	integer :: N,LIM,ILOW,IHIGH,ISELEC,NIV,MBLOCK,NLOOPS,NMV,ierror,MAXITER
	real(kind=8) :: CRITE,CRITC,CRITR,ORTHO
	logical :: HIEND
	external op
	real(kind=8),allocatable :: HDIAG(:),DavidWORK(:),dummycoeff(:),dummynewcoeff(:)
	integer :: smadim,IWRSZ
	character(len=1) :: direction

! check how many states fullfill good quantum number
! every process do it
! N is the number of good quantum number states
	N=0
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs+ncharges) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
		N=N+1
		end if
	end do
	end do
	ngoodstates=N
!-----------------------------------------------------
	!IWRSZ=2*N*LIM+LIM*LIM+(nstate+10)*LIM+nstate
	IWRSZ=LIM*(2*N+LIM+9)+LIM*(LIM+1)/2+nstate+100

	if(myid==0) then
		write(*,*) "number of good quantum number states",N
		allocate(HDIAG(N),stat=error)
		if(error/=0) stop
		allocate(DavidWORK(IWRSZ),stat=error)
		if(error/=0) stop
	end if

		call GetHDiag(HDIAG)
! get the intital vector
! here we only consider
! 1. nstate=1 finit
! 2. nstate>1 finit exscheme=1
! 3  infinit
! we can add exscheme=2 and logic_spinreversal=1 later
	if(myid==0) then
		if(direction/='i' .and. NIV==1 ) then
			call Initialfinit(DavidWORK)
		else
			call Initialrandomweight(DavidWORK,NIV)
		end if
	end if

	
	if(myid==0) then
		call DVDSON(op,N,LIM,HDIAG,ILOW,IHIGH,ISELEC &
		    ,NIV,MBLOCK,CRITE,CRITC,CRITR,ORTHO,MAXITER,DavidWORK,&
		    IWRSZ,HIEND,NLOOPS,NMV,ierror)
		    smadim=0
		call MPI_bcast(smadim,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	end if

	if(myid/=0) then
		do while(.true.)
			call MPI_bcast(smadim,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
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


	if(myid==0) then
		if(HIEND/=.false.) then
			write(*,*) "---------------------------"
			write(*,*) "didn't get the lowest state"
			write(*,*) "---------------------------"
			stop
		end if
		allocate(coeffIF(4*Lrealdim,4*Rrealdim,IHIGH),stat=error)
		if(error/=0) stop
		coeffIF=0.0D0
! the DavidWORK only contains the ngoodstates coeff.
! other nongoodstates should be set to 0
		m=1
		do k=1,IHIGH,1
		do i=1,4*Rrealdim,1
		do j=1,4*Lrealdim,1
			if((quantabigL(j,1)+quantabigR(i,1)==nelecs+ncharges) .and. &
				quantabigL(j,2)+quantabigR(i,2)==totalSz) then
				coeffIF(j,i,k)=DavidWORK(m)
				m=m+1
			end if
		end do
		end do
		end do

! check if the good quantum number Sz and nelecs if fullfilled
		do k=1,IHIGH,1
		do i=1,4*Rrealdim,1
		do j=1,4*Lrealdim,1
			if((quantabigL(j,1)+quantabigR(i,1)/=nelecs+ncharges) .or. &
				(quantabigL(j,2)+quantabigR(i,2)/=totalSz)) then
				if(abs(coeffIF(j,i,k))>relazero) then
					write(*,*) "------------------------------------"
					write(*,*) "did not fullfill good quantum number",coeffIF(j,i,k)
					write(*,*) "------------------------------------"
				end if
			end if
		end do
		end do
		end do
!--------------------------------------------------------------

		write(*,*) "low state energy"
		do i=1,IHIGH,1
		write(*,*) DavidWORK((i-1)*ngoodstates+1:i*ngoodstates)
		write(*,*) i,"th energy=",DavidWORK(IHIGH*ngoodstates+i)
		end do
		write(*,*) "energy converge:",DavidWORK(IHIGH*ngoodstates+IHIGH+1)
		write(*,*) "residual norm"
		write(*,*) DavidWORK(IHIGH*ngoodstates+IHIGH+2:IHIGH*ngoodstates+2*IHIGH+1)
		write(*,*) "NLOOPS=",NLOOPS
		Write(*,*) "IERROR=",IERROR
		write(*,*) "NMV=",NMV


		deallocate(HDIAG)
		deallocate(DavidWORK)
	end if
return
end subroutine
