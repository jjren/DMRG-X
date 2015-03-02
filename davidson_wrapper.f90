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
	logical :: done

! check how many states fullfill good quantum number
! every process do it
! N is the number of good quantum number states

	if(myid==0) then
		write(*,*) "enter in davidson diagonalization subroutine"
	end if

	N=0
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
		N=N+1
		end if
	end do
	end do
	ngoodstates=N
	
if(myid==0 .and. (logic_spinreversal/=0 .or. logic_C2/=0)) then
	allocate(symmlinkgood(ngoodstates,2),stat=error)
	if(error/=0) stop
! in the good quantum number states space
! get the symmetry link information
! symmlinkgood(m,1) means the left space index in 4M basis
! symmlinkgood(m,2) means the right space index in 4M basis

	m=1
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
		if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
		symmlinkgood(m,1)=j
		symmlinkgood(m,2)=i
		m=m+1
		end if
	end do
	end do
		if(m-1/=ngoodstates) then
			write(*,*) "-----------------------------"
			write(*,*) "symmlinkgood m-1/=ngoodstates"
			write(*,*) "-----------------------------"
			stop
		end if
end if

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
! we can add exscheme=2 and later
	if(myid==0) then
		if(direction/='i' .and. NIV==1 .and. logic_C2==0) then
			call Initialfinit(DavidWORK,direction)
		else
		!	call Initialunivector(HDIAG,DavidWORK,NIV)
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
		
		if(logic_spinreversal/=0) then
			call spincorrect(DavidWORK(1:ngoodstates*nstate))
		end if
		
		coeffIF=0.0D0
! the DavidWORK only contains the ngoodstates coeff.
! other nongoodstates should be set to 0
		m=1
		do k=1,IHIGH,1
		do i=1,4*Rrealdim,1
		do j=1,4*Lrealdim,1
			if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
				quantabigL(j,2)+quantabigR(i,2)==totalSz) then
				coeffIF(j,i,k)=DavidWORK(m)
				m=m+1
			end if
		end do
		end do
		end do

! check if the good quantum number Sz and nelecs if fullfilled
!		do k=1,IHIGH,1
!		do i=1,4*Rrealdim,1
!		do j=1,4*Lrealdim,1
!			if((quantabigL(j,1)+quantabigR(i,1)/=nelecs) .or. &
!				(quantabigL(j,2)+quantabigR(i,2)/=totalSz)) then
!				if(abs(coeffIF(j,i,k))>relazero) then
!					write(*,*) "------------------------------------"
!					write(*,*) "did not fullfill good quantum number",coeffIF(j,i,k)
!					write(*,*) "------------------------------------"
!				end if
!			end if
!		end do
!		end do
!		end do
!--------------------------------------------------------------

		write(*,*) "low state energy"
		do i=1,IHIGH,1
!		write(*,*) DavidWORK((i-1)*ngoodstates+1:i*ngoodstates)
		write(*,*) nleft+1,norbs-nright,i,"th energy=",DavidWORK(IHIGH*ngoodstates+i)
		write(*,*) "energy converge:",DavidWORK(IHIGH*ngoodstates+IHIGH+i)
		write(*,*) "residual norm:",DavidWORK(IHIGH*ngoodstates+2*IHIGH+IHIGH+i)
		end do
		write(*,*) "NLOOPS=",NLOOPS
		Write(*,*) "IERROR=",IERROR
		write(*,*) "NMV=",NMV

		if(IERROR/=0) then
			write(*,*) "---------------------------"
			write(*,*) "failed! IERROR=",IERROR
			write(*,*) "---------------------------"
			stop
		end if
		
! update the sweepenergy
		do i=1,IHIGH,1
			if(DavidWORK(IHIGH*ngoodstates+i)<sweepenergy(isweep,i)) then
				sweepenergy(isweep,i)=DavidWORK(IHIGH*ngoodstates+i)
			end if
		end do

		
		deallocate(HDIAG)
		deallocate(DavidWORK)
		if(logic_spinreversal/=0 .or. logic_C2/=0) then
			deallocate(symmlinkgood)
		end if
	end if
return
end subroutine
