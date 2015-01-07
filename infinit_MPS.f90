	subroutine infinit_MPS

	USE Variables
	USE MPI
	
	integer :: isystem,i,j

	if(myid==0) then 
		write(*,*) "enter subroutine infinit_MPS"
	end if

! realnelecs is the real electrons in the system
	realnelecs=nelecs+ncharges

	do isystem=1,norbs/2-1,1

	! when doing infinit MPS we use half filled system until arrive the realnelecs
		if((isystem+1)*2>realnelecs) then
			nelecs=realnelecs
		else
			nelecs=(isystem+1)*2
		end if
	!--------------------------
		nleft=isystem
		nright=isystem
		if(4**nleft<subM) then
			Lrealdim=4**nleft
		else
			Lrealdim=subM
		end if
		if(4**nright<subM) then
			Rrealdim=4**nleft
		else
			Rrealdim=subM
		end if
		
		if(Lrealdim/=Rrealdim .and. myid==0) then
			write(*,*) "---------------------------------------"
			write(*,*) "infinit DMRG Lrealdim/=Rrealdim failed!"
			write(*,*) "---------------------------------------"
		end if
! sigmaL subspace operator matrix
		call onesitematrix(nleft+1)
! L subspace initial
		if(nleft==1) then
			call infinit_smallL
		end if
! construct the L+sigmaL subspace operator matrix
		call system_bigL
! sigmaR subspace operator matrix
		call onesitematrix(norbs-nright)
! R subspace initial
		if(nright==1) then
			call infinit_smallR
		end if
! construct the R+sigmaR subspace operator matrix
		call system_bigR
! construct the good quantum number Sz and occpuation
		call system_constructquantaL
		call system_constructquantaR
! construct the spin_reversal adapted matrix
		if(logic_spinreversal/=0) then
			call Spin_reversalmatL
			call Spin_reversalmatR
		end if
! store the operator matrix and the good quantum number
		call store_operator(nleft+1,norbs-nright)
		if(4*Lrealdim>subM) then
! construct the total H(direct method) and davidson diagnalization
		call hamiltonian('i')
		end if
! Renormalization all the operator matrix
		call Renormalization(nleft+1,norbs-nright,'i')
		
!	call MPI_barrier(MPI_COMM_WORLD,ierr)
	end do

! when the norbs is odd. Then in the last process of infinit MPS.
! only add 1 orbital in the left space
	if(MOD(norbs,2)/=0) then


	return
	end subroutine infinit_MPS



