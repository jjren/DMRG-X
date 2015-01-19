	subroutine infinit_MPS

	USE Variables
	USE MPI
	
	implicit none
	integer :: isystem,i,j,error
	real(kind=8),allocatable :: treal(:,:)
	integer(kind=4),allocatable :: bondlinkreal(:,:)

	if(myid==0) then 
		write(*,*) "enter subroutine infinit_MPS"
	end if
! isweep means the initial finit-MPS stage
	isweep=0
! realnelecs is the real electrons in the system 
	realnelecs=nelecs+ncharges
! when construct the infinit MPS, let the small left space and 
! right space to link together and simulate the real condition

	allocate(treal(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(bondlinkreal(norbs,norbs),stat=error)
	if(error/=0) stop
	treal=t
	bondlinkreal=bondlink

	do isystem=1,norbs/2-1,1
! be careful that the norbs may be odd

	! when doing infinit MPS we use half filled system until arrive the realnelecs
		
		if((isystem+1)*2>realnelecs) then
			nelecs=realnelecs
		else
			nelecs=(isystem+1)*2
		end if
	! in the last step we set the total electron to be real nelectrons
		if(mod(norbs,2)==0 .and. isystem==(norbs/2-1)) then
			nelecs=realnelecs
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
! add the quasi link between the left and right small space
		if(nleft+1<norbs/2 .or. mod(norbs,2)/=0) then 
		t=treal
		bondlink=bondlinkreal
		bondlink(nleft+1,norbs-nright)=1
		bondlink(norbs-nright,nleft+1)=1
		t(nleft+1,norbs-nright)=-2.4D0
		t(norbs-nright,nleft+1)=-2.4D0
		end if


		if(Lrealdim/=Rrealdim .and. myid==0) then
			write(*,*) "---------------------------------------"
			write(*,*) "infinit DMRG Lrealdim/=Rrealdim failed!"
			write(*,*) "---------------------------------------"
			stop
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
	!	if(logic_spinreversal/=0) then
	!		call Spin_reversalmatL
	!		call Spin_reversalmatR
	!	end if
! store the operator matrix and the good quantum number
		call store_operatorL(nleft+1)
		call store_operatorR(norbs-nright)
! direct diagonalization
!		call fullmat

		if(4*Lrealdim>subM) then
! construct the total H(direct method) and davidson diagnalization
		call hamiltonian('i')
		end if
! Renormalization all the operator matrix
		call Renormalization(nleft+1,norbs-nright,'i')
		
!	call MPI_barrier(MPI_COMM_WORLD,ierr)
	end do

	t=treal
	bondlink=bondlinkreal
! when the norbs is odd. Then in the last process of infinit MPS.
! only add 1 orbital in the left space
! and the nelecs set to the the realnelecs
	if(MOD(norbs,2)/=0) then
		nelecs=realnelecs
		nleft=norbs/2
		nright=norbs/2-1
		Lrealdim=subM
		Rrealdim=subM
! caution here may be some problem ? because the right space using the
! last step operamatbig
		call onesitematrix(nleft+1)
		call system_bigL
		call system_constructquantaL
		!if(logical_spinreversal/=0) then
		!	call Spin_reversalmatL
		!end if
		call store_operatorL(nleft+1)
		call hamiltonian('i')
		call Renormalization(nleft+1,norbs-nright,'i')
	end if

	
	deallocate(treal)
	deallocate(bondlinkreal)

	call MPI_barrier(MPI_COMM_WORLD,ierr)

	return
	end subroutine infinit_MPS



