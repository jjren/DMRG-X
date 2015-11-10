subroutine Infinit_MPS
! this subroutine do infinit MPS

	use variables
	use mpi
	use communicate
	use exit_mod
	use Renormalization_mod
	use OnesiteMatrix
	use hamiltonian_mod
	
	implicit none
	! local
	real(kind=r8),allocatable :: treal(:,:)
	integer(kind=i4),allocatable :: bondlinkreal(:,:)
	! treal bondlinkreal store the initial real value
	integer :: error,ierr
	integer :: isystem,i,j
	integer :: logic_C2real
	
	call master_print_message("enter subroutine infinit_MPS")
	
	! in the infinite procedure didnot constrain the logic_C2
	! only in the diagonalization process
	logic_C2real=logic_C2

	! the exactsite refer to the space that L space or R space that can be accurately discribe
	! (without sigmaL and sigmaR)
	exactsite=1
	do while(.true.)
		if(4**exactsite<=subM) then
			exactsite=exactsite+1
		else
			exactsite=exactsite-1
			exit
		end if
	end do

	! if mode==r and isweep/=0 means the infinit DMRG is finished
	! isweep means the initial finit-MPS stage
	if(mode=='r' .and. isweep/=0) then
		nleft=(norbs+1)/2-1
		nright=norbs-nleft-2
		return
	else
		isweep=0
	end if

	! when construct the infinit MPS, let the small left space and 
	! right space to link together and simulate the real condition

	allocate(treal(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(bondlinkreal(norbs,norbs),stat=error)
	if(error/=0) stop
	treal=t
	bondlinkreal=bondlink


	! be careful that the norbs may be odd
	! when doing infinit MPS we use half filled system until arrive the realnelecs
	do isystem=1,norbs/2-1,1
		if(mode=='r') then
			! in the infinite processs if restart the bondorder matrix and
			! local spin matrix will not recover
			! only in the finite process
			if(isystem<nleft) then
				cycle
			else if(isystem==nleft) then
				call Enviro_Big('L')
				call Enviro_Big('R')
				call Hamiltonian('i')
				call Renormalization('i')
				cycle
			end if
		end if

		if((isystem+1)*2>realnelecs) then
			nelecs=realnelecs
		else
			nelecs=(isystem+1)*2
		end if
		call master_print_message(nelecs,"nelecs=")
		
		nleft=isystem
		nright=isystem
		
		if(isystem==1) then
			Lrealdim=1
			Rrealdim=1
		end if
		if(4*Lrealdim>subM) then
			Lrealdim=subM
		else
			Lrealdim=Lrealdim*4
		end if
		if(4*Rrealdim>subM) then
			Rrealdim=subM
		else
			Rrealdim=Rrealdim*4
		end if
		if(Lrealdim/=Rrealdim) then
			call exit_DMRG(sigAbort,"infinit DMRG Lrealdim/=Rrealdim failed!")
		end if

	! add the quasi link between the left and right small space
	! to let the boundry be more real
		t=treal
		bondlink=bondlinkreal
		if(nleft+1<norbs/2 .or. mod(norbs,2)/=0) then 
			bondlink(nleft+1,norbs-nright)=1
			bondlink(norbs-nright,nleft+1)=1
			t(nleft+1,norbs-nright)=-1.0D0
			t(norbs-nright,nleft+1)=-1.0D0
		end if
!====================L space=================================
	! sigmaL subspace operator matrix
		call ConstructOnesiteMatrix(nleft+1)
	! L subspace initial
		if(nleft==1) then
			call Infinit_InitMat('L')
		end if
	! construct the L+sigmaL subspace operator matrix
		call System_Big('L')
	! construct the good quantum number Sz and occpuation
		call System_Constructquanta('L')
	! store the operator matrix and the good quantum number
		call Store_Operator('L')
!============================================================

!===================R space==================================
		if(logic_C2==0) then
			! R subspace initial
			if(nright==1) then
				call Infinit_InitMat('R')
			end if
			! sigmaR subspace operator matrix
			call ConstructOnesiteMatrix(norbs-nright)
			! construct the R+sigmaR subspace operator matrix
			call System_Big('R')
			call System_Constructquanta('R')
		else
			call C2_copy('i')
		end if
		call Store_Operator('R')
!============================================================

	! construct the total H(direct method) and davidson diagnalization
		if(4*Lrealdim>subM) then
			logic_C2=0
			call Hamiltonian('i')
			logic_C2=logic_C2real
		end if
	! Renormalization all the operator matrix
		call Renormalization('i')
	end do

	! set t and bondlink to real value
	t=treal
	bondlink=bondlinkreal

	! when the norbs is odd. Then in the last process of infinit MPS.
	! only add 1 orbital in the left space
	if(MOD(norbs,2)/=0) then
		! add nelecs 2 by 2
		nelecs=nelecs+2
		if(nelecs>realnelecs) then
			nelecs=realnelecs
		end if

		call master_print_message(nelecs,"nelecs=")
		
		nleft=norbs/2
		nright=norbs/2-1
		if(4*Lrealdim>subM) then
			Lrealdim=subM
		else
			Lrealdim=Lrealdim*4
		end if

	! caution here may be some problem ? because the right space using the
	! last step operamatbig quantabigR and so on HbigR
		call ConstructOnesiteMatrix(nleft+1)
		call System_Big('L')
		call System_Constructquanta('L')
		call Store_Operator('L')
		logic_C2=0
		call Hamiltonian('i')
		logic_C2=logic_C2real
		call Renormalization('i')
	end if
	
	deallocate(treal)
	deallocate(bondlinkreal)

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

return
end subroutine Infinit_MPS



