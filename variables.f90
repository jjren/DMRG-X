module variables
	
	use kinds_mod
	implicit none

!=================================================
	character(len=1) :: mode
	integer(kind=i4) :: modeindex
	! mode=s standard 
	! mode=d debug
	! mode=r restart
!=================================================

	! molecule part
	integer(kind=i4) :: &
	norbs, &                   ! how many orbitals
	nelecs, &                  ! how many electrons
	formernelecs, &            ! the last step electron number
	natoms, &                  ! how many atoms
	ncharges, &                ! how many extra charges (+1,add one electron)
	realnelecs, &              ! ncharges+nelecs
	totalSz, &                 ! totalSz component up nelecs - down nelecs
	nbonds , &                 ! how many bonds
	nstate                     ! how many targeted states
	real(kind=r8),allocatable :: &    
	nuclQ(:) , &                      ! nuclear Q :: chemical potential :: e.g carbon +1
	coord(:,:) , &                    ! coordinate(natoms,1:3)
	nweight(:)                        ! excited states average method nweight
	integer(kind=i4),allocatable :: bondlink(:,:)  ! bondlink information

!========================================================
	
	! control part
	integer(kind=i4) :: &
	logic_PPP, &            ! if use PPP model    
	logic_meanfield , &     ! if do meanfield scf procedure
	logic_tree, &           ! if use tree tensor algorithm, how many blocks                
	blocks, &
	logic_spinreversal,&    ! if use spinreversal symmetry :: +1 singlet -1 triplet 0 none
	logic_C2                ! if use C2 like symmetry :: +1 A ;-1 B ; 0 none
	integer(kind=i4),allocatable :: treelink(:,:)  ! treelink information
    
	logical  :: logic_fullmat = .false.  ! whether do direct diagonalization, default not

	integer(kind=i4),parameter :: &
	logic_bondorder= 0       ! if calculate bond order
!=========================================================
	
	! loadbalance part
	integer(kind=i4),allocatable ::  &
	orbid1(:,:) , &        ! the orbid(norbs,2) is the process id every orbital; 
	orbid2(:,:,:)         ! the orbid2(norbs,norbs,2) is the process id every 2 electron operator 

!=========================================================

	! DMRG part
	integer(kind=i4) :: exscheme         ! target excited state scheme
	integer(kind=i4) :: &
	subM , &                             ! DMRG subspace M 
	maxSweeps, &                         ! maximum finite DMRG sweeps number
	sweeps , &                           ! how many finit sweeps have actually been done
	exactsite , &                        ! the number of exact discrible sites
	stepPerSweep, &                      ! how many steps in one finite DMRG sweep
	isweep                               ! the at present isweep
	real(kind=r8) :: energythresh        ! energy convegence threshold in the middle of every sweep
	logical       :: reachedEnergyThresh   = .false.  !during finite sweep, whether the crite of davidson reached 0.1*energythresh
	real(kind=r8),allocatable :: sweepenergy(:,:)     ! store every sweep energy in the middle
	integer(kind=i4) :: Lrealdim,Rrealdim   ! L/R space real dimension
	integer(kind=i4) :: nleft,nright        ! L/R space site number L+sigmaL+sigmaR+R
	integer(kind=i4) :: ngoodstates         ! the number of basis fullfill the good quantum number

!=========================================================

	! Hamiltonian part
	real(kind=r8),allocatable :: & 
	t(:,:) , &                    ! transfer integral in PPP/one electron integral  
	v(:,:,:,:) , &                ! two electron integral full Quantum Chemistry
	hubbardU(:) , &               ! hubbard term
	pppV(:,:)                     ! PPP term
	integer(kind=i4),allocatable :: &
	quantasmaL(:,:) , &           ! L space good quantum number (N and Sz)in M basis
	quantasmaR(:,:) , &           ! R space good quantum number (N and Sz)in M basis
	quantabigL(:,:) , &           ! L space good quantum number (N and Sz)in 4M basis
	quantabigR(:,:)               ! R space good quantum number (N and Sz)in 4M basis
!============================================================
    ! state overlap information 
	logical          :: startedStateSpecific = .false. ! whether conduct state-specific DMRG by maximum overlap algorithm
	integer(kind=4)  :: formerStateIndex               ! the excited state of last step
	integer(kind=4)  :: targetStateIndex               ! the excited state to be targetted (traced) in this step
	character(len=10):: targetStateFlag  = 'none'           ! 'none': not started yet or already finished
	integer          :: realTargetStateIndex           ! user specified index which to be traced
	integer(kind=4)  :: highestStateIndex              ! highest state considered when tracing excited state
	real(kind=r8)     :: overlapThresh                  ! If two states have a bigger overlap than this value,
                                                            ! 'trysame': target the state with the same index of last step
                                                            ! 'getsame': the state above is correct(overlap is large)
                                                            ! 'ngetsame': the state above is incorrect(overlap is small)
                                                            ! 'trylower': target the states with lower indices
                                                            ! 'getlower': found correct state
                                                            ! 'ngetlower': didn't found correct state
                                                            ! 'tryhigher': target a state with higher index
                                                            ! 'gethigher': the state above is correct
                                                            ! 'ngethigher': the state above is incorrect
                                                            ! 'reachedmax': all states below highestStateIndex are incorrect
!============================================================

	! symmetry part
	integer(kind=i2),allocatable :: symmlinksma(:,:,:),symmlinkbig(:,:,:),symmlinkgood(:,:)
	! symmetrylink represent the symmetry link of every state
	! the first variable means the states 4M or M
	! the second means different symmetry spin_reversal/electron-hole
	! the third is the L/1 and R/2 space
	! symmlinkgood means the good quantum number states symmetry link information

!============================================================

	! constant
	real(kind=r8),parameter :: relazero=1.0D-8   ! relative zero

!============================================================

end Module variables
