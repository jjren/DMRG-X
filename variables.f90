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

!=========================================================
	
	! loadbalance part
	integer(kind=i4),allocatable :: orbid(:)         ! the orbid(norbs) is the process id every orbital

!=========================================================

	! DMRG part
	integer(kind=i4) :: exscheme         ! target excited state scheme
	integer(kind=i4) :: &
	subM , &                             ! DMRG subspace M 
	sweeps , &                           ! finite DMRG sweeps
	exactsite , &                        ! the number of exact discrible sites
	isweep                               ! the at present isweep
	real(kind=r8) :: energythresh        ! energy convegence threshold in the middle of every sweep
    logical       :: reachedEnergyThresh   = .false. !during finite sweep, whether the crite of davidson reached 0.1*energythresh
	real(kind=r8),allocatable :: sweepenergy(:,:)  ! store every sweep energy in the middle
	integer(kind=i4) :: Lrealdim,Rrealdim   ! L/R space real dimension
	integer(kind=i4) :: nleft,nright        ! L/R space site number L+sigmaL+sigmaR+R
	integer(kind=i4) :: ngoodstates         ! the number of basis fullfill the good quantum number

!=========================================================

	! Hamiltonian part
	real(kind=r8),allocatable :: & 
	t(:,:) , &                    ! transfer integral in PPP/one electron integral  
	v(:,:,:,:) , &                ! two electron integral full Quantum Chemistry
	hubbardU(:) , &               ! hubbard term
	pppV(:,:) , &                 ! PPP term
	operamatbig(:,:,:) , &        ! operator matrix in 4M basis a+(up),a+(down),n  
	Hbig(:,:,:) , &               ! subspace HL and HR in 4M basis
	operamatsma(:,:,:) , &        ! operator matrix in M basis a+(up),a+(down),n
	Hsma(:,:,:)                   ! subspace HL and HR in 4M basis
	integer(kind=i4),allocatable :: &
	quantasmaL(:,:) , &           ! L space good quantum number (N and Sz)in M basis
	quantasmaR(:,:) , &           ! R space good quantum number (N and Sz)in M basis
	quantabigL(:,:) , &           ! L space good quantum number (N and Sz)in 4M basis
	quantabigR(:,:)               ! R space good quantum number (N and Sz)in 4M basis
	real(kind=r8) :: onesitemat(4,4,5)            ! one site matrix in 4*4 basis 
	real(kind=r8),allocatable :: coeffIF(:,:,:)   ! coeffIF is the inital and final wavefunction coefficient 

!============================================================
    !He Ma  max overlap
    logical         :: startedMaxOverlap     = .false.     !whether conduct state-specific DMRG by maximum overlap
    integer(kind=4) :: targetStateIndex      = 1           !the excited state to be targetted (traced)
    
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
