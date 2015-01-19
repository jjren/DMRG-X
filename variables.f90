	Module Variables
	USE MPI

	implicit none


	! MPI part
	integer(kind=4) :: ierr,myid,nprocs,version,subversion

	character(len=1) :: mode
	integer(kind=4) :: modeindex
	


	! system part
	integer(kind=4) :: norbs,nelecs,natoms,ncharges,blocks,realnelecs,totalSz
	integer(kind=4),allocatable :: orbid(:)
	real(kind=8),allocatable :: nuclQ(:)
	! the orbid(:) is the process id every orbital
	! totalSz is the up nelecs - down nelecs

	integer(kind=4) :: logic_PPP,nbonds

	integer(kind=4) :: logic_tree
	! if use tree tensor algorithm

	real(kind=8),allocatable :: coord(:,:)
	integer(kind=4),allocatable :: bondlink(:,:)

! we use a new schema to include symmetry
! logic_spinreversal=+-1
	integer(kind=4) :: logic_spinreversal
! symmetrylink represent the symmetry link of every state
! the first variables means the states
! the second means different symmetry spin_reversal electron-hole
! the third is the L and R space
! symmlinkgood means the good quantum number states symmetry link information
	integer(kind=2),allocatable :: symmlinksma(:,:,:),symmlinkbig(:,:,:),symmlinkgood(:,:)
	!real(kind=8),allocatable :: adaptedsma(:,:,:),adaptedbig(:,:,:)

	! this is the onesite spin_reversal operator matrix
	!real(kind=8) :: parityonesitemat(4,4)
	
	

	! DMRG part
	integer(kind=4) :: subM,sweeps,exactsite
	integer(kind=4),allocatable :: treelink(:,:)

	! Hamiltonian part

	real(kind=8),allocatable :: t(:,:),v(:,:,:,:) !full Quantum Chemistry
	real(kind=8),allocatable :: hubbardU(:),pppV(:,:) ! PPP model 

	real(kind=8),allocatable :: operamatbig(:,:,:),Hbig(:,:,:)
	real(kind=8),allocatable :: operamatsma(:,:,:),Hsma(:,:,:)
	integer(kind=4),allocatable :: quantasmaL(:,:),quantasmaR(:,:),quantabigL(:,:),quantabigR(:,:)
	! quanta is the good quantum number such as number of electron and number
	! of Sz
	real(kind=8) :: onesitemat(4,4,5)
	integer(kind=4) :: Lrealdim,Rrealdim

	integer(kind=4) :: nleft,nright,ngoodstates
	! nleft is the number of orbs in the L space
	! without sigmaL
	integer(kind=4) :: nstate,exscheme
	real(kind=8),allocatable :: nweight(:)
	real(kind=8),allocatable :: coeffIF(:,:,:)
	! coeffIF is the inital and final coefficient 

! constant
	!real(kind=8),parameter :: quantaconst=0.1
	real(kind=8),parameter :: relazero=1.0D-8



	end Module
