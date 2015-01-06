	Module Variables
	USE MPI

	implicit none


	! MPI part
	integer(kind=4) :: ierr,myid,nprocs,version,subversion
	


	! system part
	integer(kind=4) :: norbs,nelecs,natoms,ncharges,blocks,realnelecs
	integer(kind=4),allocatable :: orbid(:)
	real(kind=8),allocatable :: nuclQ(:)
	! the orbid(:) is the process id every orbital

	real(kind=8) :: totalSz

	integer :: logic_PPP,nbonds

	integer :: logic_tree
	! if use tree tensor algorithm

	real(kind=8),allocatable :: coord(:,:)
	integer,allocatable :: bondlink(:,:)
	
	integer :: logic_spinreversal
	real(kind=8),allocatable :: adaptedsma(:,:,:),adaptedbig(:,:,:)
	! this is the onesite spin_reversal operator matrix
	real(kind=8) :: parityonesitemat(4,4)
	
	

	! DMRG part
	integer :: subM,sweeps,exactsite
	integer,allocatable :: treelink(:,:)

	! Hamiltonian part

	real(kind=8),allocatable :: t(:,:),v(:,:,:,:) !full Quantum Chemistry
	real(kind=8),allocatable :: hubbardU(:),pppV(:,:) ! PPP model 

	real(kind=8),allocatable :: operamatbig(:,:,:),Hbig(:,:,:)
	real(kind=8),allocatable :: operamatsma(:,:,:),Hsma(:,:,:)
	integer(kind=4) :: operanum
	! operanum is the max site operator every process have
	integer(kind=4),allocatable :: quantasmaL(:,:),quantasmaR(:,:),quantabigL(:,:),quantabigR(:,:)
	! quanta is the good quantum number such as number of electron and number
	! of Sz
	real(kind=8) :: onesitemat(4,4,5)
	integer(kind=4) :: Lrealdim,Rrealdim

	integer(kind=4) :: nleft,nright,ngoodstates
	integer(kind=4) :: nstate,exscheme
	real(kind=8),allocatable :: nweight(:)
	real(kind=8),allocatable :: coeffIF(:,:,:)
	! coeffIF is the inital and final coefficient 

! constant
	real(kind=8),parameter :: quantaconst=0.1
	real(kind=8),parameter :: relazero=1.0D-8



	end Module
