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
    nstate , &                 ! how many targeted states
    C2state                    ! if logic_C2/=0, C2state=2*nstate
    real(kind=r8),allocatable :: &    
    nuclQ(:) , &                      ! nuclear Q :: chemical potential :: e.g carbon +1
    coord(:,:) , &                    ! coordinate(natoms,1:3)
    nweight(:) , &                    ! excited states average method nweight
    atommass(:)                       ! atomic mass  
    integer(kind=i4),allocatable :: bondlink(:,:) , & ! bondlink information
                                      atomindex(:)      ! atomindex such as carbon 6
    real(kind=r8) :: cntofmass(3)  ! center of mass used in transition dipole

!========================================================
    
    ! control part
    integer(kind=i4) :: &
    logic_PPP, &            ! if use PPP model    
    logic_meanfield , &     ! if do meanfield scf procedure
    logic_tree, &           ! if use tree tensor algorithm, how many blocks                
    blocks, &
    logic_spinreversal,&    ! if use spinreversal symmetry :: +1 singlet -1 triplet 0 none
    logic_C2          , &   ! if use C2 like symmetry :: +1 A ;-1 B ; 0 none
    logic_C2real     , &    ! the real logic_C2
    logic_bondorder  , &    ! if calculate bond order
    logic_localspin         ! if calculate local spin
    integer(kind=i4),allocatable :: treelink(:,:)  ! treelink information

    character(len=20) :: diagmethod  ,&    ! the diagonalization method
                         PPPpot      ,&    ! PPP potential model
                         opmethod    ,&    ! HC method
                         C2method          ! C2 method mix/independent
!=========================================================
    
    ! loadbalance part
    integer(kind=i4),allocatable ::  &
    orbid1(:,:)    , &        ! the orbid(norbs,2) is the process id every orbital; 
    orbid2(:,:,:)  , &        ! the orbid2(norbs,norbs,2) is the process id every 2 electron operator 
    orbid3(:,:,:)             ! the orbid3(norbs,norbs,2) is the process contain local spin operator

!=========================================================

    ! DMRG part
    integer(kind=i4) :: exscheme         ! target excited state scheme
    integer(kind=i4) :: &
    subM , &                             ! DMRG subspace M 
    subMp, &                             ! the perturbation space M
    sweeps , &                           ! finite DMRG sweeps
    exactsite , &                        ! the number of exact discrible sites
    isweep                               ! the at present isweep
    real(kind=r8) :: energythresh,&        ! energy convegence threshold in the middle of every sweep
                     hopthresh
    real(kind=r8),allocatable :: sweepenergy(:,:)  ! store every sweep energy in the middle
    integer(kind=i4) :: Lrealdim,Rrealdim , &  ! L/R space real dimension
                      Lrealdimp,Rrealdimp    ! L/R space perturbation space dimension
    integer(kind=i4) :: nleft,nright        ! L/R space site number L+sigmaL+sigmaR+R
    integer(kind=i4) :: ngoodstates         ! the number of basis fullfill the good quantum number

    logical :: IfOpenperturbation   ! control the perturbation
    integer :: logic_perturbation   ! if do perturbation in sweep
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
    quantabigR(:,:) , &           ! R space good quantum number (N and Sz)in 4M basis
    quantasmaLp(:,:) , &          ! perturbation mode : L space good quantum number (N and Sz)in M basis
    quantasmaRp(:,:) , &          ! perturbation mode : R space good quantum number (N and Sz)in M basis
    quantabigLp(:,:) , &          ! perturbation mode : L space good quantum number (N and Sz)in 4M basis
    quantabigRp(:,:)              ! perturbation mode : R space good quantum number (N and Sz)in 4M basis

!============================================================

    ! symmetry part
    integer(kind=i2),allocatable :: symmlinksma(:,:,:),symmlinkbig(:,:,:)
    ! symmetrylink represent the symmetry link of every state
    ! the first variable means the states 4M or M
    ! the second means different symmetry spin_reversal/electron-hole
    ! the third is the L/1 and R/2 space

!============================================================

    ! output part
    real(kind=r8),allocatable :: dmrgenergy(:)

!============================================================

    ! constant
    real(kind=r8),parameter :: relazero=1.0D-8   ! relative zero
    real(kind=r8),parameter :: eAtodebye=4.8032038D0,AutoAngstrom=0.5217721092D0

!============================================================

end Module variables
