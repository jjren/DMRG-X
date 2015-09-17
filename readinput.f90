Subroutine ReadInput
! this subroutine is to read in the input parameters
! process 0 read the input file and distribute to other process

	USE variables
	USE PPP_term
	USE MPI
	use communicate
	implicit none

	integer :: i,j,error,ierr 
	integer :: link1,link2           ! link1 and link2 is the bondlink atom information
	integer :: junk_natoms           ! to compare if the natoms is right
	character(len=2) :: symbol       ! the element symbol not used

	character(len=1),allocatable :: packbuf(:)  
	integer :: position1
	integer :: packsize       
	real(kind=r8) :: dummyt           ! transfer integral intermediate variable

	if(myid==0) then
	call master_print_message("enter in readinput subroutine")
!===============================================================
	! default value
	exscheme=0
	! exscheme is the calculate excited states scheme
	! if only ground state is need,exscheme==0
	blocks=2
	! blocks is the total number of blocks
	! when blocks/=2 tree tensor
	modeindex=0
	! modeindex==0 is the standard mode
	mode='s'
	! mode=s means standard mode start from scratch
	! mode=d means debug mode
	! mode=r means restart mode
!==============================================================

	open(unit= 10,file="inp",status="old")
	read(10,*) mode,modeindex ! which mode you want to use
	read(10,*) nthreads(1:4)
! including standard/restart/debug/
! in the restart mode we need the isweep,nleft,nright
! modeindex=1 the fromleftsweep 
! modeindex=-1 the from right sweep
! isweep is the now sweep 
! nleft and nright is the L space and R space that the operator matrix is in the disc
	
	if(mode=='r') then
		read(10,*) isweep,nleft,nright
	end if

	read(10,*) norbs                 ! how many orbitals
	read(10,*) junk_natoms           ! how many atoms
	read(10,*) nelecs                ! how many electrons
	read(10,*) ncharges              ! how many extra charges +1 means add 1 electron
	read(10,*) totalsz               ! total Sz of the system
	read(10,*) logic_PPP             ! if do PPP model logic_PPP=1
	read(10,*) logic_MeanField       ! if do meanfield SCF calculation
	read(10,*) logic_spinreversal    ! if do spin reversal logic_spinreversal=+-1 
	read(10,*) logic_C2              ! if do C2 symmetry or the same mirror reflection and center reflection
	read(10,*) logic_tree            ! if do tree tensor algorithm logic_tree=1
	read(10,*) subM                  ! DMRG Sub space  M
	read(10,*) sweeps                ! DMRG how many sweeps
	read(10,*) nstate                ! how many state wanted to get
	read(10,*) energythresh          ! the threshold of the total energy you want to get
	read(10,*) diagmethod            ! the diagonalization method

! sweepenergy is the total energy of every sweep(in the middle of the chain)
	allocate(sweepenergy(0:sweeps,nstate),stat=error)
	if(error/=0) stop
	sweepenergy=0.0D0
	allocate(dmrgenergy(nstate))
! 
	allocate(nweight(nstate),stat=error)
	if(error/=0) stop
	nweight=1.0D0         !default value
	
	if(nstate/=1) then
		read(10,*) exscheme ! exschemem=1 average method =2 my new specifc method =3 the dipole operator average method
		if(exscheme==1 .or. exscheme==3) then
			read(10,*) nweight(1:nstate) ! nweight is the average DMRG excited state
			nweight=nweight/sum(nweight)
		end if
	end if

	if(logic_tree==1) then
		read(10,*) blocks ! DMRG how many blocks . Tree tensor network
	end if

	if (logic_spinreversal/=0 .and. totalsz/=0) then
		call master_print_message("spin reversal needs Sz=0, failed!")
		stop
	end if

!-------------------------------------------------------------
	! if logic_ppp=1 read how many bonds and link information  
	if(logic_PPP==1) then
		read(10,*) nbonds
		allocate(bondlink(norbs,norbs),stat=error)
		if(error/=0) stop

	! coordiantes of the system
		open(unit= 12,file="coord.xyz",status="old")
		read(12,*) natoms
	  
	! if logic_PPP=1 read nuclear q (physically, chemical potential)
		allocate(nuclQ(natoms),stat=error)
		if(error/=0) stop
	  
		if(natoms/=junk_natoms) then
			call master_print_message("!the natoms has some problem natoms/=junk_natoms!")
			stop
		end if
  
		allocate(coord(3,1:natoms),stat=error)
		if(error/=0) stop
		allocate(atomindex(natoms),stat=error)
		if(error/=0) stop
		allocate(atommass(natoms),stat=error)
		if(error/=0) stop

		read(12,*)
		do i=1,natoms,1   
			read(12,*) symbol,coord(1:3,i),nuclQ(i)
			call symboltoatomindex(symbol,atomindex(i),atommass(i))
		end do
		! get the center of mass
		call centerofmass
		close(12)
	end if

! integral of the system
	open(unit=14,file="integral.inp",status="old")
	if(logic_PPP==1) then
		if(natoms/=norbs) then
			call master_print_message("!Use PPP model, the norbs/=natoms failed!")
			stop
		end if
		allocate(hubbardU(norbs),stat=error)
		if(error/=0) stop
		allocate(t(norbs,norbs),stat=error)
		if(error/=0) stop
	! be careful about the structure of the integral formatted
		t=0.0D0
		bondlink=0
		do i=1,nbonds,1
			read(14,*) link1,link2,dummyt
			bondlink(link1,link2)=1  ! if linked , bondlink=1
			bondlink(link2,link1)=1
			t(link1,link2)=dummyt
			t(link2,link1)=dummyt
		end do

		do i=1,norbs,1
			bondlink(i,i)=2    ! bondlink=2 means that it is the site energy
			read(14,*) t(i,i)  ! t(i,i) is the site energy
		end do
		
		hubbardU=0.0D0
		do i=1,norbs,1
			read(14,*) hubbardU(i)  ! read in every hubbard U
		end do
		
	else
		write(*,*) "--------------------------------------------------"
		write(*,*) "in the QC-DMRG case readin the FCIDUMP integrals"
		write(*,*) "--------------------------------------------------"
	end if  
	
	if(logic_tree==1) then
		open(unit=16,file="tree.inp",status="old")
		allocate(treelink(blocks+1,norbs),stat=error)
		if(error/=0) stop
		do i=1,norbs,1
		read(16,*) treelink(1:blocks+1,i)
		end do
		close(16)
	end if
	
	close(14)
	close(10)
	end if
!===================================================================================
! broadcast to other process

	packsize=100000
	allocate(packbuf(packsize),stat=error)
	if(error/=0) stop

	if(myid==0) then
		position1=0
		call MPI_PACK(mode,1,MPI_character,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(modeindex,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nthreads,4,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		
		if(mode=='r') then
		call MPI_PACK(isweep,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nleft,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nright,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		end if
		
		call MPI_PACK(norbs,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(natoms,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nelecs,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(ncharges,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(totalSz,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_PPP,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_meanfield,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_spinreversal,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_C2,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_tree,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(subM,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(sweeps,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nstate,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(exscheme,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nweight(1),nstate,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(blocks,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nbonds,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(bondlink(1,1),norbs*norbs,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nuclQ,natoms,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(coord(1,1),3*natoms,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(t(1,1),norbs*norbs,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(hubbardU(1),norbs,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(diagmethod,20,MPI_CHARACTER,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		write(*,*) "packsizedefine=",packsize,"packbufsize=",position1
	end if
	
	call MPI_BCAST(position1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(packbuf,position1,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

	if(myid/=0) then
		position1=0
		call MPI_UNPACK(packbuf,packsize,position1,mode,1,MPI_character,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,modeindex,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nthreads,4,MPI_integer4,MPI_COMM_WORLD,ierr)
		if(mode=='r') then
			call MPI_UNPACK(packbuf,packsize,position1,isweep,1,MPI_integer4,MPI_COMM_WORLD,ierr)
			call MPI_UNPACK(packbuf,packsize,position1,nleft,1,MPI_integer4,MPI_COMM_WORLD,ierr)
			call MPI_UNPACK(packbuf,packsize,position1,nright,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		end if
		call MPI_UNPACK(packbuf,packsize,position1,norbs,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,natoms,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nelecs,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,ncharges,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,totalSz,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,logic_PPP,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,logic_meanfield,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,logic_spinreversal,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,logic_C2,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,logic_tree,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,subM,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,sweeps,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nstate,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,exscheme,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		
		allocate(nweight(nstate),stat=error)
		if(error/=0) stop
		call MPI_UNPACK(packbuf,packsize,position1,nweight,nstate,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,blocks,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nbonds,1,MPI_integer4,MPI_COMM_WORLD,ierr)

		allocate(bondlink(norbs,norbs),stat=error)
		if(error/=0) stop
		allocate(nuclQ(natoms),stat=error)
		if(error/=0) stop
		allocate(coord(3,natoms),stat=error)
		if(error/=0) stop
		allocate(hubbardU(norbs),stat=error)
		if(error/=0) stop
		allocate(t(norbs,norbs),stat=error)
		if(error/=0) stop
		
		call MPI_UNPACK(packbuf,packsize,position1,bondlink(1,1),norbs*norbs,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nuclQ(1),natoms,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,coord(1,1),3*natoms,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,t(1,1),norbs*norbs,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,hubbardU(1),norbs,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,diagmethod,20,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
		
		write(*,*) myid,"getpacksize=",position1
	end if

	allocate(pppV(norbs,norbs),stat=error)
	if(error/=0) stop
	pppV=0.0D0
	call Ohno_Potential
	
	! realnelecs is the real electrons in the system 
	realnelecs=nelecs+ncharges

!================================================================
	if(myid==1) then
		write(*,*) "----------the input information---------"
		write(*,*) "mode,modeindex=",mode,modeindex
		write(*,*) "nthreads=",nthreads
		write(*,*) "norbs=",norbs
		write(*,*) "natoms=",natoms
		write(*,*) "nelectrons=",nelecs
		write(*,*) "nextracharges=",ncharges
		write(*,*) "totalSz=",totalsz
		write(*,*) "logic_PPP=",logic_PPP
		write(*,*) "logic_MeanField=",logic_meanfield
		write(*,*) "logic_spinreversal=",logic_spinreversal
		write(*,*) "logic_C2=",logic_C2
		write(*,*) "logic_tree=",logic_tree
		write(*,*) "subM=",subM
		write(*,*) "sweeps=",sweeps
		write(*,*) "nstates=",nstate
		write(*,*) "exscheme=",exscheme
		write(*,*) "Diagonalization method=",diagmethod
		write(*,*) "nweight=",nweight
		write(*,*) "energythresh",energythresh ! other process do not know this number
		write(*,*) "nbonds=",nbonds,"bondlink="
		do i=1,norbs,1
		write(*,*) bondlink(:,i)
		end do
		write(*,*) "coord=","nuclearQ="
		do i=1,natoms,1
		write(*,*) coord(1:3,i),nuclQ(i)
		end do
		write(*,*) "hubbardU="
		write(*,*) hubbardU
		write(*,*) "oneEterm="
		write(*,*) t
		write(*,*) "pppV="
		write(*,*) pppV
		write(*,*) "----------------------------------------"
	end if
!================================================================

	deallocate(packbuf)
return
	
end Subroutine ReadInput

subroutine symboltoatomindex(symbol,atomindex,atommass)
	implicit none
	character(len=2) :: symbol
	integer :: atomindex
	real(kind=8) :: atommass
	
	if(symbol=="C") then
		atomindex=6
		atommass=12.01D0
	else if(symbol=="S") then
		atomindex=16
		atommass=32.07D0
	else if(symbol=="O") then
		atomindex=8
		atommass=16.00D0
	else if(symbol=="N") then
		atomindex=7
		atommass=14.01D0
	else 
		write(*,*) "No such atom symbol"
		stop
	end if
return

end subroutine Symboltoatomindex

!=======================================================================
!=======================================================================

subroutine centerofmass
	use variables
	use kinds_mod
	implicit none
	
	integer :: i
	real(kind=r8) :: totalmass

	do i=1,natoms,1
		cntofmass=coord(:,i)*atommass(i)+cntofmass
		totalmass=atommass(i)+totalmass
	end do
	cntofmass=cntofmass/totalmass
return

end subroutine centerofmass

!=======================================================================
!=======================================================================
