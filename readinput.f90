	Subroutine READINPUT
! this subroutine is to read in the input parameters
! process 0 read the input file and distribute to other process

	USE variables
	USE PPP_term
	USE MPI
	implicit none

	integer :: i,j,error,link1,link2
	!link1 and link2 is the bondlink atom information

	integer :: junk_natoms  
	character(len=2) :: symbol
	character(len=1),allocatable :: packbuf(:)
	integer(kind=4) :: position1
	integer(kind=4) :: packsize
	! to compare if the natoms is right 

	! read the input file
	! read file even unit
	! write file odd unit

	if(myid==0) then
		write(*,*) "enter in readinput subroutine"
!------------------------------------------------------------
	! default value
	! default value
	exscheme=0
	blocks=2
	modeindex=0
	mode='s'

	open(unit= 10,file="inp",status="old")
	read(10,*) mode,modeindex ! which mode you want to use
! including restart/debug/
	read(10,*) norbs     !how many orbitals
	read(10,*) junk_natoms   !how many atoms
	read(10,*) nelecs  !how many electrons
	read(10,*) ncharges ! how many extra charges +1 means add 1 electron
	read(10,*) totalsz !total Sz of the system
	read(10,*) logic_PPP ! if do PPP model logic_PPP=1
	read(10,*) logic_spinreversal ! if do spin reversal logic_spinreversal=+-1 
	read(10,*) logic_tree ! if do tree tensor algorithm logic_tree=1
	read(10,*) subM  ! DMRG SUB M
	read(10,*) sweeps ! DMRG how many sweeps
	read(10,*) nstate ! how many state wanted to get
	read(10,*) energythresh ! the threshold of the total energy you want to get

! sweepenergy is the total energy of every sweep
	allocate(sweepenergy(0:sweeps,nstate),stat=error)
	if(error/=0) stop
	sweepenergy=0.0D0
! 
	allocate(nweight(nstate),stat=error)
	if(error/=0) stop
	!default value
	nweight=1.0D0
	
	if(nstate/=1) then
	  read(10,*) exscheme ! exschemem=1 average method =2 my new method
	  if(exscheme==1) then
	  	read(10,*) nweight(1:nstate) ! nweight is the average DMRG excited state
	  end if
	end if

	if(logic_tree==1) then
	  read(10,*) blocks ! DMRG how many blocks . Tree tensor net work
	end if


	! if logic_spinreversal/=0 and totalsz/=0 then stop
	if (logic_spinreversal/=0 .and. totalsz/=0 .and. myid==0) then
		write(*,*) "---------------------------------"
		write(*,*) "spin reversal needs Sz=0, failed!"
		write(*,*) "---------------------------------"
		stop
	end if



!-------------------------------------------------------------
	! if logic_ppp=1 read how many bonds and link information  
	if(logic_PPP==1) then
	  read(10,*) nbonds
	  allocate(bondlink(norbs,norbs),stat=error)
	  if(error/=0) stop
	  bondlink=0
	  do i=1,nbonds,1
	    read(10,*) link1,link2
	    bondlink(link1,link2)=1
	    bondlink(link2,link1)=1
	  end do
	

	! coordiantes of the system
	  open(unit= 12,file="coord.xyz",status="old")
	  read(12,*) natoms
	  
	!if logic_PPP=1 read nuclear q (physically, chemical potential)
	  allocate(nuclQ(natoms),stat=error)
	  if(error/=0) stop
	  
	  if(natoms/=junk_natoms .and. myid==0) then
	    write(*,*) "--------------------------------------------------"
	    write(*,*) "!the natoms has some problem natoms/=junk_natoms!"
	    write(*,*) "--------------------------------------------------"
	    stop
	  end if
  
	  allocate(coord(3,1:natoms),stat=error)
	  if(error/=0) stop
  
	  read(12,*)
	  do i=1,natoms,1   
	    read(12,*) symbol,coord(1:3,i),nuclQ(i)
	  end do
	  close(12)
	end if

! integral of the system
	open(unit=14,file="integral.inp",status="old")
	if(logic_PPP==1) then
		if(natoms/=norbs .and. myid==0) then
			write(*,*) "--------------------------------------------------"
			write(*,*) "!Use PPP model, the norbs/=natoms Caution!"
			write(*,*) "--------------------------------------------------"
			stop
		end if
		allocate(hubbardU(norbs),stat=error)
		if(error/=0) stop
		allocate(t(norbs,norbs),stat=error)
		if(error/=0) stop
	! be careful about the structure of the integral formatted
		t=0.0D0
		do i=1,norbs,1
			do j=1,norbs,1
			if(bondlink(j,i)==1) then
				read(14,*) t(j,i)    ! read in every transfer integral
			end if                 
			end do
		end do
		
		do i=1,norbs,1
			bondlink(i,i)=2  ! bondlink=2 means that it is the site energy
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
	close(14)
	
	if(logic_tree==1) then
		open(unit=16,file="tree.inp",status="old")
		allocate(treelink(blocks+1,norbs),stat=error)
		if(error/=0) stop
		do i=1,norbs,1
		read(16,*) treelink(1:blocks+1,i)
		end do
		close(16)
	end if

	close(10)
	end if



	packsize=3500000
	allocate(packbuf(packsize),stat=error)
	if(error/=0) stop

	if(myid==0) then
		position1=0
		call MPI_PACK(mode,1,MPI_character,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(modeindex,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(norbs,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(natoms,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nelecs,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(ncharges,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(totalSz,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_PPP,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_spinreversal,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(logic_tree,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(subM,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(sweeps,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nstate,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(exscheme,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nweight(1:nstate),nstate,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(blocks,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nbonds,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		do i=1,norbs,1
		call MPI_PACK(bondlink(:,i),norbs,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		end do
		call MPI_PACK(nuclQ,natoms,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		do i=1,natoms,1
		call MPI_PACK(coord(:,i),3,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		end do
		do i=1,norbs,1
		call MPI_PACK(t(:,i),norbs,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		end do
		call MPI_PACK(hubbardU,norbs,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		write(*,*) "packbufsize=",position1
	end if
	
	call MPI_bcast(packbuf,packsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

	if(myid/=0) then
		position1=0
		call MPI_UNPACK(packbuf,packsize,position1,mode,1,MPI_character,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,modeindex,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,norbs,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,natoms,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nelecs,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,ncharges,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,totalSz,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,logic_PPP,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,logic_spinreversal,1,MPI_integer4,MPI_COMM_WORLD,ierr)
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
		do i=1,norbs,1
		call MPI_UNPACK(packbuf,packsize,position1,bondlink(:,i),norbs,MPI_integer4,MPI_COMM_WORLD,ierr)
		end do
		call MPI_UNPACK(packbuf,packsize,position1,nuclQ,natoms,MPI_real8,MPI_COMM_WORLD,ierr)
		do i=1,natoms,1
		call MPI_UNPACK(packbuf,packsize,position1,coord(:,i),3,MPI_real8,MPI_COMM_WORLD,ierr)
		end do
		do i=1,norbs,1
		call MPI_UNPACK(packbuf,packsize,position1,t(:,i),norbs,MPI_real8,MPI_COMM_WORLD,ierr)
		end do
		call MPI_UNPACK(packbuf,packsize,position1,hubbardU,norbs,MPI_real8,MPI_COMM_WORLD,ierr)
		
		write(*,*) myid,"getpacksize=",position1
	end if


	allocate(pppV(norbs,norbs),stat=error)
	if(error/=0) stop
	pppV=0.0D0
	call ohno_potential
	! here can add many different potentials in Module PPP term

!-----------------------------------------------------
	if(myid==0) then
		write(*,*) "----------the input information---------"
		write(*,*) "norbs=",norbs
		write(*,*) "natoms=",natoms
		write(*,*) "nelectrons=",nelecs
		write(*,*) "nextracharges=",ncharges
		write(*,*) "totalSz=",totalsz
		write(*,*) "logic_PPP=",logic_PPP
		write(*,*) "logic_spinreversal=",logic_spinreversal
		write(*,*) "logic_tree=",logic_tree
		write(*,*) "subM=",subM
		write(*,*) "sweeps=",sweeps
		write(*,*) "energythresh",energythresh
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

	deallocate(packbuf)
	return
	
	end Subroutine
