module Symmetry
	
	use variables
	use Communicate
	implicit none

	real(kind=r8),allocatable :: symmat(:)
	integer,allocatable :: columnindex(:),rowindex(:)
	integer :: nsymmstate
	! symmat(rowindex(nsymmstate+1)-1),columnindex(rowindex(nsymmstate+1)-1)
	! rowindex(nsymmstate+1) :: the symmetry matrix in 3-array CSR sparse form
	! nsymmstate :: the number of symmetry state


contains
!===============================================================================

subroutine CreatSymmlinkbig(realdim,domain,Hindex)
	implicit none

	integer :: realdim,Hindex
	character(len=1) :: domain

	! local
	integer :: i
	
	if(domain=='R' .and. logic_C2==0) then
		do i=1,realdim,1
			symmlinkbig((i-1)*4+1,1,Hindex)=((abs(symmlinksma(i,1,Hindex))-1)*4+1)*&
				sign(1,symmlinksma(i,1,Hindex))
			symmlinkbig((i-1)*4+2,1,Hindex)=((abs(symmlinksma(i,1,Hindex))-1)*4+3)*&
				sign(1,symmlinksma(i,1,Hindex))
			symmlinkbig((i-1)*4+3,1,Hindex)=((abs(symmlinksma(i,1,Hindex))-1)*4+2)*&
				sign(1,symmlinksma(i,1,Hindex))
			symmlinkbig((i-1)*4+4,1,Hindex)=((abs(symmlinksma(i,1,Hindex))-1)*4+4)*&
				sign(1,symmlinksma(i,1,Hindex))*(-1)
		end do
	else 
		symmlinkbig(1:realdim,1,Hindex)=symmlinksma(1:realdim,1,Hindex)
		symmlinkbig(realdim+1:2*realdim,1,Hindex)=(abs(symmlinksma(1:realdim,1,Hindex))+2*&
			realdim)*sign(1,symmlinksma(1:realdim,1,Hindex))
		symmlinkbig(2*realdim+1:3*realdim,1,Hindex)=(abs(symmlinksma(1:realdim,1,Hindex))+&
			realdim)*sign(1,symmlinksma(1:realdim,1,Hindex))
		symmlinkbig(3*realdim+1:4*realdim,1,Hindex)=(abs(symmlinksma(1:realdim,1,Hindex))+3*&
			realdim)*sign(1,symmlinksma(1:realdim,1,Hindex))*(-1)
	end if
return
end subroutine CreatSymmlinkbig

!===============================================================================
!===============================================================================

subroutine CreatSymmlinkgood
	
	use exit_mod
	implicit none
	
	integer :: m,i,j
	integer :: error

	allocate(symmlinkgood(ngoodstates,2),stat=error)
	if(error/=0) stop
	! in the good quantum number states space
	! get the symmetry link information
	! symmlinkgood(m,1) means the left space index in 4M basis
	! symmlinkgood(m,2) means the right space index in 4M basis

	m=0
	do i=1,4*Rrealdim,1
	do j=1,4*Lrealdim,1
	if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
		quantabigL(j,2)+quantabigR(i,2)==totalSz) then
		m=m+1
		symmlinkgood(m,1)=j
		symmlinkgood(m,2)=i
		end if
	end do
	end do
	if(m/=ngoodstates) then
		call master_print_message(m,"symmlinkgood m/=ngoodstates")
		call exit_DMRG(sigAbort,"symmlinkgood m/=ngoodstates")
	end if
return

end subroutine CreatSymmlinkgood

!===============================================================================
!================================================================================

subroutine SymmetryMat
	! this subroutine is to contruct the symmetry matrix in sparse matrix
	! form using the symmlink information
	
	use mpi
	implicit none

	! the XXline is the link information
	integer,allocatable :: spinline(:,:),C2line(:,:)
	integer :: fulline(4,2)
	! 3-array CSR sparse form workarray(enough to store)
	real(kind=r8),allocatable :: symmat1(:)
	integer,allocatable :: columnindex1(:),rowindex1(:)
	
	integer,allocatable ::  havetouchindex(:)
	real(kind=r8) :: norm
	integer :: i,j,m,total,havetouch,k,exist
	logical :: done,iftouch,ifexist
	
	integer :: error
	integer :: ierr ! MPI flag
	
	call master_print_message("enter symmetrymat subroutine")

! allocate the work memory
	if(myid==0) then
		allocate(symmat1(ngoodstates),stat=error)
		if(error/=0) stop
		allocate(columnindex1(ngoodstates),stat=error)
		if(error/=0) stop
		allocate(rowindex1(ngoodstates+1),stat=error)
		if(error/=0) stop
	
		if(logic_spinreversal/=0) then
			allocate(spinline(ngoodstates,2),stat=error)
			if(error/=0) stop
			spinline=0
		end if
		if(logic_C2/=0 .and. nleft==nright) then
			allocate(C2line(ngoodstates,2),stat=error)
			if(error/=0) stop
			C2line=0
		end if
		
		allocate(havetouchindex(ngoodstates),stat=error)
		if(error/=0) stop
		havetouchindex=0
	
	! when L space Sz>0 then we need to make the symmetry pair Sz<0 using same coefficient
	! when L space Sz=0 then we need to make the L+R space basis to have the specfic spin parity. Then set others to zero

		if(logic_spinreversal/=0) then
			do i=1,ngoodstates,1
				! the symmlink is himself
				if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1) .and. &
					abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(i,2)) then
					if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
					/=logic_spinreversal*((-1)**(mod(nelecs/2,2)))) then
						spinline(i,1)=i
						spinline(i,2)=0   ! the coeff should be 0
					else
						spinline(i,1)=i
						spinline(i,2)=1   ! can have a random coeff
					end if
				! either of the L or R space symmlink is not himself
				else if((abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1) .or. &
						abs(symmlinkbig(symmlinkgood(i,2),1,2))/=symmlinkgood(i,2)) &
						.and. spinline(i,1)==0) then
						! spinline(i,1)==0 means it is not be touched
					done=.false.
					do j=1,ngoodstates,1
						if(spinline(j,1)==0) then     ! find the link information
							if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(j,1) .and. &
								abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(j,2)) then
								done=.true.
								if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
								*logic_spinreversal*((-1)**(mod(nelecs/2,2)))==-1) then
									spinline(i,1)=j
									spinline(j,1)=i
									spinline(i,2)=-1
									spinline(j,2)=-1
								else if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
								*logic_spinreversal*((-1)**(mod(nelecs/2,2)))==1) then
									spinline(i,1)=j
									spinline(j,1)=i
									spinline(i,2)=1
									spinline(j,2)=1
								else
									call master_print_message("failed! spinreversal symmetry problem!")
									write(*,*) sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
									*logic_spinreversal*((-1)**(mod(nelecs/2,2)))
									stop
								end if
								exit
							end if
						end if
					end do
					if(done==.false.) then  ! not find the another link part
						write(*,*) "-----------------------------------------------------------------------------"
						write(*,*) "in symmetrymat subroutine did't find the correspond state",&
							i,symmlinkgood(i,1),symmlinkgood(i,2),"corrsponds",&
							symmlinkbig(symmlinkgood(i,1),1,1),symmlinkbig(symmlinkgood(i,2),1,2)
						write(*,*) "-----------------------------------------------------------------------------"
						stop
					end if
				end if
			end do
		end if

! C2 symmetry part
! the rule is if L space fai1 fai2, corresponding R space fai3,fai4
! C2(fai1*fai4)=fai3*fai2=-1^(num3*num2)*fai2*fai3
		if(logic_C2/=0 .and. nleft==nright) then
			do i=1,ngoodstates,1
				if(symmlinkgood(i,1)==symmlinkgood(i,2)) then
					if(logic_C2*(-1)**mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)==-1) then
						C2line(i,1)=i
						C2line(i,2)=0  ! the coeff should be 0
					else
						C2line(i,1)=i
						C2line(i,2)=1  ! the coeff can be random
					end if
				else
					done=.false.
					do j=1,ngoodstates,1  ! find the link info
						if(symmlinkgood(j,1)==symmlinkgood(i,2) .and. &
							symmlinkgood(j,2)==symmlinkgood(i,1)) then
							if(logic_C2*(-1)**mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)==-1) then
								C2line(i,1)=j
								C2line(j,1)=i
								C2line(i,2)=-1
								C2line(j,2)=-1
							else
								C2line(i,1)=j
								C2line(j,1)=i
								C2line(i,2)=1
								C2line(j,2)=1
							end if
							done=.true.
							exit
						end if
					end do
					if(done==.false.) then
						call master_print_message("in symmetrymat C2 symmetry fails!")
						stop
					end if
				end if
			end do
		end if
	
		nsymmstate=0
		rowindex1(1)=1
		havetouch=0
	! the core part! find the link diagram subgroup
		do i=1,ngoodstates,1
			iftouch=.false.
			do k=1,havetouch,1
				if(i==havetouchindex(k)) then
					iftouch=.true.
					exit
				end if
			end do
			if(iftouch==.false.) then
				m=1
				total=1
				fulline(1,1)=i   ! index   ! at most 4 state share the link info
				fulline(1,2)=1   ! phase
				do while(.true.)
					if(logic_C2/=0 .and. nleft==nright) then
						ifexist=.false.
						do j=1,total,1
							if(C2line(fulline(m,1),1)==fulline(j,1)) then
								ifexist=.true.
								exist=j
								exit
							end if
						end do
						if(ifexist==.false.) then
							fulline(total+1,1)=C2line(fulline(m,1),1)
							fulline(total+1,2)=fulline(m,2)*C2line(fulline(m,1),2)
							total=total+1
						else if(C2line(fulline(m,1),2)*fulline(m,2)/=fulline(exist,2)) then
							! not conform in the link group every coeff set to be zero
							do j=1,total,1
								fulline(j,2)=0
							end do
						end if
					end if

					if(logic_spinreversal/=0) then
						ifexist=.false.
						do j=1,total,1
							if(spinline(fulline(m,1),1)==fulline(j,1)) then
								ifexist=.true.
								exist=j
								exit
							end if
						end do
						if(ifexist==.false.) then
							fulline(total+1,1)=spinline(fulline(m,1),1)
							fulline(total+1,2)=fulline(m,2)*spinline(fulline(m,1),2)
							total=total+1
						else if(spinline(fulline(m,1),2)*fulline(m,2)/=fulline(exist,2)) then
							do j=1,total,1
								fulline(j,2)=0
							end do
						end if
					end if
					
					m=m+1  ! update the number have touched 
					if(m>total) exit
				end do

				! update the havetouch index
				havetouchindex(havetouch+1:havetouch+total)=fulline(1:total,1)
				havetouch=havetouch+total
				! update the symmmat
				if(fulline(1,2)/=0) then
					nsymmstate=nsymmstate+1
					norm=sqrt(1.0D0/DBLE(total))
					do k=1,total,1
						symmat1(rowindex1(nsymmstate)-1+k)=norm*DBLE(fulline(k,2)/fulline(1,2))
						columnindex1(rowindex1(nsymmstate)-1+k)=fulline(k,1)
						rowindex1(nsymmstate+1)=rowindex1(nsymmstate)+total
					end do
				end if
			end if
		end do
		! copy the bigger space to the correct space of the sparse matrix
		allocate(symmat(rowindex1(nsymmstate+1)-1),stat=error)
		if(error/=0) stop
		allocate(columnindex(rowindex1(nsymmstate+1)-1),stat=error)
		if(error/=0) stop
		allocate(rowindex(nsymmstate+1),stat=error)
		if(error/=0) stop
		rowindex=rowindex1(1:nsymmstate+1)
		columnindex=columnindex1(1:rowindex1(nsymmstate+1)-1)
		symmat=symmat1(1:rowindex1(nsymmstate+1)-1)

		!deallocate the memory space
		if(logic_spinreversal/=0) then
			deallocate(spinline)
		end if
		if(logic_C2/=0 .and. (nleft==nright)) then
			deallocate(C2line)
		end if
		deallocate(symmat1)
		deallocate(columnindex1)
		deallocate(rowindex1)
		deallocate(havetouchindex)
	end if  ! 0 process end

! broadcast sparse matrix to other process
	call MPI_BCAST(nsymmstate,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	if(myid/=0) then
		allocate(rowindex(nsymmstate+1),stat=error)
		if(error/=0) stop
	end if
	call MPI_BCAST(rowindex,nsymmstate+1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	if(myid/=0) then
		allocate(columnindex(rowindex(nsymmstate+1)-1),stat=error)
		if(error/=0) stop
		allocate(symmat(rowindex(nsymmstate+1)-1),stat=error)
		if(error/=0) stop
	end if
	call MPI_BCAST(columnindex,rowindex(nsymmstate+1)-1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(symmat,rowindex(nsymmstate+1)-1,MPI_real8,0,MPI_COMM_WORLD,ierr)

return
end subroutine SymmetryMat

!==========================================================================================

subroutine SymmHDiag(HDIAG)
! This subroutine is to get the S+HS matrix diagonal element
!    out :: HDIAG(nsymmstate) :: the symmetrized diagonal element 

	use mathlib
	use blas95
	use F95_precision
	use mpi
	implicit none

	real(kind=r8) :: HDIAG(nsymmstate)
	! local
	integer :: operaindex,index1(4),dim1,havesend,&
		sender,isymm,nlbasis,nrbasis,indexbuf(2),operanum

	integer :: is,i,j,m,mr,mc,k,il,ir
	integer :: error
	integer :: ierr,status(MPI_STATUS_SIZE)  ! MPI_flag

	real(kind=r8),allocatable :: bufmat0(:,:,:,:,:),bufmat(:,:,:,:),bufH(:,:,:,:),&
	workarray(:,:),workdummy(:,:),Idummyr(:,:),Idummyl(:,:),bufmatdumm(:,:),&
	bufmatdummyr(:,:),bufmatdummyl(:,:),&
	vecdummy(:),nosymmat(:,:)

	integer,allocatable :: indexall(:,:,:),indexsymm(:,:)
	logical :: ifexist,iffind
	real(kind=r8) :: output

	call master_print_message("enter symmHDiag subroutine")
	
	! calculate the operator number at every process at most
	if(mod(norbs,nprocs-1)==0) then
		operanum=norbs/(nprocs-1)
	else
		operanum=norbs/(nprocs-1)+1
	end if
	
	! to store symmlink small group matrix of every process
	! every orbital operator in every nsymmstate block basis
	allocate(bufmat(4,4,3*operanum,nsymmstate),stat=error)
	if(error/=0) stop
	bufmat=0.0D0

	! store the index of every nsymmstate in ngoodstates basis
	allocate(indexall(5,2,nsymmstate),stat=error)  
	! first number 1:4 at most 4 state; 5 is the in fact number states
	! second number 1 means L spae; 2 means R space
	if(error/=0) stop
	indexall=0

	if(myid==0) then
		! 0 procee store every orbital operator in every nsymmstate block basis
		allocate(bufmat0(4,4,3,norbs,nsymmstate),stat=error)
		if(error/=0) stop
		bufmat0=0.0D0

		allocate(bufH(4,4,2,nsymmstate),stat=error)
		if(error/=0) stop
		bufH=0.0D0
		! bufH store HL and HR in every nsymmstate in process 0
	end if
!===================================================================================
	! get L space block matrix correspond to one single nsymmstate, at most 4*4
	do is=1,nleft+1,1
	if(myid==orbid(is) .or. (myid==0 .and. is==1)) then
		if(mod(is,nprocs-1)==0) then
			operaindex=is/(nprocs-1)
		else
			operaindex=is/(nprocs-1)+1
		end if
		do i=1,nsymmstate,1
			index1=0
			m=0
			do j=rowindex(i),rowindex(i+1)-1,1   ! i state element index in sparse form
				ifexist=.false.
				do k=1,m,1
					if(symmlinkgood(columnindex(j),1)==index1(k)) then
						ifexist=.true.
						exit
					end if
				end do
				if(ifexist==.false.) then
					m=m+1
					index1(m)=symmlinkgood(columnindex(j),1)
				end if
			end do
			
			do j=1,m,1
				indexall(j,1,i)=index1(j)
			end do
			indexall(5,1,i)=m
			
			do mc=1,m,1
			do mr=1,m,1
				! copy the operamatbig matrix to bufmat and bufH
				if(myid/=0) then
					bufmat(mr,mc,3*operaindex-2:3*operaindex,i)=&
						operamatbig(index1(mr),index1(mc),3*operaindex-2:3*operaindex)
				else
					bufH(mr,mc,1,i)=Hbig(index1(mr),index1(mc),1)
				end if
			end do
			end do
		end do
	end if
	end do

	! R space the same as L space algrithom
	do is=norbs,norbs-nright,-1
	if(myid==orbid(is) .or. (myid==0 .and. is==norbs)) then
		if(mod(is,nprocs-1)==0) then
			operaindex=is/(nprocs-1)
		else
			operaindex=is/(nprocs-1)+1
		end if
		do i=1,nsymmstate,1
			index1=0
			m=0
			do j=rowindex(i),rowindex(i+1)-1,1
				ifexist=.false.
				do k=1,m,1
					if(symmlinkgood(columnindex(j),2)==index1(k)) then
						ifexist=.true.
						exit
					end if
				end do
				if(ifexist==.false.) then
					m=m+1
					index1(m)=symmlinkgood(columnindex(j),2)
				end if
			end do

			do j=1,m,1
				indexall(j,2,i)=index1(j)
			end do
			indexall(5,2,i)=m

			do mc=1,m,1
			do mr=1,m,1
				if(myid/=0) then
					bufmat(mr,mc,3*operaindex-2:3*operaindex,i)=&
						operamatbig(index1(mr),index1(mc),3*operaindex-2:3*operaindex)
				else
					bufH(mr,mc,2,i)=Hbig(index1(mr),index1(mc),2)
				end if
			end do
			end do
		end do
	end if
	end do
!===============================================================================================
	
! broadcast the indexall information from process 0
	call MPI_BCAST(indexall,10*nsymmstate,MPI_Integer,0,MPI_COMM_WORLD,ierr)

	if(myid/=0) then
		call MPI_SEND(bufmat,48*operanum*nsymmstate,mpi_real8,0,1,MPI_COMM_WORLD,ierr)
	end if

	if(myid==0) then
		do i=1,nprocs-1,1
			call MPI_RECV(bufmat,48*operanum*nsymmstate,mpi_real8,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,status,ierr)
			m=0
			do j=status(MPI_SOURCE),norbs,nprocs-1
				m=m+1
				do k=1,nsymmstate,1
					bufmat0(:,:,:,j,k)=bufmat(:,:,3*m-2:3*m,k)
				end do
			end do
		end do
	end if

	deallocate(bufmat)

	allocate(bufmat(4,4,3,norbs+1),stat=error)
	if(error/=0) stop

	if(myid/=0) then
		allocate(bufmatdumm(4,4),stat=error)
		if(error/=0) stop
	end if

! master-slaver mode to calculate the diagnol element, from Huang Xiaomeng's PPT
! 0 process send every symmstate calculation to slaver process
	if(myid==0) then
		HDIAG=0.0D0
		do i=1,nprocs-1,1
			bufmat(:,:,:,1:norbs)=bufmat0(:,:,:,:,i)
			bufmat(:,:,1:2,norbs+1)=bufH(:,:,:,i)
			call MPI_SEND(bufmat,48*(norbs+1),mpi_real8,i,i,MPI_COMM_WORLD,ierr)
		end do

		havesend=nprocs-1
		do i=1,nsymmstate,1
			call MPI_RECV(output,1,MPI_real8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
			sender=status(MPI_SOURCE)
			HDIAG(status(MPI_TAG))=output+HDIAG(status(MPI_TAG))

			if(havesend<nsymmstate) then
				bufmat(:,:,:,1:norbs)=bufmat0(:,:,:,:,havesend+1)
				bufmat(:,:,1:2,norbs+1)=bufH(:,:,:,havesend+1)
				call MPI_SEND(bufmat,48*(norbs+1),mpi_real8,sender,havesend+1,MPI_COMM_WORLD,ierr)
				havesend=havesend+1
			else
				call MPI_SEND(1.0,0,mpi_real8,sender,0,MPI_COMM_WORLD,ierr)
			end if
		end do
	end if


	! caculate Sij * Hjk *Ski
	if(myid/=0) then
   300	call MPI_RECV(bufmat,48*(norbs+1),mpi_real8,0,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
		
		isymm=status(MPI_TAG)
	! if isymm==0 no work
	if(isymm/=0) then
		nlbasis=indexall(5,1,isymm)  
		nrbasis=indexall(5,2,isymm)  
		! nlbasis and nrbasis is the number of basis in L and R space in isymm block
		dim1=rowindex(isymm+1)-rowindex(isymm)

		allocate(nosymmat(dim1,dim1),stat=error)  ! store Hjk
		if(error/=0) stop
		allocate(vecdummy(dim1),stat=error)   ! store Hjk*Ski
		if(error/=0) stop

		allocate(bufmatdummyl(nlbasis,nlbasis),stat=error)
		if(error/=0) stop
		bufmatdummyl=bufmat(1:nlbasis,1:nlbasis,1,norbs+1)  ! HL 
		allocate(bufmatdummyr(nrbasis,nrbasis),stat=error)
		if(error/=0) stop
		bufmatdummyr=bufmat(1:nrbasis,1:nrbasis,2,norbs+1)  ! HR

		allocate(Idummyr(nrbasis,nrbasis),stat=error) ! identity matrix Ir
		if(error/=0) stop
		Idummyr=0.0D0
		do i=1,nrbasis,1
			Idummyr(i,i)=1.0D0
		end do              
		allocate(Idummyl(nlbasis,nlbasis),stat=error) ! identity matrix Il
		if(error/=0) stop
		Idummyl=0.0D0
		do i=1,nlbasis,1
			Idummyl(i,i)=1.0D0
		end do  

		allocate(workdummy(nlbasis*nrbasis,nlbasis*nrbasis),stat=error)
		if(error/=0) stop
		allocate(workarray(nlbasis*nrbasis,nlbasis*nrbasis),stat=error)
		if(error/=0) stop
		workarray=0.0D0  ! store Hjk
		workdummy=0.0D0

		! HL and HR contribute
		call directproduct(bufmatdummyl,nlbasis,Idummyr,nrbasis,&
			workdummy)
		workarray=workdummy
		call directproduct(Idummyl,nlbasis,bufmatdummyr,nrbasis,&
			workdummy)
		workarray=workdummy+workarray

		! transfer integeral and PPPV
		do i=1,nleft+1,1
			do j=norbs,norbs-nright,-1
				if(bondlink(i,j)==1) then
					do k=1,2,1
						bufmatdumm=transpose(bufmat(:,:,k,j))
						call directproduct(bufmat(1:nlbasis,1:nlbasis,k,i),nlbasis,bufmatdumm(1:nrbasis,1:nrbasis),nrbasis,&
							workdummy)
						! add phase
						do ir=1,nrbasis,1
						do il=1,nlbasis,1
							workdummy(:,(ir-1)*nlbasis+il)=workdummy(:,(ir-1)*nlbasis+il)&
								*((-1.0D0)**(mod(quantabigL(indexall(il,1,isymm),1),2)))
						end do
						end do
						workarray=(workdummy+transpose(workdummy))*t(i,j)+workarray
					end do
				end if
				! PPP term
				call directproduct(bufmat(1:nlbasis,1:nlbasis,3,i),nlbasis,bufmat(1:nrbasis,1:nrbasis,3,j),nrbasis,&
					workdummy)
				workarray=workdummy*pppV(i,j)+workarray
			end do
		end do
		
		allocate(indexsymm(2,nrbasis*nlbasis),stat=error)
		if(error/=0) stop
		! store the initial index of every state
		do j=1,nrbasis,1
		do k=1,nlbasis,1
			indexsymm(1,(j-1)*nlbasis+k)=k
			indexsymm(2,(j-1)*nlbasis+k)=j
		end do
		end do

		do i=rowindex(isymm),rowindex(isymm+1)-1,1
			iffind=.false.
			! find the correspond state to the symmmat
			! other state is no use (nlbasis*nrbasis>dim1)
			do j=1,nrbasis*nlbasis,1
				if((symmlinkgood(columnindex(i),1)==indexall(indexsymm(1,j),1,isymm)) &
				.and. (symmlinkgood(columnindex(i),2)==indexall(indexsymm(2,j),2,isymm))) then
					if(j/=i-rowindex(isymm)+1) then
						call swap(workarray(j,:),workarray(i-rowindex(isymm)+1,:))
						call swap(workarray(:,j),workarray(:,i-rowindex(isymm)+1))
						! swap the index
						indexbuf=indexsymm(:,j)
						indexsymm(:,j)=indexsymm(:,i-rowindex(isymm)+1)
						indexsymm(:,i-rowindex(isymm)+1)=indexbuf
					end if
					iffind=.true.
					exit
				end if
			end do
			if(iffind/=.true.) then
				write(*,*) "======================================="
				write(*,*) "symmdiag not correspond!"
				write(*,*) "======================================="
				stop
			end if
		end do
		deallocate(indexsymm)
		nosymmat=workarray(1:dim1,1:dim1)  !Hjk
		call gemv(nosymmat,symmat(rowindex(isymm):rowindex(isymm+1)-1),vecdummy,1.0D0,0.0D0,'N')  ! Hjk*Ski
		output=dot(symmat(rowindex(isymm):rowindex(isymm+1)-1),vecdummy) !Sij *Hjk *Ski
		
		deallocate(vecdummy)
		deallocate(nosymmat)
		deallocate(Idummyr)
		deallocate(Idummyl)
		deallocate(bufmatdummyr)
		deallocate(bufmatdummyl)
		deallocate(workarray)
		deallocate(workdummy)
		call MPI_SEND(output,1,MPI_real8,0,isymm,MPI_COMM_WORLD,ierr)
		goto 300
	end if
	end if

	deallocate(indexall)
	deallocate(bufmat)
	if(myid==0) then
		deallocate(bufH)
		deallocate(bufmat0)
	else 
		deallocate(bufmatdumm)
	end if

return
end subroutine SymmHDiag

!====================================================================

subroutine SymmetrizeState(nnosymmstate,bigvector,smallvector,operation)
! this subroutine is to symmetrize and unsymmetrize a vector
	! in/out :: bigvectorosymmstate),smallvector(nsymmstate) :: the target vector
	! in operation :: specfic the operation, symmetrize or unsymmetrize

	implicit none
	include "mkl_spblas.fi"
	integer :: nnosymmstate
	real(kind=r8) :: bigvector(nnosymmstate),smallvector(nsymmstate)
	character(len=1) :: operation
	character(len=1) :: formation(6)
	
	formation(1)='G'
	formation(2)='L'
	formation(3)='L'
	formation(4)='F'
	
	if(operation=='s') then
		call mkl_dcsrmv('N',nsymmstate,nnosymmstate,1.0D0,formation,symmat,columnindex,rowindex(1:nsymmstate),&
			rowindex(2:nsymmstate+1),bigvector,0.0D0,smallvector)
	else if(operation=='u') then
		call mkl_dcsrmv('T',nsymmstate,nnosymmstate,1.0D0,formation,symmat,columnindex,rowindex(1:nsymmstate),&
			rowindex(2:nsymmstate+1),smallvector,0.0D0,bigvector)
	else
		write(*,*) "----------------------------------"
		write(*,*) "In symmetrizestate subroutine the operation is &
		wrong",operation
		write(*,*) "----------------------------------"
	end if

return

end subroutine SymmetrizeState

!====================================================================
!====================================================================

subroutine SymmetrizeMatrix(nnosymmstate,bigmat,smallmat)
! this subroutine is to symmetrize a matrix  S+HS
	! in :: bigmat(noosymmstate**2)
	! out :: smallmat(nsymmstate**2)
! used in fullmat

	implicit none
	include "mkl_spblas.fi"
	integer :: nnosymmstate
	real(kind=r8) :: bigmat(nnosymmstate,nnosymmstate),smallmat(nsymmstate,nsymmstate)
	real(kind=r8),allocatable :: bufmat(:,:),bufmat2(:,:)
	character(len=1) :: operation(6)
	integer :: error
	
	write(*,*) "enter SymmtrizeMatrix subroutine"
	
	allocate(bufmat(nsymmstate,nnosymmstate),stat=error)
	if(error/=0) stop
	allocate(bufmat2(nnosymmstate,nsymmstate),stat=error)
	if(error/=0) stop
	
	operation(1)='G'
	operation(2)='L'
	operation(3)='L'
	operation(4)='F'
	
	call mkl_dcsrmm('N',nsymmstate,nnosymmstate,nnosymmstate,1.0D0,operation,symmat,columnindex,rowindex(1:nsymmstate),&
	rowindex(2:nsymmstate+1),bigmat,nnosymmstate,0.0D0,bufmat,nsymmstate)
	bufmat2=transpose(bufmat)
	call mkl_dcsrmm('N',nsymmstate,nsymmstate,nnosymmstate,1.0D0,operation,symmat,columnindex,rowindex(1:nsymmstate),&
	rowindex(2:nsymmstate+1),bufmat2,nnosymmstate,0.0D0,smallmat,nsymmstate)

return

end subroutine SymmetrizeMatrix

!====================================================================
subroutine destorysymm

	implicit none
	
	deallocate(symmat)
	deallocate(rowindex)
	deallocate(columnindex)
	deallocate(symmlinkgood)
return
end subroutine destorysymm

!====================================================================
end module Symmetry
