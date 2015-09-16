module Symmetry
! this module constrain the symmetry of the wavefunction
	
	use variables
	use Communicate
	implicit none

	real(kind=r8),allocatable :: symmat(:)
	integer,allocatable :: columnindex(:),rowindex(:),symmlinkcol(:)
	integer :: nsymmstate
	! symmat(rowindex(nsymmstate+1)-1),columnindex(rowindex(nsymmstate+1)-1)
	! rowindex(nsymmstate+1) :: the symmetry matrix in 3-array CSR sparse form
	! nsymmstate :: the number of symmetry state
	! symmlinkcol is the LRcoeff CSC format
	integer :: guessngoodstates
	! guess number of ngoodstates
	integer :: maxdim


contains
!===============================================================================
!===============================================================================

subroutine CreatSymmlinkbig(realdim,domain,Hindex)
! this subroutine used in system_big
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

subroutine CreatSymmlinkgood(igoodstate,leftindex,rightindex)
! get the symmlinkgood infomation

	implicit none
	
	integer :: igoodstate,leftindex,rightindex
	
	! check 
	if(igoodstate>guessngoodstates) then
		call master_print_message(igoodstate,"igoodstate>guessgoodstates")
		stop
	end if

	! in the good quantum number states space
	! get the symmetry link information
	! symmlinkgood(m,1) means the left space index in 4M basis
	! symmlinkgood(m,2) means the right space index in 4M basis
	
	! every process should know this
	symmlinkgood(igoodstate,1)=leftindex
	symmlinkgood(igoodstate,2)=rightindex
return

end subroutine CreatSymmlinkgood

!===============================================================================
!================================================================================

subroutine SymmAllocateArray
	implicit none

	integer :: error
	! guess number of ngoodstates
	guessngoodstates=16*subM*subM/5
	
	allocate(symmlinkgood(guessngoodstates,2),stat=error)
	if(error/=0) stop
	allocate(symmlinkcol(4*Rrealdim+1),stat=error)
	if(error/=0) stop
	symmlinkcol(1)=1

return
end subroutine SymmAllocateArray

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
	
	integer,allocatable ::  havetouchindex(:)
	! check every state have touched 1/0

	real(kind=r8) :: norm
	integer :: i,j,m,k
	integer :: total,phase,exist,tmp,tmp2,symmlink(ngoodstates)
	logical :: ifexist
	integer :: error

	! MPI flag
	integer :: ierr,packsize,position1,nonzero 
	character(len=1),allocatable :: packbuf(:)
	

	call master_print_message("enter symmetrymat subroutine")

	! allocate the symmat work memory
	allocate(symmat(ngoodstates),stat=error)
	if(error/=0) stop
	allocate(columnindex(ngoodstates),stat=error)
	if(error/=0) stop
	allocate(rowindex(ngoodstates+1),stat=error)
	if(error/=0) stop
	
	
	if(myid==0) then
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
	
	! when L space Sz>0 then we need to make the symmetry pair Sz<0 using same coefficient
	! when L space Sz=0 then we need to make the L+R space basis to have the specfic spin parity. Then set others to zero

		if(logic_spinreversal/=0) then
			! symmlink is the whole space ngoodstates symmlink not subspace
			symmlink=0 
			do i=1,ngoodstates,1
				if (symmlink(i)==0) then
					if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1) .and. &
						abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(i,2)) then
						symmlink(i)=i
						spinline(i,1)=i
					else
						! use the symmlinkcol to accelerate the search process
						tmp=symmlinkcol(abs(symmlinkbig(symmlinkgood(i,2),1,2)))
						tmp2=symmlinkcol(abs(symmlinkbig(symmlinkgood(i,2),1,2))+1)-1
						
						if(i<=tmp2) then
							do j=tmp,tmp2,1
								if(symmlink(j)==0 .and. &
									abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(j,1) ) then
									symmlink(i)=j
									symmlink(j)=i
									spinline(i,1)=j
									spinline(j,1)=i
									exit
								end if
							end do
						end if
					end if
				end if
			end do

			! the symmetry phase
			phase=logic_spinreversal*((-1)**(mod(nelecs/2,2)))

			do i=1,ngoodstates,1
				! the symmlink is himself
				if(symmlink(i)==i) then
					
					if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
						/=phase) then
						spinline(i,2)=0   ! the coeff should be 0
					else
						spinline(i,2)=1   ! can have a random coeff
					end if

				! either of the L or R space symmlink is not himself
				else if(symmlink(i)>i) then
					tmp=sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(symmlink(i),2),1,2))*phase
					spinline(i,2)=tmp
					spinline(symmlink(i),2)=tmp
				end if
			end do
		end if

! C2 symmetry part
! the rule is if L space fai1 fai2, corresponding R space fai3,fai4
! C2(fai1*fai4)=fai3*fai2=-1^(num3*num2)*fai2*fai3
		if(logic_C2/=0 .and. nleft==nright) then
			do i=1,ngoodstates,1
			if(C2line(i,1)==0) then
				if(symmlinkgood(i,1)==symmlinkgood(i,2)) then
					if(logic_C2*(-1)**mod(quantabigL(symmlinkgood(i,1),1)*quantabigR(symmlinkgood(i,2),1),2)==-1) then
						C2line(i,1)=i
						C2line(i,2)=0  ! the coeff should be 0
					else
						C2line(i,1)=i
						C2line(i,2)=1  ! the coeff can be random
					end if
				else 
					! use the symmlinkcol to accelerate the search process
					tmp=symmlinkcol(symmlinkgood(i,1))    ! i's left index equals j's right index
					tmp2=symmlinkcol(symmlinkgood(i,1)+1)-1

					if(i<=tmp2) then
						do j=tmp,tmp2,1  ! find the link info
							if(C2line(j,1)==0 .and. &
								symmlinkgood(j,1)==symmlinkgood(i,2) .and. &
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
								exit
							end if
						end do
					end if
				end if
			end if
			end do
		end if
		
		nsymmstate=0
		rowindex(1)=1
		havetouchindex=0  ! every state is not touched

	! the core part! find the link diagram subgroup
		do i=1,ngoodstates,1
			if(havetouchindex(i)==0) then
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
				do j=1,total,1
					havetouchindex(fulline(j,1))=1
				end do

				! update the symmmat
				if(fulline(1,2)/=0) then
					nsymmstate=nsymmstate+1
					norm=sqrt(1.0D0/DBLE(total))
					do k=1,total,1
						symmat(rowindex(nsymmstate)-1+k)=norm*DBLE(fulline(k,2)/fulline(1,2))
						columnindex(rowindex(nsymmstate)-1+k)=fulline(k,1)
						rowindex(nsymmstate+1)=rowindex(nsymmstate)+total
					end do
				end if
			end if
		end do

		!deallocate the memory space
		if(logic_spinreversal/=0) then
			deallocate(spinline)
		end if
		if(logic_C2/=0 .and. (nleft==nright)) then
			deallocate(C2line)
		end if
		deallocate(havetouchindex)
		
	end if  ! 0 process end

	! broadcast sparse matrix to other process
	call MPI_BCAST(nsymmstate,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	packsize=(ngoodstates*2+1)*4+ngoodstates*8
	allocate(packbuf(packsize))
	if(myid==0) then
		position1=0
		call MPI_PACK(rowindex,nsymmstate+1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		nonzero=rowindex(nsymmstate+1)-1
		call MPI_PACK(columnindex,nonzero,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(symmat,nonzero,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
	end if

	call MPI_BCAST(position1,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(packbuf,position1,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

	if(myid/=0) then
		position1=0
		call MPI_UNPACK(packbuf,packsize,position1,rowindex,nsymmstate+1,MPI_integer4,MPI_COMM_WORLD,ierr)
		nonzero=rowindex(nsymmstate+1)-1
		call MPI_UNPACK(packbuf,packsize,position1,columnindex,nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,symmat,nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
	end if
	
	deallocate(packbuf)

return

end subroutine SymmetryMat

!==========================================================================================
!==========================================================================================

subroutine SymmHDiag(HDIAG)
! This subroutine is to get the S+HS matrix diagonal element
!    out :: HDIAG(nsymmstate) :: the symmetrized diagonal element 

	use mathlib
	use blas95
	use F95_precision
	use mpi
	use module_sparse
	implicit none

	real(kind=r8) :: HDIAG(nsymmstate)
	! local
	integer :: dim1,isymm,nlbasis,nrbasis,indexbuf(2),operanum,&
		localnsymmstate,localnsymmstate0

	integer :: i,j,m,k,il,ir,l,ibasis
	integer :: error
	integer :: ierr,status(MPI_STATUS_SIZE),recvrequest  ! MPI_flag

	real(kind=r8),allocatable :: bufmat0(:,:,:,:,:),bufmat(:,:,:,:),bufmatlocal(:,:,:,:,:),&
	workarray(:,:),workdummy(:,:),Idummyr(:,:),Idummyl(:,:),bufmatdumm(:,:),&
	bufmatdummyr(:,:),bufmatdummyl(:,:),&
	vecdummy(:),nosymmat(:,:),localdiag(:)

	integer,allocatable :: indexall(:,:,:),indexsymm(:,:),displs(:),sendcount(:)
	logical :: iffind
	
	call master_print_message("enter symmHDiag subroutine")
	
	if(logic_spinreversal/=0 .and. logic_C2/=0 .and. nleft==nright) then
		maxdim=4
	else
		maxdim=2
	end if

	! calculate the operator number at every process at most
	if(mod(norbs,nprocs-1)==0) then
		operanum=norbs/(nprocs-1)
	else
		operanum=norbs/(nprocs-1)+1
	end if
	
	! to store symmlink small group matrix of every process
	! every orbital operator in every nsymmstate block basis
	allocate(bufmat(maxdim,maxdim,3*operanum,nsymmstate),stat=error)
	if(error/=0) stop

	! store the index of every nsymmstate in ngoodstates basis
	allocate(indexall(maxdim+1,2,nsymmstate),stat=error)  
	! first number 1:maxdim at most maxdim state; maxdim+1 is the in fact number states
	! second number 1 means L spae; 2 means R space
	if(error/=0) stop
	indexall=0

	if(myid==0) then
		! 0 procee store every orbital operator in every nsymmstate block basis
		allocate(bufmat0(maxdim,maxdim,3,norbs+1,nsymmstate),stat=error)
		if(error/=0) stop
		! norbs+1 array  store HL and HR in every nsymmstate in process 0
	end if

	! get the subspace element index that contribute to the diagonal element
	call subspaceSymmHDIAG('L',indexall,operanum,bufmat,bufmat0)
	call subspaceSymmHDIAG('R',indexall,operanum,bufmat,bufmat0)
	
	! broadcast the indexall information from process 0
	call MPI_BCAST(indexall,2*(maxdim+1)*nsymmstate,MPI_Integer,0,MPI_COMM_WORLD,ierr)

	if(myid/=0) then
		call MPI_SEND(bufmat,3*maxdim*maxdim*operanum*nsymmstate,mpi_real8,0,1,MPI_COMM_WORLD,ierr)
	end if

	if(myid==0) then
		do i=1,nprocs-1,1
			call MPI_RECV(bufmat,3*maxdim*maxdim*operanum*nsymmstate,mpi_real8,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,status,ierr)
			do j=status(MPI_SOURCE),norbs,nprocs-1
				m=(j-1)/(nprocs-1)+1
				do k=1,nsymmstate,1
				do l=1,3,1
				do ir=1,maxdim,1
					call copy(bufmat(:,ir,3*m-3+l,k),bufmat0(:,ir,l,j,k))
				end do
				end do
				end do
			end do
		end do
	end if
	deallocate(bufmat)
	
!==================================================================
! calculate the diagonal element
! every process do part of the diagonal element

	localnsymmstate0=nsymmstate/nprocs
	if(myid==nprocs-1) then   ! nprocs-1 process do the left diagonal element
		localnsymmstate=nsymmstate-localnsymmstate0*(nprocs-1)
	else
		localnsymmstate=localnsymmstate0
	end if
	
	allocate(bufmatlocal(maxdim,maxdim,3,norbs+1,localnsymmstate),stat=error)
	if(error/=0) stop
	allocate(localdiag(localnsymmstate),stat=error)
	if(error/=0) stop
	allocate(sendcount(nprocs),stat=error)
	if(error/=0) stop
	allocate(displs(nprocs),stat=error)
	if(error/=0) stop
	allocate(bufmatdumm(maxdim,maxdim),stat=error)
	if(error/=0) stop
	
	do i=1,nprocs,1
		sendcount(i)=localnsymmstate0*3*maxdim*maxdim*(norbs+1)
		displs(i)=localnsymmstate0*3*maxdim*maxdim*(norbs+1)*(i-1)
	end do
	sendcount(nprocs)=(nsymmstate-localnsymmstate0*(nprocs-1))*3*maxdim*maxdim*(norbs+1)

	call MPI_SCATTERV(bufmat0,sendcount,displs,MPI_REAL8,&
		bufmatlocal,localnsymmstate*3*maxdim*maxdim*(norbs+1),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
	
	if(myid==0) then
		deallocate(bufmat0)
	end if

	! caculate Sij * Hjk *Ski
	! every process do the local part
	do ibasis=1,localnsymmstate,1
		! the index of the global nsymmstate
		isymm=myid*localnsymmstate0+ibasis
		
		nlbasis=indexall(maxdim+1,1,isymm)  
		nrbasis=indexall(maxdim+1,2,isymm)  
		! nlbasis and nrbasis is the number of basis in L and R space in isymm block
		dim1=rowindex(isymm+1)-rowindex(isymm)

		allocate(nosymmat(dim1,dim1),stat=error)  ! store Hjk
		if(error/=0) stop
		allocate(vecdummy(dim1),stat=error)   ! store Hjk*Ski
		if(error/=0) stop

		allocate(bufmatdummyl(nlbasis,nlbasis),stat=error)
		if(error/=0) stop
		bufmatdummyl=bufmatlocal(1:nlbasis,1:nlbasis,1,norbs+1,ibasis)  ! HL 
		allocate(bufmatdummyr(nrbasis,nrbasis),stat=error)
		if(error/=0) stop
		bufmatdummyr=bufmatlocal(1:nrbasis,1:nrbasis,2,norbs+1,ibasis)  ! HR

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
						bufmatdumm=transpose(bufmatlocal(:,:,k,j,ibasis))
						call directproduct(bufmatlocal(1:nlbasis,1:nlbasis,k,i,ibasis),nlbasis,bufmatdumm(1:nrbasis,1:nrbasis),nrbasis,&
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
				call directproduct(bufmatlocal(1:nlbasis,1:nlbasis,3,i,ibasis),nlbasis,bufmatlocal(1:nrbasis,1:nrbasis,3,j,ibasis),nrbasis,&
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
		localdiag(ibasis)=dot(symmat(rowindex(isymm):rowindex(isymm+1)-1),vecdummy) !Sij *Hjk *Ski
		
		deallocate(vecdummy)
		deallocate(nosymmat)
		deallocate(Idummyr)
		deallocate(Idummyl)
		deallocate(bufmatdummyr)
		deallocate(bufmatdummyl)
		deallocate(workarray)
		deallocate(workdummy)
	end do
	
	! gather the diagonal element to 0 process
	do i=1,nprocs,1
		sendcount(i)=localnsymmstate0
		displs(i)=localnsymmstate0*(i-1)
	end do
	sendcount(nprocs)=nsymmstate-localnsymmstate0*(nprocs-1)

	call MPI_GATHERV(localdiag,localnsymmstate,MPI_REAL8,&
		Hdiag,sendcount,displs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

	deallocate(indexall)
	deallocate(bufmatlocal,localdiag)
	deallocate(displs,sendcount)
	deallocate(bufmatdumm)
	

return
end subroutine SymmHDiag

!====================================================================
!====================================================================

subroutine subspaceSymmHDIAG(domain,indexall,operanum,bufmat,bufmat0)
	
	use mathlib
	use module_sparse
	implicit none

	character(len=1) :: domain
	integer :: indexall(maxdim+1,2,nsymmstate),operanum
	real(kind=r8) :: bufmat(maxdim,maxdim,3*operanum,nsymmstate),bufmat0(maxdim,maxdim,3,norbs+1,nsymmstate)
	! local
	integer :: orbstart,orbend,Hindex,dim1,operaindex
	integer :: index1(4),m,k,is,i,j,mr,mc
	logical :: ifexist


	if(domain=='L') then
		orbstart=1
		orbend=nleft+1
		Hindex=1
		dim1=Lrealdim
	else if(domain=='R') then
		orbstart=norbs-nright
		orbend=norbs
		Hindex=2
		dim1=Rrealdim
	end if
	
	do is=orbstart,orbend,1
	if(myid==orbid1(is,1) .or. (myid==0 .and. (is==1 .or.is==norbs))) then
		operaindex=orbid1(is,2)
		do i=1,nsymmstate,1
			index1=0
			m=0
			do j=rowindex(i),rowindex(i+1)-1,1   ! i state element index in sparse form
				ifexist=.false.
				do k=1,m,1
					if(symmlinkgood(columnindex(j),Hindex)==index1(k)) then
						ifexist=.true.
						exit
					end if
				end do
				if(ifexist==.false.) then
					m=m+1
					index1(m)=symmlinkgood(columnindex(j),Hindex)
				end if
			end do
			
			do j=1,m,1
				indexall(j,Hindex,i)=index1(j)
			end do
			indexall(maxdim+1,Hindex,i)=m
			
			do mc=1,m,1
			do mr=1,m,1
				! copy the operamatbig matrix to bufmat and
				! bufmat0(:,:,1:2,norbs+1,:)
				if(myid/=0) then
					do j=3*operaindex-2,3*operaindex,1
						call SpMatIJ(4*dim1,index1(mr),index1(mc),operamatbig1(:,j),bigcolindex1(:,j),bigrowindex1(:,j),&
							bufmat(mr,mc,j,i))
					end do
				else
					call SpMatIJ(4*dim1,index1(mr),index1(mc),Hbig(:,Hindex),Hbigcolindex(:,Hindex),Hbigrowindex(:,Hindex),&
						bufmat0(mr,mc,Hindex,norbs+1,i))
				end if
			end do
			end do
		end do
	end if
	end do
return

end subroutine subspaceSymmHdiag

!==========================================================================================
!==========================================================================================

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
	formation(3)='N'
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
	operation(3)='N'
	operation(4)='F'
	
	call mkl_dcsrmm('N',nsymmstate,nnosymmstate,nnosymmstate,1.0D0,operation,symmat,columnindex,rowindex(1:nsymmstate),&
	rowindex(2:nsymmstate+1),bigmat,nnosymmstate,0.0D0,bufmat,nsymmstate)
	bufmat2=transpose(bufmat)
	call mkl_dcsrmm('N',nsymmstate,nsymmstate,nnosymmstate,1.0D0,operation,symmat,columnindex,rowindex(1:nsymmstate),&
	rowindex(2:nsymmstate+1),bufmat2,nnosymmstate,0.0D0,smallmat,nsymmstate)

return

end subroutine SymmetrizeMatrix

!====================================================================

subroutine destroysymm

	implicit none
	
	deallocate(symmat)
	deallocate(rowindex)
	deallocate(columnindex)
	deallocate(symmlinkgood)
	deallocate(symmlinkcol)
return
end subroutine destroysymm

!====================================================================
end module Symmetry
