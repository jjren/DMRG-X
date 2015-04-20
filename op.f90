subroutine op(bigdim,smadim,coeff,newcoeff)
! this is the core subroutine to calculate the S*H*S*C or H*C
! the parallel schema follow JCP 12 3174(2004) garnet chan
! if want to save memory, then can write a wrapper, to send one coeff every time

!--------------------------------------------------------
! input bigdim,smadim,coeff
! bigdim is the totaldim 16M*M
! bigdim may be < 16M*M because we use good quantum number and spin symmetry,
! and the dimension  may be half
! if groud state smadim=1
! if gs+ex smadim may be >1
!--------------------------------------------------------
! output newcoeff
! coeff is the input coefficient and in the 1-d arrary format
! new coeff is H cross C result
!---------------------------------------------------------

	use mpi
	use variables
	use symmetry
	use mathlib
	use BLAS95
	use F95_PRECISION

	implicit none

	integer :: bigdim,smadim
	real(kind=r8) :: coeff(bigdim*smadim),newcoeff(bigdim*smadim)
	
	! local
	integer :: operaindex
	real(kind=r8),allocatable :: LRcoeffin(:,:,:),LRcoeffout(:,:,:),coeffnosymm(:),coeffnosymmreduce(:)
	real(kind=r8),allocatable :: hopmat(:,:,:,:),pppVmat(:,:,:),buffmat(:,:)
	real(kind=r8),allocatable :: phase(:)
	integer :: error,i,j,k,l,m
	
	logical :: ifhop,ifhopsend,ifpppVsend
	integer :: hoptouched(nprocs-1),hopntouched,pppVtouched(nprocs-1),pppVntouched
	
	! MPI flag
	integer :: status(MPI_STATUS_SIZE),hopsendrequest(nprocs-1),pppVsendrequest(nprocs-1),hoprecvrequest
	integer :: ierr

!============================================================
! allocate workspace
	! store nosymmetry coeff
	allocate(coeffnosymm(ngoodstates*smadim),stat=error)
	if(error/=0) stop
	allocate(buffmat(4*Lrealdim,4*Rrealdim),stat=error)
	if(error/=0) stop
	
	if(myid/=0) then
		do i=norbs,norbs-nright,-1
			if(myid==orbid(i)) then
				if( .not. allocated(LRcoeffin)) then
					! transform the 1-array to 4M*4M form 
					! 256M if nstate=2 M=1000
					allocate(LRcoeffin(4*Lrealdim,4*Rrealdim,smadim),stat=error)   ! coeff to LR format
					if(error/=0) stop
					exit
				end if
			end if
		end do
	
		do i=1,nleft+1,1
			if(myid==orbid(i)) then
				if( .not. allocated(LRcoeffout)) then
					allocate(LRcoeffout(4*Lrealdim,4*Rrealdim,smadim),stat=error)  ! newcoeff to LR format
					if(error/=0) stop
					exit
				end if
			end if
		end do

		do i=1,nleft+1,1
		do j=norbs,norbs-nright,-1
			if(bondlink(i,j)==1) then
				if(myid==orbid(i) .or. myid==orbid(j)) then
					if(.not. allocated(hopmat)) then
						allocate(hopmat(4*Lrealdim,4*Rrealdim,4,smadim),stat=error) ! store the hopping matrix
						if(error/=0) stop
					end if
				end if
			end if
		end do
		end do
		allocate(pppVmat(4*Lrealdim,4*Rrealdim,smadim),stat=error) ! store the pppV matrix
		if(error/=0) stop
	else 
		allocate(LRcoeffin(4*Lrealdim,4*Rrealdim,smadim),stat=error)   ! coeff to LR format
		if(error/=0) stop
		allocate(LRcoeffout(4*Lrealdim,4*Rrealdim,smadim),stat=error)   ! coeff to LR format
		if(error/=0) stop
	end if

!=================================================================================================
	
	if( myid==0 ) then
		! if symmetry==.true. then transform the symmetry coeff to the unsymmetry coeff
		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
			if(bigdim/=nsymmstate) then
				call master_print_message(bigdim,"In symmetry, op bigdim/=nsymmstate wrong!")
				stop
			end if
			do i=1,smadim,1
				call symmetrizestate(ngoodstates,coeffnosymm(ngoodstates*(i-1)+1:i*ngoodstates),&
					coeff(bigdim*(i-1)+1:i*bigdim),'u')
			end do
		else
			if(bigdim/=ngoodstates) then
				call master_print_message(bigdim,"op bigdim/=ngoodstates wrong!")
				stop
			end if
			coeffnosymm=coeff
		end if
	end if

	! send coeffnosymm to R space process
	do i=1,nprocs-1,1
		if(myid==0 .or. myid==i) then
			do j=norbs,norbs-nright,-1
				if(i==orbid(j)) then
					if(myid==0) then
						call MPI_SEND(coeffnosymm,ngoodstates*smadim,MPI_real8,i,1,MPI_COMM_WORLD,ierr)
					else if(myid==i) then
						call MPI_RECV(coeffnosymm,ngoodstates*smadim,MPI_real8,0,1,MPI_COMM_WORLD,status,ierr)
					end if
					exit
				end if
			end do
		end if
	end do

!-----------------------------------------------------------------------------

! R space process do it
! to transform the 16M*M coeff to 4M*4M(L*R) format ; coeff(16M^2,n) to coeff(4M,4M,n) 
! since the input coeff is ngoodstates and other nongoodstates sets to 0
	if(allocated(LRcoeffin)) then
		m=0
		LRcoeffin=0.0D0
		do k=1,smadim,1
			do i=1,4*Rrealdim,1
			do j=1,4*Lrealdim,1
				if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
					quantabigL(j,2)+quantabigR(i,2)==totalSz) then
					m=m+1
					LRcoeffin(j,i,k)=coeffnosymm(m)
				end if
			end do
			end do
		end do
		if(m/=smadim*ngoodstates) then
			call master_print_message("m/=smadim*ngoodstates op good quantum states number wrong!")
			stop
		end if
	end if

! L space process and 0 process initializaiton
! L space process have LRcoeffout and send to 0 process at last
	if(allocated(LRcoeffout)) then
		LRcoeffout=0.0D0   
	end if

!  calculate HL*1 and 1*HR 
	if(myid==0) then 
		do i=1,smadim,1
			do k=1,2,1
				if(k==1) then ! HL*1
					call gemm(Hbig(1:4*Lrealdim,1:4*Lrealdim,1),LRcoeffin(:,:,i),&
						buffmat,'N','N',1.0D0,0.0D0)
				else ! 1*HR
					call gemm(LRcoeffin(:,:,i),Hbig(1:4*Rrealdim,1:4*Rrealdim,2),&
					buffmat,'N','N',1.0D0,0.0D0)
				end if
				LRcoeffout(:,:,i)=buffmat+LRcoeffout(:,:,i)
			end do
		end do
	end if

!------------------------------------------------
! vlr=Hlrl'r'*Cl'r'=sum(opt,l',r')=sum(Lopt,l') parity*Oll'*sum(Ropt,r') Orr'cl'r'
! the parallel schema is that 0 process bcast the coeff matrix to other process
! and 0 process gather the result
	if(myid/=0) then
	do i=norbs,norbs-nright,-1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
!=====================================================================================
			
			! construct the pppVmat
			do j=1,smadim,1
				call gemm(LRcoeffin(:,:,j),operamatbig(1:4*Rrealdim,1:4*Rrealdim,operaindex*3)&
					,pppVmat(:,:,j),'N','N',1.0D0,0.0D0)
			end do
			
			! send the pppVmat
			pppVtouched=0
			pppVntouched=0
			do l=1,nleft+1,1
				if(orbid(l)/=myid) then ! if l==myid need not send hopmat
					ifpppVsend=.false.
					do m=1,pppVntouched,1
						if(orbid(l)==pppVtouched(m)) then
							ifpppVsend=.true.   ! have send pppVmat to this process
							exit
						end if
					end do
					if(ifpppVsend==.false.) then
						pppVntouched=pppVntouched+1
						pppVtouched(pppVntouched)=orbid(l)
						call MPI_ISEND(pppVmat,16*Lrealdim*Rrealdim*smadim,MPI_real8,orbid(l),i,MPI_COMM_WORLD,pppVsendrequest(pppVntouched),ierr)
					end if
				end if
			end do
!====================================================================
			! check if need hopping matrix
			ifhop=.false.
			do j=1,nleft+1,1
				if(bondlink(i,j)==1) then
					ifhop=.true.
					exit
				end if
			end do
			
			if(ifhop==.true.) then
				! construct hopmat
				do j=1,smadim,1
					do k=1,4,1
						if(k<=2) then
						! firstly needed to transfer from a(+) to a, so 'T' and 'T' is 'N', (ni-1)^+=(ni-1)
						! k=1 a up,k=2 a down,k=3 n,k=4 a+ up,k=5 a+ down;
							call gemm(LRcoeffin(:,:,j),operamatbig(1:4*Rrealdim,1:4*Rrealdim,(operaindex-1)*3+k)&
							,hopmat(:,:,k,j),'N','N',1.0D0,0.0D0)
						else
							call gemm(LRcoeffin(:,:,j),operamatbig(1:4*Rrealdim,1:4*Rrealdim,(operaindex-1)*3+k-2)&
							,hopmat(:,:,k,j),'N','T',1.0D0,0.0D0)
						end if
					end do
				end do
			!---------------------------------------------------------
				! the +1 -1 phase added to l' of hopmat
				allocate(phase(4*Lrealdim),stat=error)
				if(error/=0) stop

				do j=1,4*Lrealdim,1
					phase(j)=(-1.0D0)**(mod(quantabigL(j,1),2))
				end do

				do m=1,smadim,1
				do l=1,2,1
				do k=1,4*Rrealdim,1
				do j=1,4*Lrealdim,1
					hopmat(j,k,l,m)=hopmat(j,k,l,m)*phase(j)
				end do
				end do
				end do
				end do

				phase=phase*(-1.0D0)
				do m=1,smadim,1
				do l=3,4,1
				do k=1,4*Rrealdim,1
				do j=1,4*Lrealdim,1
					!transfer from al*ar^(+) to ar^(+)*al
					hopmat(j,k,l,m)=hopmat(j,k,l,m)*phase(j)
				end do
				end do
				end do
				end do

				deallocate(phase)
			!-----------------------------------------------------------

				! send the hopmat
				hoptouched=0
				hopntouched=0
				do l=1,nleft+1,1
					if(bondlink(i,l)==1 .and. orbid(l)/=myid) then ! if orbid(l)==myid need not send hopmat
						ifhopsend=.false.
						do m=1,hopntouched,1
							if(orbid(l)==hoptouched(m)) then
								ifhopsend=.true.   ! have send hopmat to this process
								exit
							end if
						end do
						if(ifhopsend==.false.) then
							hopntouched=hopntouched+1
							hoptouched(hopntouched)=orbid(l)
							call MPI_ISEND(hopmat,64*Lrealdim*Rrealdim*smadim,MPI_real8,orbid(l),i,MPI_COMM_WORLD,hopsendrequest(hopntouched),ierr)
						end if
					end if
				end do
			end if
		end if
!===============================================================================
	
! L space recv pppVmat hopmat---------------------------------------
		if(myid/=orbid(i)) then
			! pppVmat recv
			do l=1,nleft+1,1
				if(myid==orbid(l)) then
					call MPI_RECV(pppVmat,16*Lrealdim*Rrealdim*smadim,MPI_real8,orbid(i),i,MPI_COMM_WORLD,status,ierr)
					exit  ! only recv once
				end if
			end do
			! hopmat recv
			do l=1,nleft+1,1
				if(bondlink(l,i)==1 .and. myid==orbid(l)) then  
					call MPI_IRECV(hopmat,64*Lrealdim*Rrealdim*smadim,MPI_real8,orbid(i),i,MPI_COMM_WORLD,hoprecvrequest,ierr)
					exit  ! only recv once
				end if
			end do
		end if
!---------------------------------------------------------------------
! pppV calculation
		do l=1,nleft+1,1
			if(myid==orbid(l)) then
				if(mod(l,nprocs-1)==0) then
					operaindex=l/(nprocs-1)
				else
					operaindex=l/(nprocs-1)+1
				end if

				do j=1,smadim,1
					! buffmat is to save the intermediate matrix
					call gemm(operamatbig(1:4*Lrealdim,1:4*Lrealdim,operaindex*3),pppVmat(:,:,j)&
						,buffmat,'N','N',1.0D0,0.0D0)
					call ScaleMatrix(buffmat,4*Lrealdim,4*Rrealdim,pppV(i,l),'N')
					LRcoeffout(:,:,j)=buffmat+LRcoeffout(:,:,j)
				end do
			end if
		end do
!---------------------------------------------------------------------
! hopping term calculation
		! wait to recv hopmat
		if(myid/=orbid(i)) then
			do l=1,nleft+1,1
				if(bondlink(l,i)==1 .and. myid==orbid(l)) then  
					call MPI_WAIT(hoprecvrequest,status,ierr)
					exit  ! every process only wait once
				end if
			end do
		end if

		do l=1,nleft+1,1
			if(bondlink(i,l)==1 .and. myid==orbid(l)) then
				if(mod(l,nprocs-1)==0) then
					operaindex=l/(nprocs-1)
				else
					operaindex=l/(nprocs-1)+1
				end if

				do j=1,smadim,1
					do k=1,4,1
						!k<=2 al^+*ar,k>2 al*ar^(+) 
						if(k<=2) then
							call gemm(operamatbig(1:4*Lrealdim,1:4*Lrealdim,(operaindex-1)*3+k),hopmat(:,:,k,j)&
								,buffmat,'N','N',1.0D0,0.0D0)
						else
							call gemm(operamatbig(1:4*Lrealdim,1:4*Lrealdim,(operaindex-1)*3+k-2),hopmat(:,:,k,j)&
								,buffmat,'T','N',1.0D0,0.0D0)
						end if
						call ScaleMatrix(buffmat,4*Lrealdim,4*Rrealdim,t(i,l),'N')
						LRcoeffout(:,:,j)=buffmat+LRcoeffout(:,:,j)
					end do
				end do
			end if
		end do
		
		! confirm that the pppVmat and hopmat can be used again without problem
		if(myid==orbid(i) .and. ifhop==.true.) then
			do j=1,pppVntouched,1
				call MPI_WAIT(pppVsendrequest(j),status,ierr)
			end do
			do j=1,hopntouched,1
				call MPI_WAIT(hopsendrequest(j),status,ierr)
			end do
		end if
	end do
	end if   ! myid/=0 ends

	! every process transfer LRcoeffout to coeffnosymm
	if(allocated(LRcoeffout)) then
		m=0
		do l=1,smadim,1
		do i=1,4*Rrealdim,1
		do j=1,4*Lrealdim,1
			if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
				quantabigL(j,2)+quantabigR(i,2)==totalSz) then
				m=m+1
				coeffnosymm(m)=LRcoeffout(j,i,l)
			end if
		end do
		end do
		end do
		if(m/=ngoodstates*smadim) then
			write(*,*) "--------------------myid",myid
			write(*,*) "m/ngoodstates*smadim",m
			stop
		end if
	else
		coeffnosymm=0.0D0   ! other process coeffnosymm does not sum up
	end if
	
	if(allocated(LRcoeffin)) deallocate(LRcoeffin)
	if(allocated(LRcoeffout)) deallocate(LRcoeffout)
	if(allocated(buffmat)) deallocate(buffmat)
	if(allocated(pppVmat)) deallocate(pppVmat)
	if(allocated(hopmat)) deallocate(hopmat)
	
	if(myid==0) then
		allocate(coeffnosymmreduce(ngoodstates*smadim),stat=error)
		if(error/=0) stop
	end if

	call MPI_REDUCE(coeffnosymm,coeffnosymmreduce,ngoodstates*smadim,mpi_real8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	
	if(myid==0) then
		newcoeff=0.0D0
		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
			do i=1,smadim,1
				call symmetrizestate(ngoodstates,coeffnosymmreduce(ngoodstates*(i-1)+1:i*ngoodstates),&
					newcoeff(bigdim*(i-1)+1:i*bigdim),'s')
			end do
		else
			newcoeff=coeffnosymmreduce
		end if
	end if
	
	if(allocated(coeffnosymm)) deallocate(coeffnosymm)
	if(allocated(coeffnosymmreduce)) deallocate(coeffnosymmreduce)

return

end subroutine op





