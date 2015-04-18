subroutine op(bigdim,smadim,coeff,newcoeff)
! this is the core subroutine to calculate the S*H*S*C or H*C
! the parallel schema follow JCP 12 3174(2004) garnet chan
! if want to save memory, then can write a wrapper, to send one coeff every time

! input bigdim,smadim,coeff
! output newcoeff

	use mpi
	use variables
	use BLAS95
	use F95_PRECISION
	use symmetry
	use communicate

	implicit none

	integer :: error,i,j,k,l,m,ii,j1,m1,l1,ierr
	integer :: bigdim,smadim
	! bigdim is the totaldim 16M*M,smadim is the davidson small block dimension
	! bigdim may be < 16M*M because we use spin reversal symmetry, and the dimension is smaller may be half
	! if groud state smadim=1
	! if gs+ex smadim may be >1
	real(kind=8) :: coeff(bigdim*smadim),coeffnosymm(ngoodstates*smadim),newcoeff(bigdim*smadim)
	! coeff is the input coefficient and in the 1-d arrary format
	! new coeff is H cross C result
	integer :: operaindex
	integer :: status(MPI_STATUS_SIZE)
	real(kind=8),allocatable :: LRcoeff(:,:,:)
	real(kind=8),allocatable :: componentmat(:,:,:,:),buffermat(:,:)
	logical :: done
	real(kind=8) :: norm
! debug
!	logical ::alive
!	real(kind=8),allocatable :: Hbuffer(:,:),newcoeffdummy(:)
!	real(kind=8) :: diff
!	integer :: tmp
!
!   transform the 1-array to 4M*4M form 
	allocate(LRcoeff(4*Lrealdim,4*Rrealdim,smadim),stat=error) 
	! 256M if nstate=2 M=1000
	if(error/=0) stop
!     
	allocate(componentmat(4*Lrealdim,4*Rrealdim,6,smadim),stat=error)
	if(error/=0) stop
	! 768M if nstate=2

	if(myid/=0) then
	allocate(buffermat(4*Lrealdim,4*Rrealdim),stat=error)
	if(error/=0) stop
	end if
	
	if(myid==0) then
	newcoeff=0.0D0
	end if

	if( myid==0 ) then
	!	write(*,*) "enter H*C OP subroutine"
		
	!!	if(bigdim/=ngoodstates) then
	!		write(*,*) "------------------------------------"
	!		write(*,*) "op bigdim/=ngoodstates wrong!"
	!		write(*,*) "------------------------------------"
	!		stop
	!	end if

		if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
		do i=1,smadim,1
			call symmetrizestate(ngoodstates,coeffnosymm(ngoodstates*(i-1)+1:i*ngoodstates),&
				coeff(bigdim*(i-1)+1:i*bigdim),'u')
		end do
		else
			coeffnosymm=coeff
		end if

!-----------------------------------------------------------------------------
!to transform the 16M*M coeff to 4M*4M(L*R) format coeff(16M^2,n) to coeff(4M,4M,n) 
! since the input coeff is ngoodstates and other nongoodstates sets to 0
		m=1
		LRcoeff=0.0D0
		do k=1,smadim,1
			do i=1,4*Rrealdim,1
			do j=1,4*Lrealdim,1
				if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
					quantabigL(j,2)+quantabigR(i,2)==totalSz) then
					LRcoeff(j,i,k)=coeffnosymm(m)
					m=m+1
				end if
			end do
			end do
		end do
	!	m=m-1
	!	if(m/=smadim*ngoodstates) then
	!		write(*,*) "------------------------------------"
	!		write(*,*) "op good quantum states number wrong!"
	!		write(*,*) "------------------------------------"
	!		stop
	!	end if
	end if


	call MPI_bcast(LRcoeff,16*Lrealdim*Rrealdim*smadim,MPI_real8,0,MPI_COMM_WORLD,ierr)
	coeffnosymm=0.0D0
!  calculate HL*1 and 1*HR 
	if(myid==0) then 
		do k=1,2,1
			do i=1,smadim
			if (k==1) then ! HL*1
				call gemm(Hbig(1:4*Lrealdim,1:4*Lrealdim,1),LRcoeff(:,:,i),&
					componentmat(:,:,6,i),'N','N',1.0D0,0.0D0)
			else ! 1*HR
				call gemm(LRcoeff(:,:,i),Hbig(1:4*Rrealdim,1:4*Rrealdim,2),&
					componentmat(:,:,6,i),'N','N',1.0D0,0.0D0)
			end if
			end do
			

			m=1
			do l=1,smadim,1
			do i=1,4*Rrealdim,1
			do j=1,4*Lrealdim,1
				if((quantabigL(j,1)+quantabigR(i,1)==nelecs) .and. &
					quantabigL(j,2)+quantabigR(i,2)==totalSz) then
					coeffnosymm(m)=componentmat(j,i,6,l)+coeffnosymm(m)
					m=m+1
				end if
			end do
			end do
			end do
		end do
	end if

!------------------------------------------------
! vlr=Hlrl'r'*Cl'r'=sum(opt,l',r')=sum(Lopt,l') parity*Oll'*sum(Ropt,r') Orr'cl'r'
! the parallel schema is that 0 process bcast the coeff matrix to other process
! and 0 process gather the result

	do i=norbs,norbs-nright,-1
		if(myid==orbid(i)) then
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			
			do j=1,smadim,1
				do k=1,5,1
					if(k<=3) then
						! firstly needed to transfer from a(+) to a, so 'T' and 'T' is 'N', (ni-1)^+=(ni-1)
						! k=1 aup,k=2 a down,k=3 n,k=4 a+ up,k=5 a+ down;
						call gemm(LRcoeff(:,:,j),operamatbig(1:4*Rrealdim,1:4*Rrealdim,(operaindex-1)*3+k)&
						,componentmat(:,:,k,j),'N','N',1.0D0,0.0D0)
					else
						call gemm(LRcoeff(:,:,j),operamatbig(1:4*Rrealdim,1:4*Rrealdim,(operaindex-1)*3+k-3)&
						,componentmat(:,:,k,j),'N','T',1.0D0,0.0D0)
					end if
				end do
			end do
		end if

		call MPI_bcast(componentmat(:,:,1:5,:),16*Lrealdim*Rrealdim*smadim*5,MPI_real8,orbid(i),MPI_COMM_WORLD,ierr)
! the +1 -1 phase added to l'
! can only do once
		do m=1,4*Lrealdim,1
			componentmat(m,:,1:2,:)=componentmat(m,:,1:2,:)*((-1.0D0)**(mod(quantabigL(m,1),2)))
			!transfer from al*ar^(+) to ar^(+)*al
			componentmat(m,:,4:5,:)=componentmat(m,:,4:5,:)*((-1.0D0)**(mod(quantabigL(m,1),2)))*(-1.0D0)
		end do

		do l=1,nleft+1,1
			if(myid==orbid(l)) then
				if(mod(l,nprocs-1)==0) then
					operaindex=l/(nprocs-1)
				else
					operaindex=l/(nprocs-1)+1
				end if
				
				componentmat(:,:,6,:)=0.0D0
				do j=1,smadim,1
					if(bondlink(i,l)==1) then
						
!							open(unit=997,file="imme2.tmp",status="replace")
						do k=1,5,1
							!k<=3 al^+*ar,(nl-1)(nr-1),k>3 al*ar^(+) 
							if(k<=3) then
							call gemm(operamatbig(1:4*Lrealdim,1:4*Lrealdim,(operaindex-1)*3+k),componentmat(:,:,k,j)&
							,buffermat,'N','N',1.0D0,0.0D0)
							else
							call gemm(operamatbig(1:4*Lrealdim,1:4*Lrealdim,(operaindex-1)*3+k-3),componentmat(:,:,k,j)&
							,buffermat,'T','N',1.0D0,0.0D0)
							! buffermat is to save the intermediate matrix
							end if
! debug
!							write(997,*) buffermat

							if(k/=3) then
								componentmat(:,:,6,j)=buffermat*t(i,l)+componentmat(:,:,6,j)
							else 
								componentmat(:,:,6,j)=buffermat*pppV(i,l)+componentmat(:,:,6,j)
							end if
						end do
!							close(997)
						! add all the five operator ai^+up*ajup,ai^_down*ajdown (ni-1)(nj-1)... together
					else ! not bond linked only the pppV term
						call gemm(operamatbig(1:4*Lrealdim,1:4*Lrealdim,operaindex*3),componentmat(:,:,3,j)&
						,componentmat(:,:,6,j),'N','N',1.0D0,0.0D0)
						componentmat(:,:,6,j)=componentmat(:,:,6,j)*pppV(i,l)
					end if
				end do
				
				call MPI_SEND(componentmat(:,:,6,:),16*Lrealdim*Rrealdim*smadim,mpi_real8,0,l,MPI_COMM_WORLD,ierr)
			end if

			if(myid==0) then
				call MPI_RECV(componentmat(:,:,6,:),16*Lrealdim*Rrealdim*smadim,mpi_real8,orbid(l),l,MPI_COMM_WORLD,status,ierr)
				m1=1
				do l1=1,smadim,1
				do ii=1,4*Rrealdim,1
				do j1=1,4*Lrealdim,1
					if((quantabigL(j1,1)+quantabigR(ii,1)==nelecs) .and. &
						quantabigL(j1,2)+quantabigR(ii,2)==totalSz) then
						coeffnosymm(m1)=componentmat(j1,ii,6,l1)+coeffnosymm(m1)
						m1=m1+1
					end if
				end do
				end do
				end do
			end if
		end do
	end do
		if(myid==0) then
			newcoeff=0.0D0
			if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
				do j=1,smadim,1
					call symmetrizestate(ngoodstates,coeffnosymm(ngoodstates*(j-1)+1:j*ngoodstates),&
						newcoeff(bigdim*(j-1)+1:j*bigdim),'s')
				end do
			else
				newcoeff=coeffnosymm
			end if
!				open(unit=200,file="tmp2.tmp",status="old",position="append")
!				write(200,*) "=============================="
!				write(200,*) newcoeff
!				close(200)
		end if
! let the non good quantum coeff to be zero
	!if(myid==0) then
	!	do j=1,4*Rrealdim,1
	!		do k=1,4*Lrealdim,1
	!			if((quantabigL(k,1)+quantabigR(j,1)/=nelecs+ncharges) .or. &
	!				(quantabigL(k,2)+quantabigR(j,2)/=totalSz)) then
	!				newcoeff((j-1)*4*Lrealdim+k:16*Rrealdim*Lrealdim*smadim:16*Rrealdim*Lrealdim)=0.0D0
	!			end if
	!		end do
	!	end do
	!end if

!	if(myid==0 .and. logic_spinreversal/=0) then
!			!do i=1,ngoodstates,1
!			do i=1,ngoodstates,1
!				if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1) .and. &
!				abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(i,2)) then
!					if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
!					/=logic_spinreversal*((-1)**mod(nelecs,2))) then
!						newcoeff(i:smadim*ngoodstates:ngoodstates)=0.0D0
!					end if
!				else if(abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1) .or. &
!				abs(symmlinkbig(symmlinkgood(i,2),1,2))/=symmlinkgood(i,2)) then
!					done=.false.
!					do j=1,ngoodstates,1
!						if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(j,1) .and. &
!						abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(j,2)) then
!							newcoeff(j:smadim*ngoodstates:ngoodstates)=newcoeff(i:smadim*ngoodstates:ngoodstates)&
!							*DBLE(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2)))&
!							*DBLE(logic_spinreversal)*((-1.0D0)**mod(nelecs,2))
!							done=.true.
!							exit
!						end if
!					end do
!					if(done==.false.) then
!						write(*,*) "-----------------------------------------------------------------------------"
!						write(*,*) "in op did't find the state",i,symmlinkgood(i,1),symmlinkgood(i,2),"corrsponds",&
!						symmlinkbig(symmlinkgood(i,1),1,1),symmlinkbig(symmlinkgood(i,2),1,2)
!						write(*,*) "-----------------------------------------------------------------------------"
!						stop
!					end if
!				end if
!			end do
!			do i=1,smadim,1
!				norm=dot(newcoeff((i-1)*ngoodstates+1:i*ngoodstates),newcoeff((i-1)*ngoodstates+1:i*ngoodstates))
!				if(norm<1.0D-10) then
!					write(*,*) "--------------------------"
!					write(*,*) "in op norm is < 1.0D-10,caution!"
!					write(*,*) "--------------------------"
!				end if
!				newcoeff((i-1)*ngoodstates+1:i*ngoodstates)=newcoeff((i-1)*ngoodstates+1:i*ngoodstates)/sqrt(norm)
!			end do
			!	if(quantabigL(symmlinkgood(i,1),2)>=0 .and. abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1)) then
			!		done=.false.
			!	do j=1,ngoodstates,1
			!		if(symmlinkgood(j,1)==abs(symmlinkbig(symmlinkgood(i,1),1,1)) &
			!			.and. symmlinkgood(j,2)==abs(symmlinkbig(symmlinkgood(i,2),1,2))) then
			!			newcoeff(j:smadim*ngoodstates:ngoodstates)=newcoeff(i:smadim*ngoodstates:ngoodstates)&
			!			*DBLE(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2)))&
			!			*DBLE(logic_spinreversal)
			!			done=.true.
			!			exit
			!		end if
			!	end do
			!		if(done/=.true.) then
			!			write(*,*) "-------------------------------------------------"
			!			write(*,*) "initialrandomweight spin reversal adapted failed!"
			!			write(*,*) "-------------------------------------------------"
			!			stop
			!		end if
			!	else if(quantabigL(symmlinkgood(i,1),2)==0 .and. &
			!	abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1)) then
			!		if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))/=logic_spinreversal) then
			!			newcoeff(i:smadim*ngoodstates:ngoodstates)=0.0D0
			!		end if
			!	end if
			!end do
!	end if

! debug
! logic :: alive
!	if(myid==0) then
!		inquire(file="H.tmp",exist=alive)
!		if(alive) then
!			open(unit=998,file="H.tmp",status="old")
!		else
!			write(*,*) "no H.tmp"
!			stop
!		end if
!		allocate(Hbuffer(ngoodstates,ngoodstates),stat=error)
!		if(error/=0) stop
!		allocate(newcoeffdummy(ngoodstates*smadim),stat=error)
!		if(error/=0) stop
!		read(998,*) Hbuffer
!		do i=1,smadim,1
!		call gemv(Hbuffer,coeff((i-1)*ngoodstates+1:i*ngoodstates),newcoeffdummy((i-1)*ngoodstates+1:i*ngoodstates),1.0D0,0.0D0,'N')
!		end do
!		do i=1,ngoodstates*smadim,1
!			if(abs(newcoeffdummy(i)-newcoeff(i))>1.0D-4) then
!				write(*,*) i,newcoeffdummy(i),newcoeff(i)
!			end if
!		end do
!		read(*,*) tmp
!		close(998)
!		open(unit=996,file="LRcoeff.tmp",status="replace")
!		write(996,*) newcoeff
!		close(996)
!	end if
		





deallocate(componentmat)
deallocate(LRcoeff)
if(myid/=0) then
	deallocate(buffermat)
end if

return

end subroutine op





