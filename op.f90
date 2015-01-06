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

	implicit none

	integer :: error,i,j,k,l,m,i1,j1,m1,l1
	integer :: bigdim,smadim
	! bigdim is the totaldim 16M*M,smadim is the davidson small block dimension
	! bigdim may be < 16M*M because we use spin reversal symmetry, and the dimension is smaller may be half
	! if groud state smadim=1
	! if gs+ex smadim may be >1
	real(kind=8) :: coeff(bigdim*smadim),newcoeff(bigdim*smadim)
	! coeff is the input coefficient and in the 1-d arrary format
	! new coeff is H cross C result
	integer :: operaindex
	integer :: status(MPI_STATUS_SIZE)
	real(kind=8),allocatable :: LRcoeff(:,:,:)
	real(kind=8),allocatable :: componentmat(:,:,:,:),buffermat(:,:)
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
		
		if(bigdim/=ngoodstates) then
			write(*,*) "------------------------------------"
			write(*,*) "op bigdim/=ngoodstates wrong!"
			write(*,*) "------------------------------------"
			stop
		end if

		!if(logic_spinreversal/=0) then
		!	write(*,*) "spinreversal symmetry. first do S*C, then do H*(SC), at last S(+)*(HSC)"
			! call S*C subroutine
		!else
		!	write(*,*) "do H*C"
		!end if
!-----------------------------------------------------------------------------
!to transform the 16M*M coeff to 4M*4M(L*R) format coeff(16M^2,n) to coeff(4M,4M,n) 
! since the input coeff is ngoodstates and other nongoodstates sets to 0
		m=1
		LRcoeff=0.0D0
		do k=1,smadim,1
			do i=1,4*Rrealdim,1
			do j=1,4*Lrealdim,1
				if((quantabigL(j,1)+quantabigR(i,1)==nelecs+ncharges) .and. &
					quantabigL(j,2)+quantabigR(i,2)==totalSz) then
					LRcoeff(j,i,k)=coeff(m)
					m=m+1
				end if
			end do
			end do
		end do
		m=m-1
		if(m/=smadim*ngoodstates) then
			write(*,*) "------------------------------------"
			write(*,*) "op good quantum states number wrong!"
			write(*,*) "------------------------------------"
			stop
		end if
	end if


	call MPI_bcast(LRcoeff,16*Lrealdim*Rrealdim*smadim,MPI_real8,0,MPI_COMM_WORLD,ierr)
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
				if((quantabigL(j,1)+quantabigR(i,1)==nelecs+ncharges) .and. &
					quantabigL(j,2)+quantabigR(i,2)==totalSz) then
					newcoeff(m)=componentmat(j,i,6,l)+newcoeff(m)
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

		do l=1,nleft+1,1
			if(myid==orbid(l)) then
				if(mod(l,nprocs-1)==0) then
					operaindex=l/(nprocs-1)
				else
					operaindex=l/(nprocs-1)+1
				end if
				do j=1,smadim,1
					if(bondlink(i,l)==1) then
						! the +1 -1 phase added to l'
						do m=1,4*Lrealdim,1
							componentmat(m,:,1:2,j)=componentmat(m,:,1:2,j)*((-1.0D0)**(mod(quantabigL(m,1),2)))
							!transfer from al*ar^(+) to ar^(+)*al
							componentmat(m,:,4:5,j)=componentmat(m,:,4:5,j)*((-1.0D0)**(mod(quantabigL(m,1),2)))*(-1.0D0)
						end do
						
						componentmat(:,:,6,:)=0.0D0
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
							if(k/=3) then
								componentmat(:,:,6,j)=buffermat*t(i,l)+componentmat(:,:,6,j)
							else 
								componentmat(:,:,6,j)=buffermat*pppV(i,l)+componentmat(:,:,6,j)
							end if
						end do
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
				do i1=1,4*Rrealdim,1
				do j1=1,4*Lrealdim,1
					if((quantabigL(j1,1)+quantabigR(i1,1)==nelecs+ncharges) .and. &
						quantabigL(j1,2)+quantabigR(i1,2)==totalSz) then
						newcoeff(m1)=componentmat(j1,i1,6,l1)+newcoeff(m1)
						m1=m1+1
					end if
				end do
				end do
				end do
			end if
		end do
	end do
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


deallocate(componentmat)
deallocate(LRcoeff)
if(myid/=0) then
	deallocate(buffermat)
end if

return

end subroutine op





