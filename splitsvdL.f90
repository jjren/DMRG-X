subroutine splitsvdL(singularvalue,leftu,statebegin,stateend,indexlp1)
! this subroutine is used to split the reduced density matrix
! to different subspace according to good quantum number
! and diagonalizaiton it to get the renormalized vector
	use mpi
	use variables
	USE blas95
	use lapack95
	USE F95_PRECISION

	real(kind=8) :: leftu(4*Lrealdim,subM),singularvalue(subM)
	integer :: statebegin,stateend,indexlp1
	real(kind=8),allocatable :: valuework(:),coeffwork(:,:),coeffbuffer(:,:)&
	,quantabigLbuffer(:,:)
	integer,allocatable :: valueindex(:)
	integer :: i,error,j,k,l,m,n,info,i1,symm,p
	integer :: szzero,szl0,szzero1
	logical :: done
! szzero meane the number of sz=0 state
! szl0 means the number of sz>0 state
! szzero1 means sz=0 and spin reversal parity=1


	write(*,*) "enter in subroutine splitsvdL"

	allocate(coeffbuffer(4*subM,4*subM),stat=error)
	if(error/=0) stop
	allocate(coeffwork(4*subM,4*subM),stat=error)
	if(error/=0) stop
	
	coeffbuffer=0.0D0


! the L+sigmaL space reduced density matrix
	do i=statebegin,stateend,1
		call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),&
		coeffbuffer(1:4*Lreadlim,1:4*Lrealdim),'N','T',1.0D0,0.0D0)
		if(exscheme==1) then
			coeffwork=coeffwork+coeffbuffer*nweight(i)
		else
			coeffwork=coeffwork+coeffbuffer
		end if
	end do
	deallocate(coeffbuffer)
! split the reduced density matrix to different good quantum
! number subspace
	
	allocate(valuework(4*subM),stat=error)
	if(error/=0) stop
	allocate(valueindex(subM),stat=error)
	if(error/=0) stop
	allocate(quantabigLbuffer(4*subM),stat=error)
	if(error/=0) stop
	
	m=0
! m is the number of good quantum state so far
! j is Sz
! i is nelecs
	do j=nelecs,-nelecs,-2
		
		if(logic_spinreversal/=0 .and. j<0) then
			do i=1,4*Lrealdim,1
				if(quantabigL(i,2)>0) then
					coeffwork(i,1:szl0)=coeffwork(abs(symmlinkbig(i,1,1)),szl0+szzero+1:4*Lrealdim)*sign(1.0D0,symmlinkbig(i,1,1))
				end if
			end do
				quantabigLbuffer(szl0+szzero+1:4*Lrealdim,:)=quantabigLbuffer(1:szl0,:)
			exit
		end if

! when Sz=0 we need to seperate spin reversal=1 or -1
			if(logic_spinreversal/=0 .and. j==0) then
				symm=2
			else
				symm=1
			end if
				

		do i1=1,symm,1
		
		do i=0,nelecs,1

		
			n=m
! n is the last total m
			
			if(symm==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			else if(symm==2 .and. i1==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==k) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			else if(symm==2 .and. i1==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==-k) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			end if

			if(m/=n) then
				call syevd(coeffwork(n+1:m,n+1:m),valuework(n:m),'V','U',info)
				if(info/=0) then
					write(*,*) "left diagnolization failed!"
					stop
				end if
				coeffwork(m+1:4*subM,n+1:m)=0.0D0
				coeffwork(1:n,n+1:m)=0.0D0
				
				p=1
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
						call copy(coeffwork(p,n+1:m),coeffwork(k,n+1:m))
						coeffwork(p,n+1:m)=0.0D0
						p=p+1
					end if
				end do
				if(p-1/=m-n) then
					write(*,*) "p-1/=m-n failed!"
					stop
				end if
				quantabigLbuffer(n+1:m,1)=i
				quantabigLbuffer(n+1:m,2)=j
			end if
		end do

		if(j==0 .and. i1==1 .and. logic_spinreversal/=0) then
			szzero1=m-szl0
		end if

		end do

		if(j>0) then
			szl0=m
		else if(j==0) then
			szzero=m-szl0
		end if
	end do
	
	if(logic_spinreversal/=0) then
	if((szl0*2+szzero)/=4*Lrealdim) then
		write(*,*) "(szl0*2+szzero)/=4*Lrealdim failed!"
		stop
	end if
	end if
	
	singularvalue=0.0D0
	valueindex=0
	if(logic_spinreversal==0) then
		do i=1,4*Lrealdim,1
			do j=1,subM,1
				if(valuework(i)>singularvalue(j)) then
					valueindex(j+1:subM)=valueindex(j:subM-1)
					singularvalue(j+1:subM)=singularvalue(j:subM-1)
					valueindex(j)=i
					singularvalue(j)=valuework(i)
					exit
				end if
			end do
		end do
	else
		do i=1,szzero+szl0,1
			do j=1,subM,1
				if(valuework(i)>singularvalue(j)) then
					if(i<=szl0) then
						valueindex(j+2:subM)=valueindex(j:subM-2)
						valueindex(j)=i
						valueindex(j+1)=szl0+szzero+i
						singularvalue(j+2:subM)=singularvalue(j:subM-2)
						singularvalue(j)=valuework(i)
						singularvalue(j+1)=valuework(i)
					else
						valueindex(j+1:subM)=valueindex(j:subM-1)
						singularvalue(j+1:subM)=singularvalue(j:subM-1)
						valueindex(j)=i
						singularvalue(j)=valuework(i)
					end if
				end if
			end do
		end do
		if(quantabigLbuffer(valueindex(subM),2)>0) then
			do i=szl0+szzero,szl0+1,-1
			do j=1,subM,1
				if(valueindex(j)==i) then
					done=.false.
					exit
				else
					done=.true.
				end if
			end do
				if(done==.true.) then
				valueindex(subM)=i
				singularvalue(subM)=valuework(i)
				exit
				end if
			end do
		end if
	end if
	
	do i=1,subM,1
		leftu(:,i)=coeffwork(1:4*Lrealdim,valueindex(i))
		quantasmaL(i,1)=quantabigLbuffer(valueindex(i),1)
		quantasmaL(i,2)=quantabigLbuffer(valueindex(i),2)
		if(logic_spinreversal/=0) then
			if(valueindex(i)<=szl0) then
				symmlinksma(i,1,1)=i+1
			else if(valueindex(i)>szl0 .and. valueindex(i)<=szl0+szzero) then
				if(valueindex(i)<=szl0+szzero1) then
					symmlinksma(i,1,1)=i
				else
					symmlinksma(i,1,1)=-i
				end if
			else 
				symmlinksma(i,1,1)=i-1
			end if
		end if
	end do


!   the total discard weight between site indexLp1,indexRm1
		discard=1.0D0-sum(singularvalue(1:subM))
		write(*,'(A20,I4,D12.5)') "totaldiscardL=",indexLp1,discard

deallocate(valuework)
deallocate(coeffwork)
deallocate(quantabigLbuffer)
deallocate(valueindex)

return
end subroutine

