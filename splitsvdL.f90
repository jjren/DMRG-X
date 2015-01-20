subroutine splitsvdL(singularvalue,leftu,statebegin,stateend,indexlp1)
! this subroutine is used to split the reduced density matrix
! to different subspace according to good quantum number
! and diagonalizaiton it to get the renormalized vector
	use mpi
	use variables
	USE blas95
	use lapack95
	USE F95_PRECISION

	implicit none

	real(kind=8) :: leftu(4*Lrealdim,subM),singularvalue(subM),discard
	integer :: statebegin,stateend,indexlp1
	real(kind=8),allocatable :: valuework(:),coeffwork(:,:),coeffbuffer(:,:)&
	,coeffresult(:,:)
	integer,allocatable :: valueindex(:),szzeroindex(:),quantabigLbuffer(:,:)
	integer :: i,error,j,k,l,m,n,info,i1,symm1,p,p1
	integer :: szzero,pair1,szzero1,szl0
	logical :: done,done2
!	integer :: tmp
	integer,allocatable :: subspacenum(:)
! szl0 means the number of sz>0 states
! szzero means the number of sz=0 state including spin parity=+-1
! pair1 means the number of sz>0 and sz=0 but the symmetry pair is not himself state
! szzero1 means the number sz=0 and spin reversal parity=1


	write(*,*) "enter in subroutine splitsvdL"

	allocate(coeffbuffer(4*subM,4*subM),stat=error)
	if(error/=0) stop
	allocate(coeffwork(4*subM,4*subM),stat=error)
	if(error/=0) stop
	allocate(coeffresult(4*subM,4*subM),stat=error)
	if(error/=0) stop
	
	coeffwork=0.0D0
	coeffbuffer=0.0D0
	coeffresult=0.0D0

!	write(*,*) "coeffIF"
!	write(*,*) coeffIF(1:4*Lrealdim,1:4*Rrealdim,:)

! the L+sigmaL space reduced density matrix
	do i=statebegin,stateend,1
		call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),&
		coeffbuffer(1:4*Lrealdim,1:4*Lrealdim),'N','T',1.0D0,0.0D0)
		if(exscheme==1) then
			coeffwork=coeffwork+coeffbuffer*nweight(i)
		else
			coeffwork=coeffwork+coeffbuffer
		end if
	end do
! debug
!	write(*,*) "coeffwork"
!	do i=1,4*Lrealdim,1
!	do j=1,4*Lrealdim,1
!	if(abs(coeffwork(j,i))>1.0D-8) then
!	write(*,*) coeffwork(j,i),j,i
!	end if
!	end do
!	end do
! split the reduced density matrix to different good quantum
! number subspace
	
	allocate(valuework(4*subM),stat=error)
	if(error/=0) stop
	valuework=0.0D0
	allocate(valueindex(subM),stat=error)
	if(error/=0) stop
	allocate(quantabigLbuffer(4*subM,2),stat=error)
	if(error/=0) stop
	allocate(subspacenum((2*(nleft+1)+1)**2+1),stat=error)
	if(error/=0) stop
	subspacenum=0
	
	if(logic_spinreversal/=0) then
	allocate(szzeroindex(4*subM),stat=error)
	if(error/=0) stop
	end if
	
	coeffbuffer=coeffwork
! m is the number of good quantum state so far
! j is Sz
! i is nelecs
	n=0
	! n is the total number of states
	do j=nleft+1,-nleft-1,-1
		
		if(logic_spinreversal/=0 .and. j<0) then
			exit
		! using pair symmetry copy
		end if

! when Sz=0 we need to seperate spin reversal=1 or -1 and symmetry pair is not
! himself
		if(logic_spinreversal/=0 .and. j==0) then
			symm1=3
		else
			symm1=1
		end if
				

		do i1=1,symm1,1
		
		do i=0,2*(nleft+1),1

		
			m=0
			coeffwork=coeffbuffer
! n is the last total m
			if(symm1==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			!	write(*,*) m
			else if(symm1==3 .and. i1==1) then
				szzeroindex=0
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. abs(symmlinkbig(k,1,1))/=k) then
						done2=.true.
						do p1=1,m,1
							if(abs(symmlinkbig(k,1,1))/=szzeroindex(p1)) then
								done2=.true.
							else
								done2=.false.
								exit
							end if
						end do
						if(done2==.true.) then
							m=m+1
							szzeroindex(m)=k
							call swap(coeffwork(:,m),coeffwork(:,k))
							call swap(coeffwork(m,:),coeffwork(k,:))
						end if
					end if
				end do
			!	write(*,*) m
			else if(symm1==3 .and. i1==2) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==k) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			!	write(*,*) m
			else if(symm1==3 .and. i1==3) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==-k) then
						m=m+1
						call swap(coeffwork(:,m),coeffwork(:,k))
						call swap(coeffwork(m,:),coeffwork(k,:))
					end if
				end do
			!	write(*,*) m
			end if

			if(m/=0) then
			!	write(*,*) i,j,i1,m
			!	write(*,*) coeffwork(1:m,1:m)
				call syevd(coeffwork(1:m,1:m),valuework(n+1:n+m),'V','U',info)
				if(info/=0) then
					write(*,*) "left diagnolization failed!"
					stop
				end if
			!	write(*,*) i,j,m
			!	write(*,*) valuework(n+1:n+m)
				
			p=0
			if(symm1==1) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
						p=p+1
						call copy(coeffwork(p,1:m),coeffresult(k,n+1:n+m))
					end if
				end do
			else if(symm1==3 .and. i1==1) then
					do p1=1,4*Lrealdim,1
						if(szzeroindex(p1)/=0) then
							p=p+1
							call copy(coeffwork(p,1:m),coeffresult(szzeroindex(p1),n+1:m))
						else
							exit
						end if
					end do
			else if(symm1==3 .and. i1==2) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==k) then
						p=p+1
						call copy(coeffwork(p,1:m),coeffresult(k,n+1:m))
					end if
				end do
			else if(symm1==3 .and. i1==3) then
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j .and. symmlinkbig(k,1,1)==-k) then
						p=p+1
						call copy(coeffwork(p,1:m),coeffresult(k,n+1:m))
					end if
				end do
			end if

				if(p/=m) then
					write(*,*) "p/=m failed!",p,m
					stop
				end if
				quantabigLbuffer(n+1:n+m,1)=i
				quantabigLbuffer(n+1:n+m,2)=j

				subspacenum(1)=subspacenum(1)+1
				subspacenum(subspacenum(1)+1)=m
			end if
		n=m+n

		end do

		if(j>0) then
			szl0=n
		else if(j==0 .and. i1==1 .and. logic_spinreversal/=0) then
			pair1=n
		else if(j==0 .and. i1==2 .and. logic_spinreversal/=0) then
			szzero1=n-pair1
		else if(j==0 .and. i1==3 .and. logic_spinreversal/=0) then
			szzero=n-pair1
		end if

		end do

	end do
	
	if(logic_spinreversal/=0) then
		write(*,*) "pair1=",pair1,"szl0=",szl0,"szzero1=",szzero1,"szzero=",szzero
		if((pair1*2+szzero)/=4*Lrealdim) then
			write(*,*) "(pair1*2+szzero)/=4*Lrealdim failed!","pair1=",pair1,"szzero",szzero
			stop
		end if
		if(sum(subspacenum(2:subspacenum(1)+1))/=pair1+szzero) then
			write(*,*) "------------------------"
			write(*,*) "subspacenum=pair1+szzero,failed!",subspacenum
			write(*,*) "------------------------"
			stop
		end if
	else
		if(n/=4*Lrealdim) then
			write(*,*) "------------------------"
			write(*,*) "n/=4*Lrealdim,failed!",n
			write(*,*) "------------------------"
			stop
		end if
		if(sum(subspacenum(2:subspacenum(1)+1))/=4*Lrealdim) then
			write(*,*) "------------------------"
			write(*,*) "subspacenum=4*Lrealdim,failed!",subspacenum
			write(*,*) "------------------------"
			stop
		end if
	end if
	


! copy the Sz>0 part and Sz=0(symmetry pair is not himself) to the symmetry pair
	if(logic_spinreversal/=0) then
		do i=1,4*Lrealdim,1
			if(quantabigL(i,2)>0) then
				coeffresult(abs(symmlinkbig(i,1,1)),szl0+szzero+1:2*szl0+szzero)=coeffresult(i,1:szl0)*DBLE(sign(1,symmlinkbig(i,1,1)))
			end if
		end do
			do p1=1,4*Lrealdim,1
				if(szzeroindex(p1)/=0) then
				coeffresult(abs(symmlinkbig(szzeroindex(p1),1,1)),pair1+szzero+szl0+1:4*Lrealdim)&
				=coeffresult(szzeroindex(p1),szl0+1:pair1)*&
					DBLE(sign(1,symmlinkbig(szzeroindex(p1),1,1)))
				end if
			end do
		quantabigLbuffer(pair1+szzero+1:4*Lrealdim,1)=quantabigLbuffer(1:pair1,1)
		quantabigLbuffer(pair1+szzero+1:4*Lrealdim,2)=-1*quantabigLbuffer(1:pair1,2)
	end if


	if(logic_spinreversal==0) then
		call selectstates(valuework,4*Lrealdim,valueindex,singularvalue,subspacenum,nleft)
	else
		call selectstates(valuework,4*Lrealdim,valueindex,singularvalue,subspacenum,nleft,szzero,pair1)
	end if
! when 4*Lrealdim<subM (in the debug mode)
	if(4*Lrealdim>subM) then
		m=subM
	else
		m=4*Lrealdim
	end if

	do i=1,m,1
		leftu(:,i)=coeffresult(1:4*Lrealdim,valueindex(i))
		quantasmaL(i,:)=quantabigLbuffer(valueindex(i),:)
		if(logic_spinreversal/=0) then
			if(valueindex(i)<=pair1) then
				symmlinksma(i,1,1)=i+1
			else if(valueindex(i)>pair1 .and. valueindex(i)<=pair1+szzero) then
				if(valueindex(i)<=pair1+szzero1) then
					symmlinksma(i,1,1)=i
				else
					symmlinksma(i,1,1)=-i
				end if
			else 
				symmlinksma(i,1,1)=i-1
			end if
		end if
	end do
	!	write(*,*) "valuework",valuework
	!	write(*,*) "singularvalue",singularvalue
	!	write(*,*) "valueindex",valueindex
!   the total discard weight between site indexLp1,indexRm1
		discard=1.0D0-sum(singularvalue(1:m))
		write(*,'(A20,I4,D12.5)') "totaldiscardL=",indexLp1,discard


deallocate(valuework)
deallocate(coeffwork)
deallocate(quantabigLbuffer)
deallocate(valueindex)
deallocate(coeffbuffer)
deallocate(coeffresult)
deallocate(subspacenum)
if(logic_spinreversal/=0) then
	deallocate(szzeroindex)
end if
!	tmp=0
!	do while(tmp==0) 
!		call sleep(2)
!	end do
return
end subroutine

