subroutine splitsvdL(singularvalue,leftu,statebegin,stateend,indexlp1)
! this subroutine is used to split the reduced density matrix
! to different subspace according to good quantum number
! and diagonalizaiton it to get the renormalized vector
	use mpi
	use variables
	USE blas95
	use lapack95
	USE F95_PRECISION
    use selectstate

	implicit none

	real(kind=8) :: leftu(4*Lrealdim,subM),singularvalue(subM),discard
	integer :: statebegin,stateend,indexlp1
	! statebegin is the index of the begin state
	! stateend is the index of the end state
	real(kind=8),allocatable :: valuework(:),coeffwork(:,:),coeffbuffer(:,:)&
	,coeffresult(:,:),transform(:,:),coeffdummy(:,:)
	! coeffbuffer is used to store the initial coeffwork
	! coeffresult stores the output of the eigenvector of the density matrix
	! transform stores the tranfer matrix when Sz==0
	integer,allocatable :: valueindex(:),szzeroindex(:),quantabigLbuffer(:,:),symmlinkbigbuffer(:)
	integer :: i,error,j,k,l,m,n,info,p,q,l1
	integer :: szzero,szl0,himp1
	logical :: done
	integer,allocatable :: subspacenum(:)
	real(kind=8) :: sum1
! szl0 means the number of sz>0 states
! szzero means the number of sz=0 state

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
    if(exscheme == 4 .and. startedMaxOverlap) then ! He Ma  
        call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,targettedStateIndex),&
                  coeffIF(1:4*Lrealdim,1:4*Rrealdim,targettedStateIndex),&
		          coeffwork(1:4*Lrealdim,1:4*Lrealdim),'N','T',1.0D0,0.0D0)
    else
	    do i=statebegin,stateend,1
		    call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),&
		    coeffbuffer(1:4*Lrealdim,1:4*Lrealdim),'N','T',1.0D0,0.0D0)
		    if(exscheme==1 .or. (exscheme == 4 .and. startedMaxOverlap==.false.)) then
			    coeffwork=coeffwork+coeffbuffer*nweight(i)
		    else
			    coeffwork=coeffwork+coeffbuffer
		    end if
        end do
    end if
	
	sum1=0.0D0
	do i=1,4*Lrealdim,1
	sum1=sum1+coeffwork(i,i)
	end do
	
	write(*,*) "the norm of the coeff"
	write(*,*) "sum1=",sum1
		
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
	allocate(symmlinkbigbuffer(4*subM),stat=error)
	if(error/=0) stop
	symmlinkbigbuffer=0
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

		do i=0,2*(nleft+1),1
			
			m=0
		!	coeffwork=coeffbuffer
			do k=1,4*subM,1
				call copy(coeffbuffer(:,k),coeffwork(:,k))
			end do
			
			if(logic_spinreversal/=0) then
				szzeroindex=0
			end if

			do k=1,4*Lrealdim,1
				if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
					m=m+1
					call swap(coeffwork(:,m),coeffwork(:,k))
					call swap(coeffwork(m,:),coeffwork(k,:))
					if(logic_spinreversal/=0 .and. j==0) then
						szzeroindex(m)=k
					end if
				end if
			end do

			if(m/=0) then
			!	write(*,*) i,j,m
			!	write(*,*) coeffwork(1:m,1:m)
			! when j==0 we first transform the basis to the new basis
			! which the symmlink is him self
				if(logic_spinreversal/=0 .and. j==0 ) then
					allocate(transform(m,m),stat=error)
					if(error/=0) stop
					transform=0.0D0
				
					do l=1,m,1
						if(abs(symmlinkbig(szzeroindex(l),1,1))==szzeroindex(l)) then
							transform(l,l)=1.0D0
						else
							do q=l+1,m,1
								if(symmlinkbig(szzeroindex(l),1,1)==szzeroindex(q)) then
								! the fisrt column the parity is 1
								! the second column the parity is -1
									transform(l,l)=sqrt(2.0D0)/2.0D0
									transform(q,l)=sqrt(2.0D0)/2.0D0
									transform(l,q)=sqrt(2.0D0)/2.0D0
									transform(q,q)=-sqrt(2.0D0)/2.0D0
								else if(symmlinkbig(szzeroindex(l),1,1)==-szzeroindex(q)) then
									transform(l,l)=sqrt(2.0D0)/2.0D0
									transform(q,l)=-sqrt(2.0D0)/2.0D0
									transform(l,q)=sqrt(2.0D0)/2.0D0
									transform(q,q)=sqrt(2.0D0)/2.0D0
								end if
							end do
						end if
					end do
					allocate(coeffdummy(m,m),stat=error)
					if(error/=0) stop
					call gemm(coeffwork(1:m,1:m),transform,coeffdummy(1:m,1:m),'N','N',1.0D0,0.0D0)
					call gemm(transform,coeffdummy(1:m,1:m),coeffwork(1:m,1:m),'T','N',1.0D0,0.0D0)
				! swap the symmlink==1 to the first few columns 
				! and the symmlink==-1 to the last few columns
					himp1=0
					do l=1,m,1
						if(symmlinkbig(szzeroindex(l),1,1)==szzeroindex(l)) then
							himp1=himp1+1
							call swap(coeffwork(:,himp1),coeffwork(:,l))
							call swap(coeffwork(himp1,:),coeffwork(l,:))
							call swap(transform(:,l),transform(:,himp1))
						else if(symmlinkbig(szzeroindex(l),1,1)==-szzeroindex(l)) then
							cycle
						else
							done=.true.
							do l1=1,l-1,1
								if(abs(symmlinkbig(szzeroindex(l),1,1))==szzeroindex(l1)) then
									done=.false.
									exit
								end if
							end do
							if(done==.true.) then
								himp1=himp1+1
								call swap(coeffwork(:,himp1),coeffwork(:,l))
								call swap(coeffwork(himp1,:),coeffwork(l,:))
								call swap(transform(:,l),transform(:,himp1))
							end if
						end if
					end do
					call syevd(coeffwork(1:himp1,1:himp1),valuework(n+1:n+himp1),'V','U',info)
					if(info/=0) then
						write(*,*) "left diagnolization failed! himp1"
						stop
					end if
					symmlinkbigbuffer(n+1:n+himp1)=1
					! symmlinkbigbuffer==1 means the symmlink is himself
					! symmlinkbigbuffer==-1 means the symmlink is -himself
					call syevd(coeffwork(himp1+1:m,himp1+1:m),valuework(n+himp1+1:n+m),'V','U',info)
					if(info/=0) then
						write(*,*) "left diagnolization failed! himm1"
						stop
					end if
					symmlinkbigbuffer(n+himp1+1:n+m)=-1
						
					call gemm(transform,coeffwork(1:m,1:m),coeffdummy(1:m,1:m),'N','N',1.0D0,0.0D0)
					coeffwork(1:m,1:m)=coeffdummy(1:m,1:m)
					
					deallocate(transform)
					deallocate(coeffdummy)
				else
					call syevd(coeffwork(1:m,1:m),valuework(n+1:n+m),'V','U',info)
				!	write(*,*) valuework(n+1:n+m)
					if(info/=0) then
						write(*,*) "left diagnolization failed!"
						stop
					end if
				end if

				p=0
				do k=1,4*Lrealdim,1
					if(quantabigL(k,1)==i .and. quantabigL(k,2)==j) then
						p=p+1
						call copy(coeffwork(p,1:m),coeffresult(k,n+1:n+m))
					end if
				end do

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
		else if(j==0 .and. logic_spinreversal/=0) then
			szzero=n-szl0
		end if

	end do
	
	if(logic_spinreversal/=0) then
		write(*,*) "szl0=",szl0,"szzero=",szzero
		if((szl0*2+szzero)/=4*Lrealdim) then
			write(*,*) "(szl0*2+szzero)/=4*Lrealdim failed!","szl0",szl0,"szzero",szzero
			stop
		end if
		if(sum(subspacenum(2:subspacenum(1)+1))/=szl0+szzero) then
			write(*,*) "------------------------"
			write(*,*) "subspacenum=szl0+szzero,failed!",subspacenum
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
				coeffresult(abs(symmlinkbig(i,1,1)),szl0+szzero+1:4*Lrealdim)=coeffresult(i,1:szl0)*DBLE(sign(1,symmlinkbig(i,1,1)))
			end if
		end do
		quantabigLbuffer(szl0+szzero+1:4*Lrealdim,1)=quantabigLbuffer(1:szl0,1)
		quantabigLbuffer(szl0+szzero+1:4*Lrealdim,2)=-1*quantabigLbuffer(1:szl0,2)
		valuework(szl0+szzero+1:4*Lrealdim)=valuework(1:szl0)
	end if

	do i=1,4*subM,1
		if(valuework(i)<-1.0D-10) then
			write(*,*) "-----------------------------"
			write(*,*) "caution valuework<0.0D0",valuework(i)
			write(*,*) "-----------------------------"
			valuework(i)=0.0D0
		end if
	end do

	if(logic_spinreversal==0) then
		call selectstates(valuework,4*Lrealdim,valueindex,singularvalue,subspacenum,nleft)
	else
		call selectstates(valuework,4*Lrealdim,valueindex,singularvalue,subspacenum,nleft,szzero,szl0)
	end if
! check if the valueindex is right
	do i=1,subM,1
		if(valueindex(i)==0) then
			write(*,*) "----------------------------------"
			write(*,*) "splitsvdL valueindex(i)==0",i
			write(*,*) "----------------------------------"
			stop
		end if
	end do

	do i=1,subM,1
		do j=i+1,subM,1
			if(valueindex(i)==valueindex(j)) then
				write(*,*) "----------------------------------"
				write(*,*) "splitsvdL valueindex(i)=valueindex(j)",i,j,valueindex(i)
				write(*,*) "----------------------------------"
				stop
			end if
		end do
	end do
		
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
			if(valueindex(i)<=szl0) then
				symmlinksma(i,1,1)=i+1
			else if(valueindex(i)>szl0 .and. valueindex(i)<=szl0+szzero) then
				if(symmlinkbigbuffer(valueindex(i))==0) then
					write(*,*) "----------------------------------------"
					write(*,*) "L symmlinkbigbuffer(valueindex(i)==0 failed!"
					write(*,*) "----------------------------------------"
					stop
				end if
				symmlinksma(i,1,1)=i*symmlinkbigbuffer(valueindex(i))
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
	deallocate(symmlinkbigbuffer)
end if
return
end subroutine

