subroutine splitsvdR(singularvalue,rightv,statebegin,stateend,indexRm1)
! this subroutine is used to split the reduced density matrix
! to different subspace according to good quantum number
! and diagonalizaiton it to get the renormalized vector
	use mpi
	use variables
	USE blas95
	use lapack95
	USE F95_PRECISION

	implicit none

	real(kind=8) :: rightv(subM,4*Rrealdim),singularvalue(subM),discard
	integer :: statebegin,stateend,indexRm1
	real(kind=8),allocatable :: valuework(:),coeffwork(:,:),coeffbuffer(:,:)&
	,coeffresult(:,:),transform(:,:),coeffdummy(:,:)
	integer,allocatable :: valueindex(:),szzeroindex(:),quantabigRbuffer(:,:),&
	symmlinkbigbuffer(:)
	integer :: i,error,j,k,l,m,n,info,symm1,p,l1,q
	integer :: szzero,szl0,himp1
	logical :: done,ifexist,iffind
	integer,allocatable :: subspacenum(:)
	real(kind=8) :: diffzero,scale1
	integer :: scalenum
! szl0 means the number of sz>0 states
! szzero means the number of sz=0 state 


	write(*,*) "enter in subroutine splitsvdR"

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
!	write(*,*) coeffIF(1:4*Rrealdim,1:4*Rrealdim,:)

! the R+sigmaR space reduced density matrix
	do i=statebegin,stateend,1
		call gemm(coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),coeffIF(1:4*Lrealdim,1:4*Rrealdim,i),&
		coeffbuffer(1:4*Rrealdim,1:4*Rrealdim),'T','N',1.0D0,0.0D0)
		if(exscheme==1) then
			coeffwork=coeffwork+coeffbuffer*nweight(i)
		else
			coeffwork=coeffwork+coeffbuffer
		end if
	end do
	!write(*,*) "coeffwork"
	!write(*,*) coeffwork
! split the reduced density matrix to different good quantum
! number subspace
	
	allocate(valuework(4*subM),stat=error)
	if(error/=0) stop
	valuework=0.0D0
	allocate(valueindex(subM),stat=error)
	if(error/=0) stop
	allocate(quantabigRbuffer(4*subM,2),stat=error)
	if(error/=0) stop
	allocate(subspacenum((2*(nright+1)+1)**2+1),stat=error)
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
	do j=nright+1,-1-nright,-1
		
		if(logic_spinreversal/=0 .and. j<0) then
			exit
		! using pair symmetry copy
		end if

		do i=0,2*(nright+1),1

			m=0
			!coeffwork=coeffbuffer
			do k=1,4*subM,1
				call copy(coeffbuffer(:,k),coeffwork(:,k))
			end do

			if(logic_spinreversal/=0) then
				szzeroindex=0
			end if

			do k=1,4*Rrealdim,1
				if(quantabigR(k,1)==i .and. quantabigR(k,2)==j) then
					m=m+1
					call swap(coeffwork(:,m),coeffwork(:,k))
					call swap(coeffwork(m,:),coeffwork(k,:))
					if(logic_spinreversal/=0 .and. j==0) then
						szzeroindex(m)=k
					end if
				end if
			end do

			if(m/=0) then
			!	irite(*,*) i,j,m
			!	write(*,*) coeffwork(1:m,1:m)
			! when j==0 we first transform the basis to the new basis
			! which the symmlink is him self
				if(logic_spinreversal/=0 .and. j==0 ) then
					allocate(transform(m,m),stat=error)
					if(error/=0) stop
					transform=0.0D0
				
					do l=1,m,1
						if(abs(symmlinkbig(szzeroindex(l),1,2))==szzeroindex(l)) then
							transform(l,l)=1.0D0
						else
							do q=l+1,m,1
								if(symmlinkbig(szzeroindex(l),1,2)==szzeroindex(q)) then
								! the fisrt column the parity is 1
								! the second column the parity is -1
									transform(l,l)=sqrt(2.0D0)/2.0D0
									transform(q,l)=sqrt(2.0D0)/2.0D0
									transform(l,q)=sqrt(2.0D0)/2.0D0
									transform(q,q)=-sqrt(2.0D0)/2.0D0
								else if(symmlinkbig(szzeroindex(l),1,2)==-szzeroindex(q)) then
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

					himp1=0
					do l=1,m,1
						if(symmlinkbig(szzeroindex(l),1,2)==szzeroindex(l)) then
							himp1=himp1+1
							call swap(coeffwork(:,himp1),coeffwork(:,l))
							call swap(coeffwork(himp1,:),coeffwork(l,:))
							call swap(transform(:,l),transform(:,himp1))
						else if(symmlinkbig(szzeroindex(l),1,2)==-szzeroindex(l)) then
							cycle
						else
							done=.true.
							do l1=1,l-1,1
								if(abs(symmlinkbig(szzeroindex(l),1,2))==szzeroindex(l1)) then
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
						write(*,*) "right diagnolization failed! himp1"
						stop
					end if
					symmlinkbigbuffer(n+1:n+himp1)=1
					call syevd(coeffwork(himp1+1:m,himp1+1:m),valuework(n+himp1+1:n+m),'V','U',info)
					if(info/=0) then
						write(*,*) "right diagnolization failed! himm1"
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
						write(*,*) "right diagnolization failed!"
						stop
					end if
				end if
				
				p=0
				do k=1,4*Rrealdim,1
					if(quantabigR(k,1)==i .and. quantabigR(k,2)==j) then
						p=p+1
						call copy(coeffwork(p,1:m),coeffresult(k,n+1:m+n))
					end if
				end do

				if(p/=m) then
					write(*,*) "p/=m failed!",p,m
					stop
				end if
				quantabigRbuffer(n+1:n+m,1)=i
				quantabigRbuffer(n+1:n+m,2)=j
				
				subspacenum(1)=subspacenum(1)+1
				subspacenum(subspacenum(1)+1)=m
			end if

		n=n+m

		end do
		
		if(j>0) then
			szl0=n
		else if(j==0 .and. logic_spinreversal/=0) then
			szzero=n-szl0
		end if

	end do
	
	if(logic_spinreversal/=0) then
		write(*,*) "szl0=",szl0,"szzero=",szzero
		if((szl0*2+szzero)/=4*Rrealdim) then
			write(*,*) "(szl0*2+szzero)/=4Rrealdim failed!","szl0=",szl0,"szzero",szzero
			stop
		end if
		if(sum(subspacenum(2:subspacenum(1)+1))/=szl0+szzero) then
			write(*,*) "------------------------"
			write(*,*) "subspacenum=szl0+szzero,failed!",subspacenum
			write(*,*) "------------------------"
			stop
		end if
	else 
		if(n/=4*Rrealdim) then
			write(*,*) "------------------------"
			write(*,*) "n/=4*Rrealdim,failed!",n
			write(*,*) "quantabigR",quantabigR(:,:)
			write(*,*) "------------------------"
			stop
		end if
		if(sum(subspacenum(2:subspacenum(1)+1))/=4*Rrealdim) then
			write(*,*) "------------------------"
			write(*,*) "subspacenum=4*Rrealdim,failed!",subspacenum
			write(*,*) "------------------------"
			stop
		end if
	end if
	


! copy the Sz>0 part and Sz=0(symmetry pair is not himself) to the symmetry pair
	if(logic_spinreversal/=0) then
		do i=1,4*Rrealdim,1
			if(quantabigR(i,2)>0) then
				coeffresult(abs(symmlinkbig(i,1,2)),szl0+szzero+1:4*Rrealdim)=coeffresult(i,1:szl0)*DBLE(sign(1,symmlinkbig(i,1,2)))
			end if
		end do
		quantabigRbuffer(szl0+szzero+1:4*Rrealdim,1)=quantabigRbuffer(1:szl0,1)
		quantabigRbuffer(szl0+szzero+1:4*Rrealdim,2)=-1*quantabigRbuffer(1:szl0,2)
		valuework(szl0+szzero+1:4*Rrealdim)=valuework(1:szl0)
	end if
	
	do i=1,4*subM,1
		if(valuework(i)<0.0D0) then
			write(*,*) "-----------------------------"
			write(*,*) "caution valuework<0.0D0",valuework(i)
			write(*,*) "-----------------------------"
			valuework(i)=0.0D0
		end if
	end do
	!if(logic_spinreversal==0) then
	!	call selectstates(valuework,4*Rrealdim,valueindex,singularvalue,subspacenum,nright)
	!else
	!	call selectstates(valuework,4*Rrealdim,valueindex,singularvalue,subspacenum,nright,szzero,pair1)
	!end if
	
	! R space select states should be corresponse to the L space states
	! so R space need to arrange again
	! when nstate==1 we can find the corresponding state in the L and R state
	! with the same eigenvalue , but when nstate/=1 there is no such condition
	if(nstate==1) then
	valueindex=0

	do i=1,subM,1

		scale1=1.0D0
		scalenum=0
		do while(.true.)
			if(singularvalue(i)>scale1) exit
			scale1=scale1*0.1D0
			scalenum=scalenum+1
			if(scalenum>=15) then
				write(*,*) "---------------------"
				write(*,*) "caution! scalenum>100"
				write(*,*) "---------------------"
				exit
			end if
		end do
		
		if(logic_spinreversal==0) then
			if(scalenum<15) then
				do j=1,subspacenum(1),1
				if((quantasmaL(i,1)+quantabigRbuffer(sum(subspacenum(2:j+1)),1)==nelecs) .and. &
				(quantasmaL(i,2)+quantabigRbuffer(sum(subspacenum(2:j+1)),2)==totalSz)) then
					diffzero=1.0D-10*(0.1D0**scalenum)
					do while(.true.)
						do k=sum(subspacenum(2:j+1)),sum(subspacenum(2:j+1))-subspacenum(j+1)+1,-1
							if(abs(valuework(k)-singularvalue(i))<diffzero) then
								Ifexist=.false.
								do p=1,i-1,1
									if(valueindex(p)==k) then
										Ifexist=.true.
										exit
									end if
								end do
								if(Ifexist==.false.) then
									valueindex(i)=k
									exit
								end if
							end if
						end do
						if(valueindex(i)/=0) exit
						diffzero=diffzero*10.0D0
					end do
					exit
				end if
				end do
			else
				do j=1,subM,1
					ifexist=.false.
					do k=1,i-1,1
						if(valueindex(k)==j) then
							ifexist=.true.
							exit
						end if
					end do
					if(ifexist==.false.) then
						valueindex(i)=j
						exit
					end if
				end do
			end if
		else
		  if(scalenum<15) then
			if(quantasmaL(i,2)<=0) then
			do j=1,subspacenum(1),1
			if((quantasmaL(i,1)+quantabigRbuffer(sum(subspacenum(2:j+1)),1)==nelecs) .and. &
			(quantasmaL(i,2)+quantabigRbuffer(sum(subspacenum(2:j+1)),2)==totalSz)) then
				diffzero=1.0D-10*(0.1D0**scalenum)
				do while(.true.)
					do k=sum(subspacenum(2:j+1)),sum(subspacenum(2:j+1))-subspacenum(j+1)+1,-1
						if(abs(valuework(k)-singularvalue(i))<diffzero) then
							Ifexist=.false.
							do p=1,i-1,1
								if(valueindex(p)==k) then
									Ifexist=.true.
									exit
								end if
							end do
							if(Ifexist==.false.) then
							if(quantasmaL(i,2)<0) then
								valueindex(i)=k
								valueindex(i-1)=k+szl0+szzero
								exit
							else if(quantasmaL(i,2)==0) then
								valueindex(i)=k
								exit
							end if
							end if
						end if
					end do
					if(valueindex(i)/=0) exit
					diffzero=diffzero*10.0D0
				end do
				exit
			end if
			end do
			end if

			if(quantasmaL(i,2)<=0 .and. valueindex(i)==0) then
			write(*,*) "----------------------------"
			write(*,*) "R space valueindex==0",i
			write(*,*) "----------------------------"
			stop
			end if
		  else
			if(i<subM .and. valueindex(i)==0) then
			do j=szl0+1,2*szl0+szzero,1
				ifexist=.false.
				do k=1,i-1,1
					if(valueindex(k)==j) then
						ifexist=.true.
						exit
					end if
				end do
				if(ifexist==.false.) then
					valueindex(i)=j
					if(j>szl0+szzero) then
					valueindex(i+1)=j-(szzero+szl0)
					end if
					exit
				end if
			end do
			else if(i==subM .and. valueindex(i)==0) then
				iffind=.false.
				do  j=szl0+1,szl0+szzero,1
					ifexist=.false.
					do k=1,i-1,1
						if(valueindex(k)==j) then
						ifexist=.true.
						exit
						end if
					end do
					if(ifexist==.false.) then
						valueindex(i)=j
						iffind=.true.
						exit
					end if
				end do
				if(iffind==.false.) then
					write(*,*) "-----------------------------------"
					write(*,*) "did not find the last index valueindex=0"
					write(*,*) "-----------------------------------"
					stop
				end if
			end if
		  end if
		end if
	end do
	else if(logic_spinreversal==0) then
		call selectstates(valuework,4*Rrealdim,valueindex,singularvalue,subspacenum,nright)
	else
		call selectstates(valuework,4*Rrealdim,valueindex,singularvalue,subspacenum,nright,szzero,szl0)
	end if
		
	!	write(*,*) valueindex(1:subM)
	!	write(*,*) valuework(valueindex(1:subM))
! check if the valueindex is right
	do i=1,subM,1
		if(valueindex(i)==0) then
			write(*,*) "----------------------------------"
			write(*,*) "splitsvdR valueindex(i)==0",i
			write(*,*) "----------------------------------"
			stop
		end if
	end do

	do i=1,subM,1
		do j=i+1,subM,1
			if(valueindex(i)==valueindex(j)) then
				write(*,*) "----------------------------------"
				write(*,*) "splitsvdR valueindex(i)=valueindex(j)",i,j,valueindex(i)
				write(*,*) "----------------------------------"
				stop
			end if
		end do
	end do


	if(4*Rrealdim>subM) then
		m=subM
	else
		m=4*Rrealdim
	end if

	do i=1,m,1
		call copy(coeffresult(1:4*Rrealdim,valueindex(i)),rightv(i,:))
		quantasmaR(i,:)=quantabigRbuffer(valueindex(i),:)
		if(logic_spinreversal/=0) then
			if(valueindex(i)<=szl0) then
				if(nstate==1) then
					symmlinksma(i,1,2)=i-1
				else
					symmlinksma(i,1,2)=i+1
				end if
			else if(valueindex(i)>szl0 .and. valueindex(i)<=szl0+szzero) then
				if(symmlinkbigbuffer(valueindex(i))==0) then
					write(*,*) "----------------------------------------"
					write(*,*) "R symmlinkbigbuffer(valueindex(i)==0 failed!"
					write(*,*) "----------------------------------------"
					stop
				end if
				symmlinksma(i,1,2)=i*symmlinkbigbuffer(valueindex(i))
			else 
				if(nstate==1) then
					symmlinksma(i,1,2)=i+1
				else
					symmlinksma(i,1,2)=i-1
				end if
			end if
		end if
	end do

		discard=1.0D0-sum(valuework(valueindex(1:subM)))
		write(*,'(A20,I4,D12.5)') "totaldiscardR=",indexRm1,discard

deallocate(valuework)
deallocate(coeffwork)
deallocate(quantabigRbuffer)
deallocate(valueindex)
deallocate(coeffbuffer)
deallocate(coeffresult)
deallocate(subspacenum)
if(logic_spinreversal/=0) then
	deallocate(szzeroindex)
	deallocate(symmlinkbigbuffer)
end if
!	tmp=0
!	do while(tmp==0) 
!		call sleep(2)
!	end do
return
end subroutine

