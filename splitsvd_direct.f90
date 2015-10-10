subroutine SplitSVD_direct(istate,leftu,rightv,singularvalue)
! this subroutine is to do svd directly from the coeffC
! only can be used in the single state case
	use kinds_mod
	use variables
	use module_sparse
	use communicate
	use blas95
	use lapack95
	use f95_precision
	
	implicit none
	
	integer :: istate
	real(kind=r8) :: leftu(4*Lrealdim,subM),rightv(subM,4*Rrealdim),&
	singularvalue(subM)

	! local
	integer :: job(8)
	real(kind=r8),allocatable :: &
		LRcoeff(:,:)        , &
		leftufull(:,:)      , &
		rightvfull(:,:)     , &
		svaluefull(:)       , &
		ww(:)               , &
		transformL(:,:)     , &
		transformR(:,:)     , &
		midmat(:,:)         , &
		midmat2(:,:)        , &
		midmat3(:,:)        , &
		midmat4(:,:)
	integer,allocatable ::  &
		quantabigLbuf(:,:) , &
		quantabigRbuf(:,:) , &
		symmlinkLbuf(:)    , &
		symmlinkRbuf(:)    , &
		lindex(:)          , &
		rindex(:)          , &
		valueindex(:)      , &
		subspacenum(:)
	integer :: error,info,tmp
	integer :: lp,ls,rp,rs,i,j,q,smadim,rhimp1,lhimp1
	integer :: mv,nv,ml,nl,mr,nr,phase,szl0,szzero
	real(kind=r8) :: discard

	call master_print_message("enter in subroutine SplitSVD_direct")

	allocate(LRcoeff(4*Lrealdim,4*Rrealdim),stat=error)
	if(error/=0) stop
	allocate(leftufull(4*Lrealdim,4*Lrealdim),stat=error)
	if(error/=0) stop
	leftufull=0.0D0
	allocate(rightvfull(4*Rrealdim,4*Rrealdim),stat=error)
	if(error/=0) stop
	rightvfull=0.0D0
	allocate(svaluefull(4*Lrealdim),stat=error)
	if(error/=0) stop
	svaluefull=0.0D0
	allocate(ww(4*Lrealdim-1),stat=error)
	if(error/=0) stop
	allocate(quantabigLbuf(4*Lrealdim,2),stat=error)
	if(error/=0) stop
	allocate(quantabigRbuf(4*Rrealdim,2),stat=error)
	if(error/=0) stop
	allocate(symmlinkLbuf(4*Lrealdim),stat=error)
	if(error/=0) stop
	allocate(symmlinkRbuf(4*Rrealdim),stat=error)
	if(error/=0) stop
	allocate(lindex(4*Lrealdim),stat=error)
	if(error/=0) stop
	allocate(rindex(4*Rrealdim),stat=error)
	if(error/=0) stop
	allocate(valueindex(subM),stat=error)
	if(error/=0) stop
	valueindex=0
	allocate(subspacenum((2*nleft+3)*(2*nleft+3)+1),stat=error)
	if(error/=0) stop
	subspacenum=0

	! recover coeffIF to its dense format
	job(1)=1
	job(2)=1
	job(3)=1
	job(4)=2
	job(5)=0
	job(6)=1
	call mkl_ddnscsr(job,4*Lrealdim,4*Rrealdim,LRcoeff(:,:),4*Lrealdim,coeffIF(:,istate),coeffIFcolindex(:,istate),coeffIFrowindex(:,istate),info)
	if(info/=0) then
		call master_print_message(info,"MoreInitFinite info/=")
		stop
	end if
	
	! record the basis index
	do i=1,4*Lrealdim,1
		lindex(i)=i
	end do
	do j=1,4*Rrealdim,1
		rindex(j)=j
	end do

	! phase represent the nonzero element in LRcoeff
	phase=logic_spinreversal*((-1)**(mod(nelecs/2,2)))

!===========================================================================================

	ml=0  ! the total number of basis have counted
	mr=0
	mv=0  ! singular value have counted

	! Sz loop
	do ls=nleft+1,-(nleft+1),-1

		! the Sz>0 rotate matrix is the same as Sz<0 rotatematrix; only need copy
		if(logic_spinreversal/=0 .and. ls<0) then
			exit
		end if

		do lp=0,2*(nleft+1),1
		
			! define R space sz and np
			rs=totalSz-ls
			rp=nelecs-lp
			
			! L space sweep
			nl=ml+1
			do i=nl,4*Lrealdim,1
				if(quantabigL(lindex(i),1)==lp .and. quantabigL(lindex(i),2)==ls) then
					ml=ml+1
					call swap(LRcoeff(ml,:),LRcoeff(i,:))
					tmp=lindex(ml)
					lindex(ml)=lindex(i)
					lindex(i)=tmp
				end if
			end do
			
			! R space sweep
			nr=mr+1  
			do j=nr,4*Rrealdim,1
				if(quantabigR(rindex(j),1)==rp .and. quantabigR(rindex(j),2)==rs) then
					mr=mr+1
					call swap(LRcoeff(:,mr),LRcoeff(:,j))
					tmp=rindex(mr)
					rindex(mr)=rindex(j)
					rindex(j)=tmp
				end if
			end do
			
			nv=mv+1
			! this subgroup have basis
			if(mr/=nr-1 .and. ml/=nl-1) then
				
				! ls=0 is very special in spin reversal symmetry
				if(logic_spinreversal/=0 .and. ls==0) then
					
					allocate(transformL(nl:ml,nl:ml),stat=error)
					if(error/=0) stop
					allocate(transformR(nr:mr,nr:mr),stat=error)
					if(error/=0) stop
					transformL=0.0D0
					transformR=0.0D0
					
					! construct the transform matrix to let the
					! basis's symmlink be themselves
					call ConstructTrans(nl,ml,1,lindex(nl:ml),transformL(nl:ml,nl:ml))
					call ConstructTrans(nr,mr,2,rindex(nr:mr),transformR(nr:mr,nr:mr))
					

					allocate(midmat(ml-nl+1,mr-nr+1),stat=error)
					if(error/=0) stop

					! transfer to the new basis
					call gemm(LRcoeff(nl:ml,nr:mr),transformR(nr:mr,nr:mr),midmat,'N','N',1.0D0,0.0D0)
					call gemm(transformL(nl:ml,nl:ml),midmat,LRcoeff(nl:ml,nr:mr),'T','N',1.0D0,0.0D0)
					
					call Sz0Swap(nl,ml,1,1,lindex(nl:ml),LRcoeff,transformL(nl:ml,nl:ml),lhimp1)
					call Sz0Swap(nr,mr,2,phase,rindex(nr:mr),LRcoeff,transformR(nr:mr,nr:mr),rhimp1)
					
					! gesvd sigular value is the descending order
					if(lhimp1/=nl-1 .and. rhimp1/=nr-1) then
						smadim=min(lhimp1-nl+1,rhimp1-nr+1)
						call gesvd(LRcoeff(nl:lhimp1,nr:rhimp1),svaluefull(mv+1:mv+smadim),leftufull(nl:lhimp1,mv+1:mv+smadim),rightvfull(mv+1:mv+smadim,nr:rhimp1),ww,'N',info)
						if(info/=0) then
							write(*,*) "Sz01 SVD failed!",info,lp,ls
						end if

						! symmlinkbigbuf==1 means the symmlink is himself
						! symmlinkbigbuf==-1 means the symmlink is -himself
						symmlinkLbuf(mv+1:mv+smadim)=1
						symmlinkRbuf(mv+1:mv+smadim)=phase
						mv=mv+smadim
					else
						write(*,*) "========================================================"
						write(*,*) "lhimp1/=nl-1 .or. rhimp1/=nr-1",lhimp1,nl-1,rhimp1,nr-1
						write(*,*) "========================================================"
					end if

					if(lhimp1/=ml .and. rhimp1/=mr) then
						smadim=min(ml-lhimp1,mr-rhimp1)
						call gesvd(LRcoeff(lhimp1+1:ml,rhimp1+1:mr),svaluefull(mv+1:mv+smadim),leftufull(lhimp1+1:ml,mv+1:mv+smadim),rightvfull(mv+1:mv+smadim,rhimp1+1:mr),ww,'N',info)
						if(info/=0) then
							write(*,*) "Sz02 SVD failed!",info,lp,ls
						end if
						! set other value to be -1.0D0
						symmlinkLbuf(mv+1:mv+smadim)=-1
						symmlinkRbuf(mv+1:mv+smadim)=phase*-1
						mv=mv+smadim
					else
						write(*,*) "=============================================="
						write(*,*) "lhimp1/=ml .or. rhimp1/=mr",lhimp1,ml,rhimp1,mr
						write(*,*) "=============================================="
					end if
					
					if(nv-1/=mv) then
						! transfer to its' inital basis representation
						allocate(midmat2(ml-nl+1,nv:mv))
						allocate(midmat3(nv:mv,mr-nr+1))
						call gemm(transformL(nl:ml,nl:ml),leftufull(nl:ml,nv:mv),midmat2,'N','N',1.0D0,0.0D0)
						leftufull(nl:ml,nv:mv)=midmat2
						call gemm(rightvfull(nv:mv,nr:mr),transformR(nr:mr,nr:mr),midmat3,'N','T',1.0D0,0.0D0)
						rightvfull(nv:mv,nr:mr)=midmat3
						deallocate(midmat2,midmat3)
					else
						write(*,*) "==========================="
						write(*,*) "nv-1==mv",lp,ls,rp,rs
						write(*,*) "==========================="
					end if
					
					deallocate(transformL,transformR,midmat)
				else
					! gesvd sigular value is the descending order
					smadim=min(mr-nr+1,ml-nl+1)
					call gesvd(LRcoeff(nl:ml,nr:mr),svaluefull(mv+1:mv+smadim),leftufull(nl:ml,mv+1:mv+smadim),rightvfull(mv+1:mv+smadim,nr:mr),ww,'N',info)
					if(info/=0) then
						write(*,*) "SVD failed!",info,lp,ls
					end if
					mv=mv+smadim
				end if

				if(nv-1/=mv) then
					! set quanta number
					quantabigLbuf(nv:mv,1)=lp
					quantabigLbuf(nv:mv,2)=ls
					quantabigRbuf(nv:mv,1)=rp
					quantabigRbuf(nv:mv,2)=rs
				
					! record the number of subgourp basis
					subspacenum(1)=subspacenum(1)+1
					subspacenum(subspacenum(1)+1)=mv-nv+1
				
					! resort to the initial basis location
					! L space
					allocate(midmat4(nl:ml,nv:mv))
					midmat4=leftufull(nl:ml,nv:mv)
					leftufull(nl:ml,nv:mv)=0.0D0
					do i=nl,ml,1
						call copy(midmat4(i,nv:mv),leftufull(lindex(i),nv:mv))
					end do
					deallocate(midmat4)

					! R space
					allocate(midmat4(nv:mv,nr:mr))
					midmat4=rightvfull(nv:mv,nr:mr)
					rightvfull(nv:mv,nr:mr)=0.0D0
					do j=nr,mr,1
						call swap(midmat4(nv:mv,j),rightvfull(nv:mv,rindex(j)))
					end do
					deallocate(midmat4)
				end if
			end if
		end do

		if(ls>0) then
			szl0=mv  ! store the number of Sz>0 Renormalized basis
		else if(ls==0) then
			szzero=mv-szl0  ! store the number of Sz=0 Renormalized basis  
		end if
	end do

	write(*,*) "szl0=",szl0,"szzero=",szzero

! copy the Sz>0 part and Sz=0(symmetry pair is not himself) to the symmetry pair

	if(logic_spinreversal/=0) then
		do i=1,4*Lrealdim,1
			if(quantabigL(i,2)>0) then
				! transfer every -1 link to 1 link
				leftufull(abs(symmlinkbig(i,1,1)),szl0+szzero+1:szl0*2+szzero)=leftufull(i,1:szl0)*DBLE(sign(1,symmlinkbig(i,1,1)))
			end if
		end do
		quantabigLbuf(szl0+szzero+1:2*szl0+szzero,1)=quantabigLbuf(1:szl0,1)
		quantabigLbuf(szl0+szzero+1:2*szl0+szzero,2)=-1*quantabigLbuf(1:szl0,2)
		
		do i=1,4*Rrealdim,1
			if(quantabigR(i,2)<0) then
				! transfer every -1 link to 1 link according to phase
				! maybe some problem here
				rightvfull(szl0+szzero+1:szl0*2+szzero,abs(symmlinkbig(i,1,2)))=rightvfull(1:szl0,i)*DBLE(sign(1,symmlinkbig(i,1,2))*phase)
			end if
		end do
		quantabigRbuf(szl0+szzero+1:2*szl0+szzero,1)=quantabigRbuf(1:szl0,1)
		quantabigRbuf(szl0+szzero+1:2*szl0+szzero,2)=-1*quantabigRbuf(1:szl0,2)
		
		svaluefull(szl0+szzero+1:2*szl0+szzero)=svaluefull(1:szl0)
	end if
	

	if(logic_spinreversal/=0) then
		call selectstates(svaluefull,szl0*2+szzero,valueindex,singularvalue,subspacenum,nleft,szzero,szl0)
	else
		call selectstates(svaluefull,mv,valueindex,singularvalue,subspacenum,nleft,szzero,szl0)
	end if

	singularvalue=singularvalue*singularvalue


!------------------------------------------------------------------------
	! check if the valueindex is right
	do i=1,subM,1
		if(valueindex(i)==0) then
			write(*,*) "----------------------------------"
			write(*,*) "splitsvd valueindex(i)==0",i
			write(*,*) "----------------------------------"
			stop
		end if
	end do
	do i=1,subM,1
		do j=i+1,subM,1
			if(valueindex(i)==valueindex(j)) then
				write(*,*) "----------------------------------"
				write(*,*) "splitsvd valueindex(i)=valueindex(j)",i,j,valueindex(i)
				write(*,*) "----------------------------------"
				stop
			end if
		end do
	end do
!------------------------------------------------------------------------

	! copy the leftufull and rightvfull to the U/V according to the valueindex
	do i=1,subM,1
		call copy(leftufull(:,valueindex(i)),leftu(:,i))
		quantasmaL(i,:)=quantabigLbuf(valueindex(i),:)
	end do
	do i=1,subM,1
		call copy(rightvfull(valueindex(i),:),rightv(i,:))
		quantasmaR(i,:)=quantabigRbuf(valueindex(i),:)
	end do

	if(logic_spinreversal/=0) then
		do i=1,subM,1
			if(valueindex(i)<=szl0) then
				symmlinksma(i,1,1)=i+1
				symmlinksma(i,1,2)=(i+1)*phase   ! Sz>0
			else if(valueindex(i)>szl0 .and. valueindex(i)<=szl0+szzero) then
				if(symmlinkLbuf(valueindex(i))==0 .or. symmlinkRbuf(valueindex(i))==0) then
					write(*,*) "----------------------------------------"
					write(*,*) "symmlinkbigbuf(valueindex(i))==0 failed!"
					write(*,*) "----------------------------------------"
					stop
				end if
				symmlinksma(i,1,1)=i*symmlinkLbuf(valueindex(i))  ! Sz=0
				symmlinksma(i,1,2)=i*symmlinkRbuf(valueindex(i))  ! Sz=0
			else 
				symmlinksma(i,1,1)=i-1
				symmlinksma(i,1,2)=(i-1)*phase  ! Sz<0
			end if
		end do
	end if

	discard=1.0D0-sum(singularvalue)
	write(*,'(A20,D12.5)') "totaldiscard=",discard

	deallocate(LRcoeff)
	deallocate(leftufull,rightvfull,svaluefull)
	deallocate(quantabigLbuf,quantabigRbuf)
	deallocate(symmlinkLbuf,symmlinkRbuf)
	deallocate(ww)
	deallocate(lindex,rindex)
	deallocate(valueindex)
	deallocate(subspacenum)
return

end subroutine splitsvd_direct


subroutine Sz0swap(n,m,Hindex,phase,index1,LRcoeff,transform,himp1)
	
	use variables,only: Lrealdim,Rrealdim,symmlinkbig
	use kinds_mod
	use blas95
	use f95_precision

	implicit none
	
	integer :: n,m,Hindex,phase
	integer :: index1(n:m)
	real(kind=r8) :: LRcoeff(4*Lrealdim,4*Rrealdim),transform(n:m,n:m)
	integer :: himp1  ! output

	! local
	integer :: l,l1
	logical :: done,thesecond

	! swap the symmlink==1 to the first few columns 
	! and the symmlink==-1 to the last few columns
	himp1=n-1  ! himself plus +1
	do l=n,m,1
		if(symmlinkbig(index1(l),1,Hindex)==phase*index1(l)) then
			himp1=himp1+1
			if(Hindex==1) then
				call swap(LRcoeff(l,:),LRcoeff(himp1,:))
			else
				call swap(LRcoeff(:,l),LRcoeff(:,himp1))
			end if
			call swap(transform(:,l),transform(:,himp1))  ! the transform matrix need to swap too
		else if(symmlinkbig(index1(l),1,Hindex)==-1*phase*index1(l)) then
			cycle
		else
			thesecond=.false.  ! the second column the symmlink is -himself
			do l1=n,l-1,1
				if(abs(symmlinkbig(index1(l),1,Hindex))==index1(l1)) then
					thesecond=.true.
					exit
				end if
			end do
			
			done=.false.
			if(phase==-1 .and. thesecond==.true.) then
				done=.true.
			end if
			if(phase==1 .and. thesecond==.false.) then
				done=.true.
			end if

			if(done==.true.) then
				himp1=himp1+1
				if(Hindex==1) then
					call swap(LRcoeff(l,:),LRcoeff(himp1,:))
				else
					call swap(LRcoeff(:,l),LRcoeff(:,himp1))
				end if
				call swap(transform(:,l),transform(:,himp1))  ! the transform matrix need to swap too
			end if
		end if
	end do
return

end subroutine Sz0Swap

subroutine ConstructTrans(n,m,Hindex,index1,transform)
	
	use variables,only : symmlinkbig
	use kinds_mod
	implicit none

	integer :: Hindex,m,n
	integer :: index1(n:m)
	real(kind=r8) :: transform(n:m,n:m)

	! local
	integer :: i,p,q

	do i=n,m,1
		if(abs(symmlinkbig(index1(i),1,Hindex))==index1(i)) then
			transform(i,i)=1.0D0  ! the symmlink is himself
		else
			do q=i+1,m,1
				if(symmlinkbig(index1(i),1,Hindex)==index1(q)) then
				! the fisrt column the parity is 1
				! the second column the parity is -1
					transform(i,i)=sqrt(2.0D0)/2.0D0
					transform(q,i)=sqrt(2.0D0)/2.0D0
					transform(i,q)=sqrt(2.0D0)/2.0D0
					transform(q,q)=-sqrt(2.0D0)/2.0D0
				else if(symmlinkbig(index1(i),1,Hindex)==-index1(q)) then
					transform(i,i)=sqrt(2.0D0)/2.0D0
					transform(q,i)=-sqrt(2.0D0)/2.0D0
					transform(i,q)=sqrt(2.0D0)/2.0D0
					transform(q,q)=sqrt(2.0D0)/2.0D0
				end if
			end do
		end if
	end do
return

end subroutine ConstructTrans



