subroutine excitedbasis(leftu,rightv,singularvalue,leftu2,rightv2,&
			quantasmaL2,quantasmaR2,symmlinksma2)

! this subroutine is to get the basis in the excited exsheme=2 algorithom
	use mpi
	use variables
	use blas95
	use f95_precision

	implicit none
	real(kind=8) :: leftu(4*Lrealdim,subM),rightv(subM,4*Rrealdim),singularvalue(subM*nstate),&
		leftu2(4*Lrealdim,subM*nstate),rightv2(subM*nstate,4*Rrealdim)
	integer :: quantasmaL2(subM*nstate,2),quantasmaR2(subM*nstate,2)
	integer,optional :: symmlinksma2(subM*nstate,1,2)
	
	integer,allocatable :: index1(:),former(:),index2(:)
	integer :: error,dim1,i,j,n,k,l,L1,L2,R1,R2
	real(kind=8),allocatable :: value(:)
	real(kind=8) :: norm
	real(kind=8) :: rid=1.0D-10

	write(*,*) "enter in excitedbasis subroutine!"
	dim1=2*subM

	allocate(index1(dim1),stat=error)
	if(error/=0) stop
	allocate(index2(dim1),stat=error)
	if(error/=0) stop
	allocate(value(dim1),stat=error)
	if(error/=0) stop
	allocate(former(dim1),stat=error)
	if(error/=0) stop

	value=0.0D0
	index1=0
	if(logic_spinreversal==0) then
		do i=1,nstate*subM,1
			do j=1,dim1,1
				if(singularvalue(i)>value(j)) then
					index1(j+1:dim1)=index1(j:dim1-1)
					value(j+1:dim1)=value(j:dim1-1)
					index1(j)=i
					value(j)=singularvalue(i)
					exit
				end if
			end do
		end do
		index2=index1
! L space Gram Schmit
		do i=(nleft+1),-nleft-1,-1
			do j=0,2*(nleft+1),1
				n=0
				former=0
				do k=1,dim1,1
				if(index1(k)/=0) then
				if(quantasmaL2(index1(k),1)==j .and. &
				quantasmaL2(index1(k),2)==i) then
					if(n>0) then
						do l=1,n,1
							norm=dot(leftu2(:,index1(k)),leftu2(:,former(l)))
							leftu2(:,index1(k))=leftu2(:,index1(k))-norm*leftu2(:,former(l))
						end do
					end if
					norm=dot(leftu2(:,index1(k)),leftu2(:,index1(k)))
					if(norm<rid) then
						index1(k)=0
					else
						leftu2(:,index1(k))=leftu2(:,index1(k))/sqrt(norm)
						n=n+1
						former(n)=index1(k)
					end if
				end if
				end if
				end do
			end do
		end do

		n=0
		do i=1,dim1,1
			if(index1(i)/=0) then
				n=n+1
				leftu(:,n)=leftu2(:,index1(i))
				quantasmaL(n,:)=quantasmaL2(index1(i),:)
				if(n==subM) then
					exit
				end if
			end if
		end do

		if(n/=subM) then
			write(*,*) "-----------------------"
			write(*,*) " L not get subM states",n
			write(*,*) "-----------------------"
			stop
		end if

! R space Gram Schmit
		do i=(nright+1),-nright-1,-1
			do j=0,2*(nright+1),1
				n=0
				former=0
				do k=1,dim1,1
				if(index2(k)/=0) then
				if(quantasmaR2(index2(k),1)==j .and. &
				quantasmaR2(index2(k),2)==i) then
					if(n>0) then
						do l=1,n,1
							norm=dot(rightv2(index2(k),:),rightv2(former(l),:))
							rightv2(index2(k),:)=rightv2(index2(k),:)-norm*rightv2(former(l),:)
						end do
					end if
					norm=dot(rightv2(index2(k),:),rightv2(index2(k),:))
					if(norm<rid) then
						index2(k)=0
					else
						rightv2(index2(k),:)=rightv2(index2(k),:)/sqrt(norm)
						n=n+1
						former(n)=index2(k)
					end if
				end if
				end if
				end do
			end do
		end do

		n=0
		do i=1,dim1,1
			if(index2(i)/=0) then
				n=n+1
				rightv(n,:)=rightv2(index2(i),:)
				quantasmaR(n,:)=quantasmaR2(index2(i),:)
				if(n==subM) then
					exit
				end if
			end if
		end do

		if(n/=subM) then
			write(*,*) "-----------------------"
			write(*,*) "R not get subM states",n
			write(*,*) "-----------------------"
			stop
		end if
	end if

	write(*,*) "index1",index1
	write(*,*) "index2",index2
	
	deallocate(index1)
	deallocate(index2)
	deallocate(former)
	deallocate(value)


return
end subroutine
