Subroutine LoadBalance
! this subroutine is to load the balance between every process
! that is to say distribute every operator between them
! in PPP model orbid is 1,2,...nprocs,1,2.....
! process 0 contains the left H and right H

! bond order operator distribute
! aiaj(L space) is allocated on orbid(i,1)

	USE variables
	use communicate
	use stateOverlap
	use module_sparse

	implicit none
	! local
	integer :: i,j,error
	integer ::  &
	        operanum1(nprocs-1)    , &  !operanum1 is the max site operator every process have
	        operanum2(nprocs-1)  
	
	call master_print_message("enter subroutine loadbalance")
	
	! the 1 index is the process id; the second index is the operator index on the specific process
	allocate(orbid1(norbs,2),stat=error)        
	if(error/=0) stop
	allocate(orbid2(norbs,norbs,2),stat=error)
	if(error/=0) stop
	
	! allocate the operators
	! PPP operator
	operanum1=0
	do i=1,norbs,nprocs-1
		do j=1,nprocs-1,1
			if((i-1+j)<=norbs) then
				operanum1(j)=operanum1(j)+1
				orbid1((i-1)+j,1)=j
				orbid1((i-1)+j,2)=operanum1(j)
			else 
				exit
			end if
		end do
	end do 
	
	! bond order operator
	operanum2=0   ! number of operators on specific process
	! L space
	do i=1,(norbs+1)/2,1
	do j=i,(norbs+1)/2,1
		if(bondlink(i,j)/=0) then
			orbid2(i,j,1)=orbid1(i,1)  ! be careful it is i here
			orbid2(j,i,1)=orbid1(i,1)
			operanum2(orbid1(i,1))=operanum2(orbid1(i,1))+1
			orbid2(i,j,2)=operanum2(orbid1(i,1))
			orbid2(j,i,2)=operanum2(orbid1(i,1))
		end if
	end do
	end do
	! R space
	do i=(norbs+1)/2+1,norbs,1
	do j=i,norbs,1
		if(bondlink(i,j)/=0) then
			orbid2(i,j,1)=orbid1(j,1)   ! be careful it is j here
			orbid2(j,i,1)=orbid1(j,1)
			operanum2(orbid1(j,1))=operanum2(orbid1(j,1))+1
			orbid2(i,j,2)=operanum2(orbid1(j,1))
			orbid2(j,i,2)=operanum2(orbid1(j,1))
		end if
	end do
	end do

	if(myid==0) then
		write(*,*) "PPP operator distribute"
		write(*,*) "orbid1=",orbid1(:,1)
		write(*,*) "index on every process",orbid1(:,2)
		write(*,*) "bond order operator distribute"
		write(*,*) "orbid2=",orbid2(:,:,2)
	end if

!============================================================================
! allocate the work space of every operator
	
	call AllocateArray(operanum1,operanum2)

	if(myid==0) then
		if(logic_spinreversal/=0) then
			allocate(symmlinksma(subM,1,2),stat=error)
			if(error/=0) stop
			allocate(symmlinkbig(4*subM,1,2),stat=error)
			if(error/=0) stop
		end if
	end if
!------------------------------------------------------
! allocate the quanta of every many body basis
! 1 means the total electron; 2 means the total Sz
	allocate(quantasmaL(subM,2),stat=error)
	if(error/=0) stop
	allocate(quantasmaR(subM,2),stat=error)
	if(error/=0) stop
	allocate(quantabigL(4*subM,2),stat=error)
	if(error/=0) stop
	allocate(quantabigR(4*subM,2),stat=error)
	if(error/=0) stop
!------------------------------------------------------
    if(myid==0) then
        allocate(stateOverlapValue(highestStateIndex),stat=error)
        if(error/=0) stop
    end if
    
return
end Subroutine LoadBalance

