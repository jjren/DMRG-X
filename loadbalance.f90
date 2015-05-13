Subroutine LoadBalance
! this subroutine is to load the balance between every process
! that is to say distribute every operator between them
! in PPP model orbid is 1,2,...nprocs,1,2.....
! process 0 contains the left H and right H

	USE variables
	use communicate
    use stateOverlap

	implicit none
	! local
	integer :: i,j,error
	integer :: operanum ! operanum is the max site operator every process have
	
	call master_print_message("enter subroutine loadbalance")
	
	allocate(orbid(norbs),stat=error)
	if(error/=0) stop
	
	!allocate the every site operator
	operanum=0
	do i=1,norbs,nprocs-1
		operanum=operanum+1
		do j=1,nprocs-1,1
			if((i-1+j)<=norbs) then
				orbid((i-1)+j)=j
			else 
				exit
			end if
		end do
	end do 
	
	if(myid==0) then
		write(*,*) "operator distribute"
		write(*,*) "orbid=",orbid
	end if

!============================================================================
! allocate the work space of every operator
	if(myid/=0) then
		if(myid<=orbid(norbs)) then
			allocate(operamatbig(4*subM,4*subM,3*operanum),stat=error)
			if(error/=0) stop
			allocate(operamatsma(subM,subM,3*operanum),stat=error)
			if(error/=0) stop
		else
			allocate(operamatbig(4*subM,4*subM,3*(operanum-1)),stat=error)
			if(error/=0) stop
			allocate(operamatsma(subM,subM,3*(operanum-1)),stat=error)
			if(error/=0) stop
		end if
	else
	! 2 means the R space ;1 means the L space
		allocate(Hbig(4*subM,4*subM,2),stat=error)
		if(error/=0) stop
		allocate(Hsma(subM,subM,2),stat=error)
		if(error/=0) stop
		allocate(coeffIF(4*subM,4*subM,nstate),stat=error)
		if(error/=0) stop
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
        allocate(stateOverlapValue(nstate),stat=error)
        if(error/=0) stop
    end if
    
return
end Subroutine LoadBalance

