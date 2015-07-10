Subroutine Infinit_InitMat(domain)
! construct the L/R subspace operator matrix initially in infinit process

	use variables
	use communicate
	use exit_mod
	use module_sparse
	use onesitematrix
	use BondOrder_mod
	use LocalSpin_mod

	implicit none
	
	character(len=1) :: domain    ! L or R
	! local
	integer :: orbindex,operaindex,Hindex
	! orbindex in orbital index
	! operaindex is on this process operator index
	! Hindex is HL/1 HR/2
	
	call master_print_message("enter in subroutine infinit_initmat")
	
	
	if(domain=='L') then
		orbindex=1
		Hindex=1
	else if(domain=='R') then
		orbindex=norbs
		Hindex=2
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
	end if
	operaindex=orbid1(orbindex,2)

	! set good quantum number
	if(domain=='R') then
		quantasmaR=0
		quantasmaR(1,1)=0
		quantasmaR(1,2)=0
		quantasmaR(2,1)=1
		quantasmaR(2,2)=1
		quantasmaR(3,1)=1
		quantasmaR(3,2)=-1
		quantasmaR(4,1)=2
		quantasmaR(4,2)=0
	else if(domain=='L') then
		quantasmaL=0
		quantasmaL(1,1)=0
		quantasmaL(1,2)=0
		quantasmaL(2,1)=1
		quantasmaL(2,2)=1
		quantasmaL(3,1)=1
		quantasmaL(3,2)=-1
		quantasmaL(4,1)=2
		quantasmaL(4,2)=0
	end if

	if(myid==orbid1(orbindex,1)) then
		! set the matrix be 0
		smarowindex1(:,3*operaindex-2:3*operaindex)=1
		call ConstructOnesiteMatrix(orbindex)
		! here only this matrix is set to zero
		! do not touch other matrix( the R space matrix )
		! store the create operator
		operamatsma1(1:4,3*operaindex-2:3*operaindex)=onesitemat(:,1:3)
		smarowindex1(1:5,3*operaindex-2:3*operaindex)=osmrowindex(:,1:3)
		smacolindex1(1:4,3*operaindex-2:3*operaindex)=osmcolindex(:,1:3)
		
	else if(myid==0) then
		! set the matrix be 0
		Hsmarowindex(:,Hindex)=1
		call ConstructOnesiteMatrix(orbindex)
		Hsma(1:4,Hindex)=onesitemat(:,6)
		Hsmarowindex(1:5,Hindex)=osmrowindex(:,6)
		Hsmacolindex(1:4,Hindex)=osmcolindex(:,6)

		if(logic_spinreversal/=0) then
			symmlinksma(:,:,Hindex)=0
			symmlinksma(1,1,Hindex)=1
			symmlinksma(2,1,Hindex)=3
			symmlinksma(3,1,Hindex)=2
			symmlinksma(4,1,Hindex)=-4
		end if
	end if
	
	! bond order onsite matrix
	if(logic_bondorder==1) then
		call init_BOmat(orbindex)
	end if

	! local spin onsite matrix
	if(logic_localspin==1) then
		call init_localspinmat(orbindex)
	end if

return
end Subroutine Infinit_InitMat






