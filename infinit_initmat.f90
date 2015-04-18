Subroutine Infinit_InitMat(domain)
! construct the L/R subspace operator matrix initially in infinit process

	use variables
	use communicate
	use exit_mod

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
		operaindex=1
		Hindex=1
	else if(domain=='R') then
		orbindex=norbs
		if(mod(norbs,nprocs-1)==0) then
			operaindex=norbs/(nprocs-1)
		else
			operaindex=norbs/(nprocs-1)+1
		end if
		Hindex=2
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/R failed!")
	end if

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

	if(myid==orbid(orbindex)) then
		operamatsma(:,:,3*operaindex-2:3*operaindex)=0.0D0
		! here only this matrix is set to zero
		! do not touch other matrix( the R space matrix )
		! store the create operator
		operamatsma(2,1,3*operaindex-2)=1.0D0
		operamatsma(4,3,3*operaindex-2)=1.0D0
		operamatsma(3,1,3*operaindex-1)=1.0D0
		operamatsma(4,2,3*operaindex-1)=-1.0D0
		operamatsma(1,1,3*operaindex)=-nuclQ(orbindex)
		operamatsma(2,2,3*operaindex)=1.0D0-nuclQ(orbindex)
		operamatsma(3,3,3*operaindex)=1.0D0-nuclQ(orbindex)
		operamatsma(4,4,3*operaindex)=2.0D0-nuclQ(orbindex)
	else if(myid==0) then
		Hsma(1:4,1:4,Hindex)=0.0D0
		Hsma(2,2,Hindex)=t(orbindex,orbindex)
		Hsma(3,3,Hindex)=t(orbindex,orbindex)
		Hsma(4,4,Hindex)=2*t(orbindex,orbindex)+hubbardU(orbindex)
		
		if(logic_spinreversal/=0) then
			symmlinksma(:,:,Hindex)=0
			symmlinksma(1,1,Hindex)=1
			symmlinksma(2,1,Hindex)=3
			symmlinksma(3,1,Hindex)=2
			symmlinksma(4,1,Hindex)=-4
		end if
	end if

return
end Subroutine Infinit_InitMat






