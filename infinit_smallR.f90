Subroutine infinit_smallR
! construct the R subspace operator matrix 

	use variables
	use mpi

	implicit none

	integer :: i,operaindex
	
	if(myid==0) then
		write(*,*) "enter in subroutine infinit_smallR"
	end if

	if(nright==1) then
		quantasmaR=0
		quantasmaR(1,1)=0
		quantasmaR(1,2)=0
		quantasmaR(2,1)=1
		quantasmaR(2,2)=1
		quantasmaR(3,1)=1
		quantasmaR(3,2)=-1
		quantasmaR(4,1)=2
		quantasmaR(4,2)=0
		if(myid==orbid(norbs)) then
			if(mod(norbs,nprocs-1)==0) then
				operaindex=norbs/(nprocs-1)
			else
				operaindex=norbs/(nprocs-1)+1
			end if
			operamatsma(:,:,3*operaindex-2:3*operaindex)=0.0D0
			operamatsma(2,1,3*operaindex-2)=1.0D0
			operamatsma(4,3,3*operaindex-2)=1.0D0
			operamatsma(3,1,3*operaindex-1)=1.0D0
			operamatsma(4,2,3*operaindex-1)=-1.0D0
			operamatsma(1,1,3*operaindex)=-nuclQ(norbs)
			operamatsma(2,2,3*operaindex)=1.0D0-nuclQ(norbs)
			operamatsma(3,3,3*operaindex)=1.0D0-nuclQ(norbs)
			operamatsma(4,4,3*operaindex)=2.0D0-nuclQ(norbs)
		else if(myid==0) then
			Hsma(1:4,1:4,2)=0.0D0
			Hsma(2,2,2)=t(norbs,norbs)
			Hsma(3,3,2)=t(norbs,norbs)
			Hsma(4,4,2)=2*t(norbs,norbs)+hubbardU(norbs)
			
			adaptedsma(:,:,2)=0.0D0
			adaptedsma(1,1,2)=1.0D0
			adaptedsma(4,4,2)=-1.0D0
			adaptedsma(2,3,2)=1.0D0
			adaptedsma(3,2,2)=1.0D0
		end if
	end if
return
end Subroutine infinit_SmallR






