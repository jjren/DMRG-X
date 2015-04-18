Subroutine infinit_smallL
! construct the L subspace operator matrix initially

	use variables
	use communicate

	implicit none
	
	! local
	integer :: i
	
	call master_print_message("enter in subroutine infinit_smallL")

	if(nleft==1) then
		quantasmaL=0
		quantasmaL(1,1)=0
		quantasmaL(1,2)=0
		quantasmaL(2,1)=1
		quantasmaL(2,2)=1
		quantasmaL(3,1)=1
		quantasmaL(3,2)=-1
		quantasmaL(4,1)=2
		quantasmaL(4,2)=0
		if(myid==orbid(1)) then
			operamatsma(:,:,1:3)=0.0D0
			! here only this matrix is set to zero
			! do not touch other matrix( the R space matrix )
			! store the create operator
			operamatsma(2,1,1)=1.0D0
			operamatsma(4,3,1)=1.0D0
			operamatsma(3,1,2)=1.0D0
			operamatsma(4,2,2)=-1.0D0
			operamatsma(1,1,3)=-nuclQ(nleft)
			operamatsma(2,2,3)=1.0D0-nuclQ(nleft)
			operamatsma(3,3,3)=1.0D0-nuclQ(nleft)
			operamatsma(4,4,3)=2.0D0-nuclQ(nleft)
		else if(myid==0) then
			Hsma(1:4,1:4,1)=0.0D0
			Hsma(2,2,1)=t(1,1)
			Hsma(3,3,1)=t(1,1)
			Hsma(4,4,1)=2*t(1,1)+hubbardU(1)
			
			if(logic_spinreversal/=0) then
				symmlinksma(:,:,1)=0
				symmlinksma(1,1,1)=1
				symmlinksma(2,1,1)=3
				symmlinksma(3,1,1)=2
				symmlinksma(4,1,1)=-4
			end if
			
		end if
	end if
return
end Subroutine infinit_SmallL






