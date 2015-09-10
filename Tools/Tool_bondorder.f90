program Toolbondorder
! this program do bondorder differece calculation after the DMRG run
	implicit none
	
	integer,parameter :: nbonds=10,norbs=10,nstates=6
	integer :: bondlink(norbs,norbs),dummybondlink(2),index1
	real(kind=8) :: bondord0(norbs,norbs,2),diffup,diffdown,dummybondord(2),bondordex(norbs,norbs,nstates)
	integer :: i,j,istate

	open(unit=10,file="bondord.out",status="old")
	open(unit=11,file="integral.inp",status="old")
	open(unit=12,file="diff_bondord.out",status="replace")
	
	bondlink=0
	do i=1,nbonds,1
		read(11,*) dummybondlink(1:2)
		bondlink(dummybondlink(1),dummybondlink(2))=1
	end do
	do i=1,norbs,1
		bondlink(i,i)=2
	end do

	! difference ex-gs
!	do istate=1,nstates*2,1
!	do i=1,norbs,1
!	do j=i,norbs,1
!		if(istate==1) then
!			read(10,*) dummybondlink(1:2),bondord0(i,j,1:2)
!		else
!			read(10,*) dummybondlink(1:2),dummybondord(1:2)
!			if(bondlink(dummybondlink(1),dummybondlink(2))/=0) then
!				diffup=dummybondord(1)-bondord0(i,j,1)
!				diffdown=dummybondord(2)-bondord0(i,j,2)
!				write(12,'(2I5,3E15.6)') i,j,diffup,diffdown,diffup+diffdown
!			end if
!		end if
!	end do
!	end do
!	end do
	
	do istate=1,nstates,1
	do i=1,norbs,1
	do j=i,norbs,1
		read(10,*) dummybondlink(1:2),dummybondord(1:2)
		bondordex(i,j,istate)=dummybondord(1)+dummybondord(2)
	end do
	end do
	end do
	
	index1=0
	! bond order
	do i=1,norbs,1
	do j=i,norbs,1
		if(bondlink(i,j)==1) then
			index1=index1+1
			write(12,'(3I5,<nstates>E15.6)') i,j,index1,bondordex(i,j,1:nstates)
		end if
	end do
	end do

	! charge density
	do i=1,norbs,1
		index1=index1+1
		write(12,'(3I5,<nstates>E15.6)') i,i,index1,bondordex(i,i,1:nstates)
	end do

	close(10)
	close(11)
	close(12)

end program



