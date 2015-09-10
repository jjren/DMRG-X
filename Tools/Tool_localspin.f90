Program Tool_localspin
! this program is to calculate the local spin of some part of the molecule
	implicit none
	
	integer,parameter :: norbs=52,nstates=6
	integer :: ibegin,iend
	real(kind=8) :: localspin(norbs,norbs,nstates),partspin(nstates)
	integer :: index1,index2
	integer :: i,j,istate

	write(*,*) "please input the begin atom and end atom index"
	read(*,*) ibegin,iend

	open(unit=10,file="spin.out",status="old")
	open(unit=11,file="partspin.out",status="replace")

	partspin=0.0D0
	do istate=1,nstates,1
		do i=1,norbs,1
		do j=1,norbs,1
			read(10,*) index1,index2,localspin(index1,index2,istate)
			if((index1>=ibegin .and. index1<=iend) .and. (index2>=ibegin .and. index2<=iend)) then
				partspin(istate)=partspin(istate)+localspin(index1,index2,istate)
			end if
		end do
		end do
		read(10,*)
		read(10,*)
		read(10,*)
		write(11,*) istate,partspin(istate)
	end do

	close(10)
	close(11)

end program


	

