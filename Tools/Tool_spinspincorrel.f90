program Tool_SpinSpinCorrel
	implicit none
	integer :: norbs,nstate
	real(kind=8),allocatable :: SpinSpinCorrelation0(:,:,:,:,:)
	integer :: npair,iorb1,iorb2
	integer :: jstate,istate,i

	open(unit=1000,file="spinspincorrelation.out",form="unformatted",status="old")

	read(1000) norbs,nstate
	allocate(SpinSpinCorrelation0(norbs,norbs,2,nstate,nstate))

	read(1000) SpinSpinCorrelation0

	close(1000)

	open(unit=1001,file="Tool_SpinSpinCorrel.inp",status="old")
	open(unit=101,file="Tool_SpinSpinCorrel.out",status="replace")
	read(1001,*) npair
	do i=1,npair,1
		read(1001,*) iorb1,iorb2
		do istate=1,nstate,1
		do jstate=1,nstate,1
			write(101,*) jstate,istate,SpinSpinCorrelation0(iorb1,iorb2,1:2,jstate,istate)
			write(101,*) jstate,istate,SpinSpinCorrelation0(iorb2,iorb1,1:2,jstate,istate)
		end do
		end do
	end do
	close(1001)
	close(101)
	
	deallocate(SpinSpincorrelation0)
	write(*,*) "Tool_SpinSpinCorrelation happy ending..."


end program 
