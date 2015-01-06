Subroutine Spin_reversalmatL

	USE mpi
	USE variables
	USE mathlib

	implicit none

	integer :: error
	! L(R) realdim is the L and R space real dimension

	if(myid==0 .and. logic_spinreversal/=0) then
		write(*,*) "enter in Spin_reversalmatL subroutine"
! contruct the parity operator in 4M basis
		call directproduct(adaptedsma(1:Lrealdim,1:Lrealdim,1),Lrealdim,parityonesitemat,4,adaptedbig(1:4*Lrealdim,1:4*Lrealdim,1))
	!	write(*,*) "left adaptedbig"
	!	write(*,*) adaptedbig(1:4*Lrealdim,1:4*Lrealdim,1)
	!	write(*,*) "right adaptedbig"
	!	write(*,*) adaptedbig(1:4*Rrealdim,1:4*Rrealdim,2)
	end if
		
return

end Subroutine Spin_reversalmatL


Subroutine Spin_reversalmatR

	USE mpi
	USE variables
	USE mathlib

	implicit none

	integer :: error
	! L(R) realdim is the L and R space real dimension

	if(myid==0 .and. logic_spinreversal/=0) then
		write(*,*) "enter in Spin_reversalmatR subroutine"
! contruct the parity operator in 4M basis
		call directproduct(parityonesitemat,4,adaptedsma(1:Rrealdim,1:Rrealdim,2),Rrealdim,adaptedbig(1:4*Rrealdim,1:4*Rrealdim,2))
	!	write(*,*) "left adaptedbig"
	!	write(*,*) adaptedbig(1:4*Lrealdim,1:4*Lrealdim,1)
	!	write(*,*) "right adaptedbig"
	!	write(*,*) adaptedbig(1:4*Rrealdim,1:4*Rrealdim,2)
	end if
		
return

end Subroutine Spin_reversalmatR

Subroutine parity_onesitematrix
	USE mpi
	USE variables

	implicit none

	parityonesitemat=0.0D0
	parityonesitemat(1,1)=1.0D0
	parityonesitemat(4,4)=-1.0D0
	parityonesitemat(2,3)=1.0D0
	parityonesitemat(3,2)=1.0D0

	return
end Subroutine parity_onesitematrix


Subroutine adaptedtrans
	USE mpi
	USE variables
	use mathlib

	implicit none
	
	real(kind=8),allocatable :: adaptedtransmat(:,:)
	integer :: i,error

	if(myid==0 .and. logic_spinreversal/=0) then
		write(*,*) "enter in adaptedtrans subroutine"
		allocate(adaptedtransmat(16*subM*subM,16*subM*subM),stat=error)
		if(error/=0) stop

		call directproduct(adaptedbig(1:Lrealdim,1:Lrealdim,1),Lrealdim,adaptedbig(1:Rrealdim,1:Rrealdim,2),Rrealdim,adaptedtransmat(1:Rrealdim*Lrealdim,1:Rrealdim*Lrealdim))
		do i=1,Lrealdim*Rrealdim,1
			adaptedtransmat(i,i)=1.0D0+adaptedtransmat(i,i)*logic_spinreversal
		end do
		write(*,*) "adaptedtran="
		write(*,*) adaptedtransmat
	end if
	
	return

end Subroutine adaptedtrans
