Subroutine System_Constructquanta(domain)

! transfer the good quantum number in M basis to 4M basis

	use variables
	use exit_mod
	use communicate

	implicit none
	character(len=1) :: domain

	call master_print_message("enter system_constructquanta subroutine!")

	if(domain=='R' .and. logic_C2==0) then
	!    R+sigmaR space
		quantabigR(1:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)
		quantabigR(2:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)+1
		quantabigR(3:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)+1
		quantabigR(4:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)+2

		quantabigR(1:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)
		quantabigR(2:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)+1
		quantabigR(3:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)-1
		quantabigR(4:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)

	else if(domain=='L') then
	!    L+sigmaL space
		quantabigL(1:Lrealdim,1)=quantasmaL(1:Lrealdim,1)
		quantabigL(Lrealdim+1:2*Lrealdim,1)=quantasmaL(1:Lrealdim,1)+1
		quantabigL(2*Lrealdim+1:3*Lrealdim,1)=quantasmaL(1:Lrealdim,1)+1
		quantabigL(3*Lrealdim+1:4*Lrealdim,1)=quantasmaL(1:Lrealdim,1)+2

		quantabigL(1:Lrealdim,2)=quantasmaL(1:Lrealdim,2)
		quantabigL(Lrealdim+1:2*Lrealdim,2)=quantasmaL(1:Lrealdim,2)+1
		quantabigL(2*Lrealdim+1:3*Lrealdim,2)=quantasmaL(1:Lrealdim,2)-1
		quantabigL(3*Lrealdim+1:4*Lrealdim,2)=quantasmaL(1:Lrealdim,2)

	else if(domain=='L') then
	!    Rreverse+sigmaR space
		quantabigR(1:Rrealdim,1)=quantasmaR(1:Rrealdim,1)
		quantabigR(Rrealdim+1:2*Rrealdim,1)=quantasmaR(1:Rrealdim,1)+1
		quantabigR(2*Rrealdim+1:3*Rrealdim,1)=quantasmaR(1:Rrealdim,1)+1
		quantabigR(3*Rrealdim+1:4*Rrealdim,1)=quantasmaR(1:Rrealdim,1)+2

		quantabigR(1:Rrealdim,2)=quantasmaR(1:Rrealdim,2)
		quantabigR(Rrealdim+1:2*Rrealdim,2)=quantasmaR(1:Rrealdim,2)+1
		quantabigR(2*Rrealdim+1:3*Rrealdim,2)=quantasmaR(1:Rrealdim,2)-1
		quantabigR(3*Rrealdim+1:4*Rrealdim,2)=quantasmaR(1:Rrealdim,2)
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/=R failed!")
	end if

return
end subroutine system_constructquanta

