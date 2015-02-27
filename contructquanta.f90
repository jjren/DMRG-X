Subroutine system_constructquantaL
	use variables
	use mpi

	implicit none
	

!    L+sigmaL space
	quantabigL(1:Lrealdim,1)=quantasmaL(1:Lrealdim,1)
	quantabigL(Lrealdim+1:2*Lrealdim,1)=quantasmaL(1:Lrealdim,1)+1
	quantabigL(2*Lrealdim+1:3*Lrealdim,1)=quantasmaL(1:Lrealdim,1)+1
	quantabigL(3*Lrealdim+1:4*Lrealdim,1)=quantasmaL(1:Lrealdim,1)+2

	quantabigL(1:Lrealdim,2)=quantasmaL(1:Lrealdim,2)
	quantabigL(Lrealdim+1:2*Lrealdim,2)=quantasmaL(1:Lrealdim,2)+1
	quantabigL(2*Lrealdim+1:3*Lrealdim,2)=quantasmaL(1:Lrealdim,2)-1
	quantabigL(3*Lrealdim+1:4*Lrealdim,2)=quantasmaL(1:Lrealdim,2)

return
end subroutine system_constructquantaL

Subroutine system_constructquantaR
	use variables
	use mpi

	implicit none
!    R+sigmaR space
	quantabigR(1:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)
	quantabigR(2:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)+1
	quantabigR(3:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)+1
	quantabigR(4:4*Rrealdim:4,1)=quantasmaR(1:Rrealdim,1)+2
	quantabigR(1:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)
	quantabigR(2:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)+1
	quantabigR(3:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)-1
	quantabigR(4:4*Rrealdim:4,2)=quantasmaR(1:Rrealdim,2)


return

end subroutine system_constructquantaR


Subroutine system_constructquantaRreverse
	use variables
	use mpi

	implicit none
	

!    R+sigmaR space
	quantabigR(1:Rrealdim,1)=quantasmaR(1:Rrealdim,1)
	quantabigR(Rrealdim+1:2*Rrealdim,1)=quantasmaR(1:Rrealdim,1)+1
	quantabigR(2*Rrealdim+1:3*Rrealdim,1)=quantasmaR(1:Rrealdim,1)+1
	quantabigR(3*Rrealdim+1:4*Rrealdim,1)=quantasmaR(1:Rrealdim,1)+2

	quantabigR(1:Rrealdim,2)=quantasmaR(1:Rrealdim,2)
	quantabigR(Rrealdim+1:2*Rrealdim,2)=quantasmaR(1:Rrealdim,2)+1
	quantabigR(2*Rrealdim+1:3*Rrealdim,2)=quantasmaR(1:Rrealdim,2)-1
	quantabigR(3*Rrealdim+1:4*Rrealdim,2)=quantasmaR(1:Rrealdim,2)

return
end subroutine system_constructquantaRreverse
