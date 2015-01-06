Subroutine finit_MPS
	USE MPI
	USE variables

	implicit none

	integer :: isweep,isystem
	
	exactsite=1
	do while(.true.)
		if(4**exactsite<=subM) then
			exactsite=exactsite+1
		else
			exit
		end if
	end do

	do isweep=1,sweeps,1
		do isystem=norbs/2,norbs-exactsite-2,1
			nleft=isystem
			nright=norbs-isystem-2
			if(4**nleft<subM) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(4**nright<subM) then
				Rrealdim=4**nleft
			else
				Rrealdim=subM
			end if
			call fromleftsweep
		end do
		
		do isystem=norbs-exactsite+1,exactsite+3,-1
			nleft=isystem
			nright=norbs-isystem-2
			if(4**nleft<subM) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(4**nright<subM) then
				Rrealdim=4**nleft
			else
				Rrealdim=subM
			end if
			call fromrightsweep
		end do

		do isystem=exactsite,norbs-1,1
			nleft=isystem
			nright=norbs-isystem-2
			if(4**nleft<subM) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(4**nright<subM) then
				Rrealdim=4**nleft
			else
				Rrealdim=subM
			end if
			call fromleftsweep
		end do
	end do
		
	
return

end subroutine
