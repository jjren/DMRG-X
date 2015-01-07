Subroutine finit_MPS
	USE MPI
	USE variables

	implicit none

	integer :: isweep,isystem,ibegin


	if(myid==0) then
		write(*,*) "enter in subroutine finit_MPS"
	end if

! the exactsite refer to the space that L space or R space(without sigmaL and sigmaR)
	exactsite=1
	do while(.true.)
		if(4**exactsite<=subM) then
			exactsite=exactsite+1
		else
			exit
		end if
	end do
! ibegin is the initial L space index(without out sigmaL)
	if(mod(norbs,2)==0) then
		ibegin=norbs/2
	else
		ibegin=norbs/2+1
	end if

	if(myid==0) then
		if(ibegin/=nleft+1) then
			write(*,*) "-----------------------------------------------------------"
			write(*,*) "finit MPS initial L space index nleft+1/=ibegin failed!"
			write(*,*) "-----------------------------------------------------------"
		end if
	end if

	do isweep=1,sweeps,1
		do isystem=ibegin,norbs-exactsite-3,1
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
		
		do isystem=exactsite,norbs-exactsite-3,1
			nleft=norbs-isystem-2
			nright=isystem
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

		do isystem=exactsite,ibegin-1,1
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
