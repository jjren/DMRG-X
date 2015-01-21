Subroutine finit_MPS
	USE MPI
	USE variables

	implicit none

	integer :: isystem,ibegin,i
	logical :: converged


	if(myid==0) then
		write(*,*) "enter in subroutine finit_MPS"
	end if

! the exactsite refer to the space that L space or R space(without sigmaL and sigmaR)
	exactsite=1
	do while(.true.)
		if(4**exactsite<=subM) then
			exactsite=exactsite+1
		else
			exactsite=exactsite-1
			exit
		end if
	end do
! ibegin is the initial L space index(without out sigmaL)
	if(mod(norbs,2)==0) then
		ibegin=norbs/2
	else
		ibegin=norbs/2+1
	end if
	
	nelecs=realnelecs
	if(myid==0) then
		if(ibegin/=nleft+1) then
			write(*,*) "-----------------------------------------------------------"
			write(*,*) "finit MPS initial L space index nleft+1/=ibegin failed!"
			write(*,*) "-----------------------------------------------------------"
		end if
	end if

	do isweep=1,sweeps,1
		do isystem=ibegin,norbs-exactsite-2,1
			nleft=isystem
			nright=norbs-isystem-2
			if(4**nleft<subM) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(4**nright<subM) then
				Rrealdim=4**nright
			else
				Rrealdim=subM
			end if
			call fromleftsweep
		end do
		
		do isystem=exactsite,norbs-exactsite-2,1
			nleft=norbs-isystem-2
			nright=isystem
			if(4**nleft<subM) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(4**nright<subM) then
				Rrealdim=4**nright
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
				Rrealdim=4**nright
			else
				Rrealdim=subM
			end if
			call fromleftsweep
		end do

		if(myid==0) then
			write(*,*) isweep,"finit MPS end!"
			write(*,*) "the lowest energy is",sweepenergy(isweep,:)
			converged=.true.
			do i=1,nstate,1
				if(abs(sweepenergy(isweep-1,i)-sweepenergy(isweep,i))>energythresh) then
					converged=.false.
					exit
				end if
			end do
			if(converged==.true.) then
				write(*,*) "energy converged! at sweep",isweep
				write(*,*) sweepenergy(0:isweep,:)
			end if
		end if
		call MPI_bcast(converged,1,MPI_logical,0,MPI_COMM_WORLD,ierr)
		
		if(converged==.true.) then
			exit
		end if
	end do
		
	
return

end subroutine
