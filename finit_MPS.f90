Subroutine Finit_MPS
! this subroutine is do finit DMRG

	USE MPI
	USE variables
	use communicate
    use stateOverlap

	implicit none

	integer :: isystem,ibegin,i,sweepbegin,sweepend
	logical :: converged
	integer :: ierr ! MPI_flag
    real(kind=r8) :: starttime,endtime

	call master_print_message("enter in subroutine finit_MPS")
    starttime=MPI_WTIME()

! the exactsite refer to the space that L space or R space that can be accurately discribe
! (without sigmaL and sigmaR)
	exactsite=1
	do while(.true.)
		if(4**exactsite<=subM) then
			exactsite=exactsite+1
		else
			exactsite=exactsite-1
			exit
		end if
    end do
    
    if(mod(norbs,2)==0) then
		stepPerSweep=2*(2*(norbs/2-exactsite)-1)
	else
		stepPerSweep=2*(2*(norbs/2-exactsite)-1)+2
	end if

! ibegin is the initial L space index(without sigmaL)
	if(mod(norbs,2)==0) then
		ibegin=norbs/2
	else
		ibegin=norbs/2+1
	end if
	
	if(myid==0) then
		if(ibegin/=nleft+1) then
			call master_print_message("finit MPS initial L space index nleft+1/=ibegin failed!")
			stop
		end if
	end if

!  only used in the restart mode
 	if(mode=='r' .and. isweep/=0 .and. (exscheme/=4 .or. startedStateSpecific==.false.)) then
 		nelecs=realnelecs
 		sweepbegin=isweep
 		isweep=isweep-1
 		nleft=ibegin-1
 		nright=norbs-ibegin-1
 		Lrealdim=subM
 		Rrealdim=subM
 		sweepenergy(0:isweep-1,:)=0.0D0
 		call Enviro_Big('L')
 		call Enviro_Big('R')
 		call Hamiltonian('i')
 		call Renormalization(nleft+1,norbs-nright,'i')
    end if
    
!=================================================================================
    if(exscheme==4 .and. startedStateSpecific == .true.) then
        call initStateSpecific()
        sweepbegin = sweeps + 1     ! sweeps is the final sweep number of state averaged DMRG
        sweepend = sweeps + maxStateSpecificSweeps
    else
        sweepbegin = 1
        sweepend = maxSweeps
    end if
!=================================================================================

    do isweep=sweepbegin,sweepend,1
		
		do isystem=ibegin,norbs-exactsite-2,1
			! add nelecs 2 by 2 if the nelecs does not reach realnelecs
			formernelecs=nelecs
			if(nelecs<realnelecs) then
				nelecs=nelecs+2
			end if
			if(nelecs>realnelecs) then
				nelecs=realnelecs
			end if
			
			call master_print_message(nelecs,"nelecs=")

			nleft=isystem
			nright=norbs-isystem-2
			
			if(nleft<=exactsite) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(nright<=exactsite) then
				Rrealdim=4**nright
			else
				Rrealdim=subM
			end if
			call Sweep('l')
		end do

! we assume the nelecs have reach realnelecs now

		do isystem=exactsite,norbs-exactsite-2,1
			nleft=norbs-isystem-2
			nright=isystem
			if(nleft<=exactsite) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(nright<=exactsite) then
				Rrealdim=4**nright
			else
				Rrealdim=subM
			end if
			call Sweep('r')
		end do

		do isystem=exactsite,ibegin-1,1
			nleft=isystem
			nright=norbs-isystem-2
			if(nleft<=exactsite) then
				Lrealdim=4**nleft
			else
				Lrealdim=subM
			end if
			if(nright<=exactsite) then
				Rrealdim=4**nright
			else
				Rrealdim=subM
			end if
			call Sweep('l')
		end do

		if(myid==0) then
            if(exscheme==4 .and. startedStateSpecific == .true.) then
                write(*,*) isweep,"state specific finit MPS end!"
                write(*,*) "the energy in the middle is",stateSpecificSweepEnergy(isweep)
                write(*,*) "the target state index in the middle is",storedStateIndex(isweep)
            else 
                write(*,*) isweep,"finit MPS end!"
			    write(*,*) "the energy in the middle is",sweepenergy(isweep,1:nstate)
            end if
			converged = .true.
            if(reachedEnergyThresh == .false.) then    !He Ma
                converged = .false.
            end if
            if(exscheme==4 .and. startedStateSpecific == .true.) then
                if(abs(stateSpecificSweepEnergy(isweep)-stateSpecificSweepEnergy(isweep-1))>energythresh) then
                    converged = .false.
                end if
            else
			    do i=1,nstate,1
				    if(abs(sweepenergy(isweep-1,i)-sweepenergy(isweep,i))>energythresh) then
					    converged=.false.
					    exit
				    end if
                end do
            end if
		end if

		call MPI_bcast(converged,1,MPI_logical,0,MPI_COMM_WORLD,ierr)
		
		if(converged==.true.) then
			exit
		end if
	end do

	if(myid==0) then
        if(exscheme==4 .and. startedStateSpecific) then
            select case(converged)
            case(.true.)
                write(*,*) "energy converged after another ",isweep - sweeps, " state specific sweeps"
            case(.false.)
                write(*,*) "max overlap maxiter reached!"
            end select
            write(*,*) "target state energy and its index at each sweep:"
            do i=sweeps+1,isweep,1
                write(*,*) stateSpecificSweepEnergy(i),storedStateIndex(i)," at sweep", i
            end do
            call checkStateSpecificResults()
            call cleanStateSpecificVariables()
        else
            sweeps = isweep         ! "sweeps" stores how many finit sweeps have been done
            select case(converged)
            case(.true.)
                write(*,*) "energy converged! at sweep",isweep
            case(.false.)
                write(*,*) "maxiter reached!"
            end select
            do i=1, nstate, 1
                write(*,*) "state ",i," energy at each sweep:"
                write(*,*) sweepenergy(0:isweep,i)
            end do
        end if
    end if
    
    call MPI_bcast(sweeps,1,MPI_integer4,0,MPI_COMM_WORLD,ierr)
    
    endtime=MPI_WTIME()
    select case(startedStateSpecific)
    case(.false.)
        call master_print_message(endtime-starttime,"Finite DMRG Runtime:")
    case(.true.)
        call master_print_message(endtime-starttime,"State-Specific Finite DMRG Runtime:")
    end select

end subroutine Finit_MPS
