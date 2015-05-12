Subroutine Finit_MPS
! this subroutine is do finit DMRG

	USE MPI
	USE variables
	use communicate

	implicit none

	integer :: isystem,ibegin,i,sweepbegin,sweepend
	logical :: converged
	integer :: ierr ! MPI_flag

	call master_print_message("enter in subroutine finit_MPS")

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

! only used in the restart mode
	if(mode=='r' .and. isweep/=0) then
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
	else
		sweepbegin=1
	end if

100 if(exscheme==4 .and. startedMaxOverlap == .true.) then  !He Ma
        sweepbegin = sweeps +1
        sweepend = sweeps + maxOverlapSweeps
    else
        sweepend = sweeps
    end if    
    
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
			write(*,*) isweep,"finit MPS end!"
            if(exscheme==4 .and. startedMaxOverlap == .true.) then
                write(*,*) "this is a max overlap sweep"
                write(*,*) "the energy in the middle is",sweepenergy(isweep,targetStateIndex)
            else 
			    write(*,*) "the energy in the middle is",sweepenergy(isweep,:)
            end if            
			converged = .true.
            if(reachedEnergyThresh == .false.) then    !He Ma
                converged = .false.
            end if
            if(exscheme==4 .and. startedMaxOverlap == .true.) then
                do i=1,nstate,1
                    if(i==targetStateIndex) then   ! for targetted state
                        if(abs(sweepenergy(isweep-1,i)-sweepenergy(isweep,i))>energythresh) then
                            converged = .false.
                            exit
                        end if
                    else     ! for other states, thresh is looser
                        if(abs(sweepenergy(isweep-1,i)-sweepenergy(isweep,i))>10*energythresh) then
                            converged = .false.
                            exit
                        end if
                    end if
                end do
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
		if(converged==.true. .and. startedMaxOverlap==.false.) then   !He Ma
			write(*,*) "energy converged! at sweep",isweep
            do i=1, nstate, 1
                write(*,*) "state ",i," energy at each sweep:"
                write(*,*) sweepenergy(0:isweep,i)
            end do
		else if(converged==.false. .and. startedMaxOverlap==.false.) then
			write(*,*) "maxiter reached!"
            do i=1, nstate, 1
                write(*,*) "state ",i," energy at each sweep:"
                write(*,*) sweepenergy(0:sweeps,i)
            end do
        else if(converged==.true. .and. startedMaxOverlap==.true.) then  
			write(*,*) "energy converged after another ",isweep - sweeps, " max overlap sweeps"
            write(*,*) "target state energy at each sweep:"
            write(*,*) sweepenergy(sweeps+1:isweep,targetStateIndex)
		else
			write(*,*) "max overlap maxiter reached!"
            write(*,*) "target state energy at each sweep:"
            write(*,*) sweepenergy(sweeps+1:isweep,targetStateIndex)
        end if
    end if
    
!=================================================================================
!He Ma  max overlap sweeps
    if(exscheme==4 .and. startedMaxOverlap == .false.) then
        if(myid==0) then
            write(*,*) "**************************"
		    write(*,*) "enter in max overlap sweep"
            write(*,*) "**************************"
        end if
        
        sweeps = isweep
        startedMaxOverlap = .true.
        reachedEnergyThresh = .false.
        call Renormalization(nleft+1,norbs-nright,'l')      !renormalize according to specific state
        
        goto 100
    end if
!=================================================================================
		
return
end subroutine Finit_MPS
