Subroutine Finit_MPS
! this subroutine is do finit DMRG

    USE MPI
    USE variables
    use communicate
    use Renormalization_mod
    use hamiltonian_mod

    implicit none

    integer :: isystem,ibegin,i,sweepbegin
    logical :: converged
    integer :: ierr ! MPI_flag

    call master_print_message("enter in subroutine finit_MPS")

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
        if(nleft==nright .and. logic_C2/=0) then
            call C2_Copy('l')
        else
            call Enviro_Big('R')
        end if
        call Hamiltonian('i')
        call Renormalization('i')
        mode='s'
    else
        sweepbegin=1
    end if
    
    converged=.false.
    formernelecs=nelecs
    do isweep=sweepbegin,sweeps,1
        
        ! open the perturbation coefficient scheme
        if(isweep==2) then
            if(logic_perturbation/=0) then
                IfOpenperturbation=.true.
                ! in the ifopenperturbation=.false. case
                ! only the big matrix is copyed from subM space to subMp space
                call pre_perturbation('L',"sma")
                call CopyCoeff2Coeffp
            end if
        end if
        
        do isystem=ibegin,norbs-exactsite-2,1
            ! add 2,0,2,0 if nelecs does not reach realnelecs
        !   if(nelecs==formernelecs .and. nelecs<realnelecs) then
        !       nelecs=nelecs+2
        !   else
        !       formernelecs=nelecs
        !   end if

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
            if(nelecs/=realnelecs) then
                call master_print_message("nelecs/=realnelecs stop!")
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                stop
            end if

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
            write(*,*) "the energy in the middle is",sweepenergy(isweep,:)
            converged=.true.
            do i=1,nstate,1
                if(abs(sweepenergy(isweep-1,i)-sweepenergy(isweep,i))>energythresh) then
                    converged=.false.
                    exit
                end if
            end do
        end if

        call MPI_bcast(converged,1,MPI_logical,0,MPI_COMM_WORLD,ierr)
        
        if(converged==.true.) then
            exit
        end if
    end do

    if(myid==0) then
        if(converged==.true.) then
            write(*,*) "energy converged! at sweep",isweep
            write(*,*) sweepenergy(0:isweep,:)
        else
            write(*,*) "maxiter reached!"
            write(*,*) sweepenergy(0:sweeps,:)
        end if
    end if

return
end subroutine Finit_MPS
