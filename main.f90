program main
! This is a DMRG_MPS program quantum chemistry
! only PPP model has been written

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% Jiajun Ren                      %
!% Zhigang Shuai's Group           %
!% Tsinghua University             %
!% Email: jiajunren0522@126.com    %   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    USE variables
    use communicate
    use exit_mod
    use mpi
    use MeanField
    use analysis_mod
    use Peierls_mod

    implicit none
    
    real(kind=r8) :: starttime,endtime
    integer :: ierr
    logical :: ifpeierlsconverge
    integer :: ipeierlsloop
    
    call init_communicate
    starttime=MPI_WTIME()
    
    if(nprocs<2) then
        call exit_DMRG(sigAbort,"nprocs<2 failed!")
    end if

    ! read the input files
    call ReadInput

    ! write FCIDUMP which is the interface to BLOCK
    call CreatFCIDUMP
    
    ! SCF mean field procedure
    if(logic_meanfield==1 .and. myid==0) then
        call SCF_driver
    end if
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if(logic_Peierls==0) then
        npeierlsloops=1
    end if

    do ipeierlsloop=1,npeierlsloops,1
        
        if(logic_Peierls/=0) then
            call master_print_message(ipeierlsloop,"peierls iloop=")
            call Peierls_init("DMRG")
        end if 

        ! allocate the operator to different process
        call LoadBalance
        
        ! do infinit DMRG process
        call Infinit_MPS
        
        ! do finit DMRG process
        call Finit_MPS
        
        ! do wave function analysis
        call Analysis

        if(logic_Peierls==1 .and. ipeierlsloop==1) then
            call SYSTEM("mkdir 1-back")
            call SYSTEM("cp *.out 1-back")
        end if
            

        if(logic_Peierls==1) then
            if(myid==0) then
                call Peierls_driver(ifpeierlsconverge)
            end if
            call MPI_BCAST(ifpeierlsconverge,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        end if

        call Free_DMRG
        
        if(logic_Peierls==1) then
            if(ifpeierlsconverge==.true.) then 
                call master_print_message("DMRG Peierls converge!")
                exit
            else
                call master_print_message("DMRG Peierls not converge!")
            end if
        end if
    end do
    
    call Free_Program

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    endtime=MPI_WTIME()
    call master_print_message(endtime-starttime,"RUNTIME:")
    call exit_DMRG(0,"Program DMRG-X end successfully")

end program main
