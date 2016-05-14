Subroutine LoadBalance
! this subroutine is to load the balance between every process
! that is to say distribute every operator between them
! in PPP model orbid is 1,2,...nprocs,1,2.....
! process 0 contains the left H and right H

! bond order operator distribute
! aiaj(L space) is allocated on orbid(i,1)

    USE variables
    use communicate
    use module_sparse
    use MKL_SERVICE
    use checkmem_mod
    use basisindex_mod

    implicit none
    ! local
    integer :: i,j,error
    
    call master_print_message("enter subroutine loadbalance")
    
    ! set every process threads
    if(nthreads(1)/=0) then ! if nthreads(1)==0 then using the enviroment variables
        if(myid==0) then
            call MKL_SET_NUM_THREADS(nthreads(1))
        else
            call MKL_SET_NUM_THREADS(nthreads(2))
        end if
    end if

    ! the 1 index is the process id; the second index is the operator index on the specific process
    allocate(operanum1(nprocs-1))
    allocate(orbid1(norbs,2),stat=error)        
    if(error/=0) stop
    operanum1=0
    orbid1=0
    if(logic_bondorder/=0) then
        allocate(operanum2(nprocs-1))
        allocate(orbid2(norbs,norbs,2),stat=error)
        if(error/=0) stop
        operanum2=0   ! number of operators on specific process
        orbid2=0
    end if
    if(logic_localspin/=0) then
        allocate(operanum3(nprocs-1))
        allocate(orbid3(norbs,norbs,2),stat=error)
        if(error/=0) stop
        operanum3=0
        orbid3=0
    end if
    
    ! allocate the operators
    ! PPP operator
    do i=1,norbs,nprocs-1
        do j=1,nprocs-1,1
            if((i-1+j)<=norbs) then
                operanum1(j)=operanum1(j)+1
                orbid1((i-1)+j,1)=j
                orbid1((i-1)+j,2)=operanum1(j)
            else 
                exit
            end if
        end do
    end do 

!====================================================================
    ! bond order operator
    ! (i,i) pair (niup-nidown)^2
    !            niup-nidown
    ! (i,j) pair ai^+up*ajup
    !            ai^+down*ajdown          

    if(logic_bondorder==1) then  ! only the bondorder term is calculated
        ! L space
        do i=1,(norbs+1)/2,1
        do j=i,(norbs+1)/2,1
            if(bondlink(i,j)/=0) then
                orbid2(i,j,1)=orbid1(i,1)  ! be careful it is i here
                orbid2(j,i,1)=orbid1(i,1)
                operanum2(orbid1(i,1))=operanum2(orbid1(i,1))+1
                orbid2(i,j,2)=operanum2(orbid1(i,1))
                orbid2(j,i,2)=operanum2(orbid1(i,1))
            end if
        end do
        end do
        ! R space
        do i=(norbs+1)/2+1,norbs,1
        do j=i,norbs,1
            if(bondlink(i,j)/=0) then
                orbid2(i,j,1)=orbid1(j,1)   ! be careful it is j here
                orbid2(j,i,1)=orbid1(j,1)
                operanum2(orbid1(j,1))=operanum2(orbid1(j,1))+1
                orbid2(i,j,2)=operanum2(orbid1(j,1))
                orbid2(j,i,2)=operanum2(orbid1(j,1))
            end if
        end do
        end do
    else if(logic_bondorder==2) then ! calculate the one partical reduced density matrix
        ! L space
        do i=1,(norbs+1)/2,1
        do j=i,(norbs+1)/2,1
            orbid2(i,j,1)=orbid1(i,1)  ! be careful it is i here
            orbid2(j,i,1)=orbid1(i,1)
            operanum2(orbid1(i,1))=operanum2(orbid1(i,1))+1
            orbid2(i,j,2)=operanum2(orbid1(i,1))
            orbid2(j,i,2)=operanum2(orbid1(i,1))
        end do
        end do
        ! R space
        do i=(norbs+1)/2+1,norbs,1
        do j=i,norbs,1
            orbid2(i,j,1)=orbid1(j,1)   ! be careful it is j here
            orbid2(j,i,1)=orbid1(j,1)
            operanum2(orbid1(j,1))=operanum2(orbid1(j,1))+1
            orbid2(i,j,2)=operanum2(orbid1(j,1))
            orbid2(j,i,2)=operanum2(orbid1(j,1))
        end do
        end do
    end if
!====================================================================
    ! local spin operator
    ! the operator order is
    ! (i,i) pair ai^+down*aiup
    !            (niup-nidown)^2
    ! (i,j) pair  ai^+down*aiup*aj^+up*aj^down
    !            (niup-nidown)*(njup-njdown)
    ! L space
    if(logic_localspin/=0) then
        do i=1,(norbs+1)/2,1
        do j=i,(norbs+1)/2,1
            orbid3(i,j,1)=orbid1(i,1)
            orbid3(j,i,1)=orbid1(i,1)

            operanum3(orbid1(i,1))=operanum3(orbid1(i,1))+2

            orbid3(i,j,2)=operanum3(orbid1(i,1))
            orbid3(j,i,2)=operanum3(orbid1(i,1))
        end do
        end do
        ! R space
        do i=(norbs+1)/2+1,norbs,1
        do j=i,norbs,1
            orbid3(i,j,1)=orbid1(j,1)
            orbid3(j,i,1)=orbid1(j,1)
            
            operanum3(orbid1(j,1))=operanum3(orbid1(j,1))+2
            
            orbid3(i,j,2)=operanum3(orbid1(j,1))    ! the orbid3(:,:,2) here is the last operator index of i,j pair
            orbid3(j,i,2)=operanum3(orbid1(j,1))
        end do
        end do
    end if
!====================================================================
        
    if(myid==0) then
        write(*,*) "PPP operator distribute"
        write(*,*) "orbid1=",orbid1(:,1)
        write(*,*) "index on every process",orbid1(:,2)
        if(logic_bondorder/=0) then
            write(*,*) "bond order operator distribute"
            write(*,*) "orbid2=",orbid2(:,:,1)
        end if
        if(logic_localspin/=0) then
            write(*,*) "local spin operator distribute"
            write(*,*) "orbid3=",orbid3(:,:,1)
        end if
    end if

!============================================================================
! allocate the work space of every operator
    
    call AllocateArray

    if(myid==0) then
        if(logic_spinreversal/=0) then
            allocate(symmlinksma(subM,1,2),stat=error)
            if(error/=0) stop
            allocate(symmlinkbig(4*subM,1,2),stat=error)
            if(error/=0) stop
        end if
    end if
!------------------------------------------------------
! allocate the quanta of every many body basis
! 1 means the total electron; 2 means the total Sz
    ! diagnolization space quanta
    allocate(quantasmaL(subM,2))
    allocate(quantasmaR(subM,2))
    allocate(quantabigL(4*subM,2))
    allocate(quantabigR(4*subM,2))
    quantasmaL=0
    quantasmaR=0
    quantabigL=0
    quantabigR=0

    ! perturbation space quanta
    if(logic_perturbation/=0) then
        allocate(quantasmaLp(subMp,2))
        allocate(quantasmaRp(subMp,2))
        allocate(quantabigLp(4*subMp,2))
        allocate(quantabigRp(4*subMp,2))
        quantasmaLp=0
        quantasmaRp=0
        quantabigLp=0
        quantabigRp=0
    end if

    call Allocate_checkmem
    call Allocate_basisindex
    !------------------------------------------------------
    return
end Subroutine LoadBalance

!================================================
!================================================

subroutine Deallocate_loadbalance
    
    use module_sparse
    use variables
    use communicate

    implicit none
    deallocate(operanum1,orbid1)
    if(logic_bondorder/=0) deallocate(operanum2,orbid2)
    if(logic_localspin/=0) deallocate(operanum3,orbid3)
    if(myid==0 .and. logic_spinreversal/=0) deallocate(symmlinksma,symmlinkbig)
    deallocate(quantasmaL,quantasmaR,quantabigL,quantabigR)
    if(logic_perturbation/=0) deallocate(quantasmaLp,quantasmaRp,quantabigLp,quantabigRp)

    return
end subroutine Deallocate_loadbalance
    
!================================================
!================================================
