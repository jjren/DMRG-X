module analysis_mod
! this subroutine is to anaylsis the wavefuntion and output the analysis result
    use communicate
    implicit none

    contains
!===================================================================
!===================================================================

subroutine Analysis
    use transmoment_mod
    use bondorder_mod
    use localspin_mod
    use module_sparse
    use blas95
    use f95_precision

    implicit none
    integer(kind=i4),allocatable :: midrowmat(:),midcolmat(:)
    integer :: i

    ! in C2 mode we can calculate the two subspace together
    ! without loss of accuracy
    if(logic_C2/=0) then
        if(logic_C2==-1 .and. myid==0) then
            allocate(midrowmat(4*subM+1))
            allocate(midcolmat(coeffIFdim))
            do i=1,nstate,1
                midrowmat=coeffIFrowindex(:,i)
                coeffIFrowindex(:,i)=coeffIFrowindex(:,i+nstate)
                coeffIFrowindex(:,i+nstate)=midrowmat
                midcolmat=coeffIFcolindex(:,i)
                coeffIFcolindex(:,i)=coeffIFcolindex(:,i+nstate)
                coeffIFcolindex(:,i+nstate)=midcolmat
                call swap(coeffIF(:,i),coeffIF(:,i+nstate))
            end do
            deallocate(midrowmat,midcolmat)
        end if
        nstate=nstate*2
    end if

    call broadcastcoeffIF

    if(logic_bondorder/=0) then
        call BondOrder
    end if
    
    if(nstate/=1) then
        call TransMoment
    end if

    if(logic_localspin==1) then
        call LocalSpin
    end if
return
end subroutine Analysis

!===================================================================
!===================================================================

subroutine broadcastcoeffIF
! broadcast the coeffIF to the other process
    use variables
    use module_sparse
    use mpi
    implicit none
    character(len=1),allocatable :: packbuf(:)
    integer :: packsize,position1
    integer :: i
    integer :: error,ierr
    
    packsize=(coeffIFdim*12+(4*subM+1)*4)*nstate+1000
    allocate(packbuf(packsize),stat=error)
    if(error/=0) stop

    if(myid==0) then
        position1=0
        do i=1,nstate,1
            call MPI_PACK(coeffIFrowindex(1,i),4*subM+1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            call MPI_PACK(coeffIF(1,i),coeffIFrowindex(4*subM+1,i)-1,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
            call MPI_PACK(coeffIFcolindex(1,i),coeffIFrowindex(4*subM+1,i)-1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
        end do
    end if

    ! 0 process broadcast the coeffIF
    call MPI_BCAST(position1,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(packbuf,position1,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    
    if(myid/=0) then
        allocate(coeffIF(coeffIFdim,nstate),stat=error)
        if(error/=0) stop
        allocate(coeffIFcolindex(coeffIFdim,nstate),stat=error)
        if(error/=0) stop
        allocate(coeffIFrowindex(4*subM+1,nstate),stat=error)
        if(error/=0) stop

        position1=0
        do i=1,nstate,1
            call MPI_UNPACK(packbuf,packsize,position1,coeffIFrowindex(1,i),4*subM+1,MPI_integer4,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position1,coeffIF(1,i),coeffIFrowindex(4*subM+1,i)-1,MPI_real8,MPI_COMM_WORLD,ierr)
            call MPI_UNPACK(packbuf,packsize,position1,coeffIFcolindex(1,i),coeffIFrowindex(4*subM+1,i)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
        end do
    end if
    
    deallocate(packbuf)
return
end subroutine broadcastcoeffIF

!===================================================================
!===================================================================

end Module analysis_mod
