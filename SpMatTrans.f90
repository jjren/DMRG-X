Module SpMatTrans_mod
    USE MPI
    use kinds_mod
    implicit none

contains
!===========================================
!===========================================

subroutine SpMatPack(mat,matcol,matrow,nrows,position1,packbuf,packsize)
    
    implicit none
    integer(kind=i4),intent(in) :: nrows,packsize,&
        matrow(nrows+1),&
        matcol(matrow(nrows+1)-1)
    real(kind=r8),intent(in) :: mat(matrow(nrows+1)-1)
    integer(kind=i4),intent(inout) :: position1
    character(len=1),intent(inout) :: packbuf(:)
    ! local
    integer :: ierr

    call MPI_PACK(matrow,nrows+1,MPI_integer4,&
        packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
    call MPI_PACK(mat(1),matrow(nrows+1)-1,MPI_real8,&
        packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
    call MPI_PACK(matcol(1),matrow(nrows+1)-1,MPI_integer4,&
        packbuf,packsize,position1,MPI_COMM_WORLD,ierr)

    return
end subroutine SpMatPack

!===========================================
!===========================================

subroutine SpMatUnPack(mat,matcol,matrow,nrows,maxnelement,position1,packbuf,packsize)
    
    implicit none
    integer(kind=i4),intent(in) :: nrows,packsize,maxnelement
    character(len=1),intent(in) :: packbuf(:)
    
    integer(kind=i4),intent(out) :: matrow(nrows+1),matcol(:)
    real(kind=r8),intent(out) :: mat(:)
    integer(kind=i4),intent(inout) :: position1
    ! local
    integer :: ierr
    
    call MPI_UNPACK(packbuf,packsize,position1,&
        matrow,nrows+1,MPI_integer4,MPI_COMM_WORLD,ierr)
    if(matrow(nrows+1)-1>maxnelement) then
        write(*,*) "SpMatUnPack matrow(nrows+1)-1>maxnelement"
        write(*,*) matrow(nrows+1)-1,maxnelement
        stop
    end if
    call MPI_UNPACK(packbuf,packsize,position1,&
        mat(1),matrow(nrows+1)-1,MPI_real8,MPI_COMM_WORLD,ierr)
    call MPI_UNPACK(packbuf,packsize,position1,&
        matcol(1),matrow(nrows+1)-1,MPI_integer4,MPI_COMM_WORLD,ierr)
    return
end subroutine SpMatUnPack

!===========================================
!===========================================
end Module SpMatTrans_mod
