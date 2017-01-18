subroutine newtql2(dim1,diagonal,offdiagonal,eigenvector,ierr)
    USE LAPACK95
    USE F95_PRECISION
    implicit none

    integer :: dim1,ierr
    real(kind=8) :: diagonal(dim1),offdiagonal(dim1-1),eigenvector(dim1,dim1)

    call rsteqr(diagonal,offdiagonal,eigenvector,"I",ierr)

    return
end subroutine newtql2
