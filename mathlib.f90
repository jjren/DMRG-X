module MathLib
    
    use kinds_mod
    use blas95
    use lapack95
    use f95_precision

    implicit none
contains

!=============================================================================
!=============================================================================

subroutine Mataxpy(a,nrow,ncol,matx,maty)
! Maty=a*matx+maty
    implicit none
    integer,intent(in) :: nrow,ncol
    real(kind=r8),intent(in) :: matx(nrow,ncol),a
    real(kind=r8),intent(inout) :: maty(nrow,ncol)
    
    ! local
    integer :: icol

    do icol=1,ncol,1
        call axpy(matx(1:nrow,icol),maty(1:nrow,icol),a)
    end do

    return
end subroutine Mataxpy
        
!=============================================================================
!=============================================================================
Subroutine DirectProduct(a,dima,b,dimb,c)
! this subroutine is to do direct product of two matrix
! for example :<sigmaL|<L|O|L>|sigmaL>
! a is the L matrix, b is the sigmaL matrix

implicit none

integer(kind=i4),intent(in) :: dima,dimb
real(kind=r8),intent(in) :: a(dima,dima),b(dimb,dimb)
real(kind=r8),intent(out) :: c(dima*dimb,dima*dimb)
integer :: iar,ial,ibr,ibl

do ibr=1,dimb,1
do ibl=1,dimb,1
    ! out-place copy
    call mkl_domatcopy('C','N',dima,dima,b(ibl,ibr),a,dima,c((ibl-1)*dima+1:ibl*dima,(ibr-1)*dima+1:ibr*dima),dima)
end do
end do

return
end Subroutine DirectProduct

!===============================================
!===============================================

subroutine SparseDirectProduct( ancols,anrows,amat,acolindex,arowindex,&
                                bncols,bnrows,bmat,bcolindex,browindex,&
                                cmat,ccolindex,crowindex,maxnelement)
! this subroutine is do Sparse Matrix direct product
! in 3-array CSR form
! mata direct product matb

    implicit none
    
    integer(kind=i4) :: ancols,anrows,bncols,bnrows
    integer(kind=i4) :: arowindex(anrows+1),browindex(bnrows+1)
    integer(kind=i4) :: acolindex(arowindex(anrows+1)-1),bcolindex(browindex(bnrows+1)-1)
    real(kind=r8) :: amat(arowindex(anrows+1)-1),bmat(browindex(bnrows+1)-1)
    
    integer :: maxnelement
    integer(kind=i4) :: crowindex(anrows*bnrows+1)
    integer(kind=i4) :: ccolindex(maxnelement)
    real(kind=r8) :: cmat(maxnelement)
    
    ! local
    integer :: ai,aj,bi,bj
    integer :: cbegin,nonzero
    
    !cmat=0.0D0
    !crowindex=0
    !ccolindex=0

    nonzero=0
    crowindex(1)=1
    do bi=1,bnrows,1
        do ai=1,anrows,1
            do bj=browindex(bi),browindex(bi+1)-1,1
                cbegin=(bcolindex(bj)-1)*ancols
                do aj=arowindex(ai),arowindex(ai+1)-1,1
                    nonzero=nonzero+1
                    if(nonzero>maxnelement) then
                        write(*,*) "===================================================="
                        write(*,*) "failed! in sparsedirectproduct nonzero>maxnelement",nonzero,maxnelement
                        write(*,*) "===================================================="
                        stop
                    end if
                    cmat(nonzero)=amat(aj)*bmat(bj)
                    ccolindex(nonzero)=cbegin+acolindex(aj)
                end do
            end do
            crowindex((bi-1)*anrows+ai+1)=nonzero+1
        end do
    end do

    return

end subroutine SparseDirectProduct

!===============================================
!===============================================

subroutine GramSchmit(nvector,lvector,vectorwork,normwork)
! Gram-Schmit Orthogonalization subroutine

! input - nvector :: the number of vectors
! input - lvector :: the length of every vector
! input - vectorwork :: the workmemory of vectors
! output - normwork :: the norm of every vector after orthogonalization


    implicit none
    integer :: nvector,lvector
    real(kind=r8) :: vectorwork(lvector*nvector),normwork(nvector)&
    ,overlap
    integer ::  i,j
    
    do i=1,nvector,1
        do j=1,i-1,1
            overlap=dot(vectorwork((i-1)*lvector+1:i*lvector),&
                vectorwork((j-1)*lvector+1:j*lvector))

            vectorwork((i-1)*lvector+1:i*lvector)=&
                vectorwork((i-1)*lvector+1:i*lvector)-&
                overlap*vectorwork((j-1)*lvector+1:j*lvector)
        end do
            normwork(i)=dot(vectorwork((i-1)*lvector+1:i*lvector),&
                vectorwork((i-1)*lvector+1:i*lvector))
            if(normwork(i)<1.0D-10) then
                write(*,*) "--------------------------"
                write(*,*) "norm is < 1.0D-10,caution!",i,"th state"
                write(*,*) "--------------------------"
            end if
            vectorwork((i-1)*lvector+1:i*lvector)=&
                vectorwork((i-1)*lvector+1:i*lvector)/sqrt(normwork(i))
    end do
    return

end subroutine GramSchmit

!=============================================

subroutine ScaleMatrix(mat,rows,cols,scaling,op)
! this subroutine use mkl mkl_?imatcopy to do scaling a matrix

    implicit none
    include "mkl_trans.fi"

    integer :: rows,cols
    real(kind=r8) :: mat(rows,cols),scaling
    character(len=1) :: op
    integer :: src_lda,dst_lda
    
    src_lda=rows
    if(op=='T' .or. op=='t' .or. op=='C' .or. op=='c') then
        dst_lda=cols
    else if(op=='N' .or. op=='n' .or. op=='R' .or. op=='r') then
        dst_lda=rows
    end if

    call mkl_dimatcopy('C',op,rows,cols,scaling,mat,src_lda,dst_lda)

return
end subroutine  ScaleMatrix

!=============================================
!=============================================

subroutine Diagsyev(dim1,mat,eigenvalue,eigenvector)
    implicit none
    
    integer :: dim1
    real(kind=r8) :: mat(dim1,dim1),eigenvalue(dim1),eigenvector(dim1,dim1)

    call syevr(mat,eigenvalue,'U',eigenvector)

return

end subroutine diagsyev

!=============================================
!=============================================

subroutine SpMatRotateBasis(Uncols,Unrows,Umat,Ucolindex,Urowindex,&
                    Oncols,Onrows,Omat,Ocolindex,Orowindex,&
                    Anrows,Amat,Acolindex,Arowindex,maxnelement)
! this subroutine is to adapted operator matrix to new basis
! UOU+ to A
! U is a sparse matrix
! O is a sparse matrix
    
    implicit none
    include "mkl_spblas.fi"

    integer(kind=i4) :: Uncols,Unrows,Oncols,Onrows,Anrows
    integer(kind=i4) :: maxnelement
    integer(kind=i4) :: Urowindex(Unrows+1),Orowindex(Onrows+1),Arowindex(Anrows+1)
    integer(kind=i4) :: &
    Ucolindex(Urowindex(Unrows+1)-1) , &
    Ocolindex(Orowindex(Onrows+1)-1) , &
    Acolindex(maxnelement)
    real(kind=r8) :: &
    Umat(Urowindex(Unrows+1)-1) , &
    Omat(Orowindex(Onrows+1)-1) , &
    Amat(maxnelement)

    integer(kind=i4),allocatable :: bufcolindex(:),bufrowindex(:)
    real(kind=r8),allocatable :: bufmat(:)
    integer :: sort,bufnrows
    integer :: error,info
    
    if(Anrows/=Uncols ) then
        write(*,*) "============="
        write(*,*) "Anrows/Uncols",Anrows,Uncols
        write(*,*) "============="
        stop
    end if
    if(Unrows/=Onrows) then
        write(*,*) "============="
        write(*,*) "Unrows/Onrows",Unrows,Onrows
        write(*,*) "============="
        stop
    end if
    bufnrows=Uncols

    allocate(bufrowindex(bufnrows+1),stat=error)
    if(error/=0) stop
    allocate(bufcolindex(Uncols*Unrows),stat=error)
    if(error/=0) stop
    allocate(bufmat(Uncols*Unrows),stat=error)
    if(error/=0) stop

    sort=8
    call mkl_dcsrmultcsr('T',0,sort,Unrows,Uncols,Oncols, &
        Umat,Ucolindex,Urowindex, &
        Omat,Ocolindex,Orowindex, &
        bufmat,bufcolindex,bufrowindex,Uncols*Unrows,info)
        
    if(info/=0) then
        write(*,*) "==============================================="
        write(*,*) "info/=0,spmatrotatebasis fist step failed!",info
        write(*,*) "==============================================="
        stop
    end if


    call mkl_dcsrmultcsr('N',0,sort,bufnrows,Unrows,Uncols, &
        bufmat,bufcolindex,bufrowindex , &
        Umat,Ucolindex,Urowindex, &
        Amat,Acolindex,Arowindex,maxnelement,info)

    if(info/=0) then
        write(*,*) "================================================="
        write(*,*) "info/=0 spmatrotatebasis second step failed!",info
        write(*,*) "================================================="
        stop
    end if

    deallocate(bufmat)
    deallocate(bufcolindex)
    deallocate(bufrowindex)
return

end subroutine SpMatRotateBasis

!=============================================
!=============================================

subroutine CSCtoCSR(operation,nleadaft,nleadbef,Amat,Acolindex,Arowindex,Arowindexnew)
! convert sparse matrix CSC format to CSR format; vice versa
! dcsrcsc only support square matrix
! CSCtoCSR have no such limitation
! nleadaft is the CSR ncols; CSC nrows
! nleadbef is the CSR nrows; CSC ncols
! Arowindex is sum 
! Acolindex is index
! in fact CSC to CSR and CSR to CSC is the same ; They both are tranpose

    implicit none
    include "mkl_spblas.fi"

    character(len=2) :: operation
    integer(kind=i4) :: nleadaft,nleadbef
    integer(kind=i4) :: &
    Arowindex(nleadbef+1) , &
    Arowindexnew(nleadaft+1) , &
    Acolindex(Arowindex(nleadbef+1)-1)
    real(kind=r8) :: Amat(Arowindex(nleadbef+1)-1)
    ! local
    integer :: job(8),dim1,info
    integer(kind=i4),allocatable :: Browindex(:),Crowindex(:),Ccolindex(:)
    real(kind=r8),allocatable :: Cmat(:)
    integer :: error
    integer :: i,j

    if(nleadaft>nleadbef) then
        dim1=nleadaft
    else
        dim1=nleadbef
    end if

    allocate(Browindex(dim1+1),stat=error)
    if(error/=0) stop
    allocate(Crowindex(dim1+1),stat=error)
    if(error/=0) stop
    allocate(Ccolindex(Arowindex(nleadbef+1)-1),stat=error)
    if(error/=0) stop
    allocate(Cmat(Arowindex(nleadbef+1)-1),stat=error)
    if(error/=0) stop

    Browindex(1:nleadbef+1)=Arowindex(1:nleadbef+1)
    if(nleadaft>nleadbef) then
        Browindex(nleadbef+2:dim1+1)=Arowindex(nleadbef+1)
    end if

    job(2)=1
    job(3)=1
    job(6)=1
    if(operation=='CR' .or. operation=='cr') then
        job(1)=1
        call mkl_dcsrcsc(job,dim1,Cmat,Ccolindex,Crowindex,Amat,Acolindex,Browindex,info)
    else if(operation=='RC' .or. operation=='rc') then
        job(1)=0
        call mkl_dcsrcsc(job,dim1,Amat,Acolindex,Browindex,Cmat,Ccolindex,Crowindex,info)
    end if

    if(info/=0) then
        write(*,*) "========================================"
        write(*,*) "info/=0 subroutine CSRtoCSC failed!",info
        write(*,*) "========================================"
        stop
    end if

    call CopySpAtoB(nleadaft,Cmat,Ccolindex,Crowindex,Amat,Acolindex,&
            Arowindexnew,Crowindex(nleadaft+1)-1)

    deallocate(Browindex)
    deallocate(Crowindex)
    deallocate(Ccolindex)
    deallocate(Cmat)
return

end subroutine CSCtoCSR

!=============================================
!=============================================

subroutine SpMatAdd(Ancols,Anrows,Amat,Amatcol,Amatrow,&
        trans,beta,Bncols,Bnrows,Bmat,Bmatcol,Bmatrow,maxnelement)
! this subroutine do sparse matrix add A+beta*B=A
    implicit none
    include "mkl_spblas.fi"

    character(len=1) :: trans
    integer(kind=i4) :: Ancols,Anrows,Bncols,Bnrows
    integer(kind=i4) :: maxnelement
    integer(kind=i4) :: Amatrow(Anrows+1),Bmatrow(Bnrows+1)
    integer(kind=i4) :: Amatcol(maxnelement),Bmatcol(Bmatrow(Bnrows+1)-1)
    real(kind=r8) :: Amat(maxnelement),Bmat(Bmatrow(Bnrows+1)-1)
    real(kind=r8) :: beta
    
    integer(kind=i4),allocatable :: Cmatcol(:),Cmatrow(:),Dmatrow(:)
    real(kind=r8),allocatable :: Cmat(:)
    integer :: error,info
    real(kind=r8) :: Dmat(1)
    integer(kind=i4) :: Dmatcol(1)
    
    allocate(Cmat(maxnelement),stat=error)
    if(error/=0) stop
    allocate(Cmatcol(maxnelement),stat=error)
    if(error/=0) stop
    allocate(Cmatrow(Anrows+1),stat=error)
    if(error/=0) stop

    call mkl_dcsradd(trans,0,0,Anrows,Ancols,Amat,Amatcol,Amatrow,&
    beta,Bmat,Bmatcol,Bmatrow,Cmat,Cmatcol,Cmatrow,maxnelement,info) 
    
    if(info/=0) then
        write(*,*) "==============================" 
        write(*,*) "SpMatAdd info/=0 failed!",info
        write(*,*) "=============================="
        stop
    end if
    
    allocate(Dmatrow(Anrows+1))
    ! zero sparse matrix intialize
    Dmatrow=1
    
    call mkl_dcsradd('N',0,0,Anrows,Ancols,Cmat,Cmatcol,Cmatrow,&
        0.0D0,Dmat,Dmatcol,Dmatrow,Amat,Amatcol,Amatrow,maxnelement,info) 
    
    if(info/=0) then
        write(*,*) "==============================" 
        write(*,*) "SpMatAdd info/=0 failed!",info
        write(*,*) "=============================="
        stop
    end if
    
    deallocate(Cmat,Cmatcol,Cmatrow,Dmatrow)

    return
end subroutine SpMatAdd
    
!=============================================
!=============================================

subroutine SpMatIJ(Anrows,irow,icol,Amat,Amatcol,Amatrow,output)
! this subroutine is to get the sparse matrix A's A(irow,icol) to ouput

    implicit none

    integer(kind=i4) :: Anrows,icol,irow
    integer(kind=i4) :: Amatrow(Anrows+1),Amatcol(Amatrow(Anrows+1)-1)
    real(kind=i8) :: Amat(Amatrow(Anrows+1)-1)
    real(kind=i8) :: output
    ! local
    integer :: i,j
    integer :: num,index1

    num=Amatrow(irow+1)-Amatrow(irow)
    output=0.0D0

    call BinarySearch(num,Amatcol(Amatrow(irow)),icol,index1)
    if(index1/=0) then
        output=Amat(index1+Amatrow(irow)-1)
    end if
    
    return

end subroutine SpMatIJ

!=============================================
!=============================================

subroutine SpMMtoDens(transA,transB,Anrows,Ancols,Bnrows,Bncols,Cnrows,Cncols ,&
    Amat,Acolindex,Arowindex,Bmat,Bcolindex,Browindex,Cmat,ldc)
! this subroutine do sparse matrix matrix product the output is a dense matrix
    implicit none
    include "mkl_spblas.fi"
    character(len=1) :: transA,transB
    integer(kind=i4) :: Anrows,Ancols,Bnrows,Bncols,Cnrows,Cncols,ldc
    integer(kind=i4) :: Arowindex(Anrows+1),Browindex(Bnrows+1)
    integer(kind=i4) :: Acolindex(Arowindex(Anrows+1)-1),Bcolindex(Browindex(Bnrows+1)-1)
    real(kind=r8) :: Amat(Arowindex(Anrows+1)-1),Bmat(Browindex(Bnrows+1)-1),Cmat(Cnrows,Cncols)
    
    ! local
    integer(kind=i4),allocatable :: Bcolindexdummy(:),Browindexdummy(:),Ccolindexdummy(:),Crowindexdummy(:)
    real(kind=r8),allocatable :: Bmatdummy(:),Cmatdummy
    integer :: error


!   if(transB=='T') then
!       allocate(Bcolindexdummy(Browindex(Bnrows+1)-1),stat=error)
!       if(error/=0) stop
!       allocate(Bmatdummy(Browindex(Bnrows+1)-1),stat=error)
!       if(error/=0) stop
!       allocate(Browindexdummy(Bncols+1),stat=error)
!       if(error/=0) stop
!       Bmatdummy=Bmat
!       Bcolindexdummy=Bcolindex
!
!       call CSCtoCSR('RC',Bncols,Bnrows,Bmatdummy,Bcolindexdummy,Browindex,Browindexdummy)
!       call mkl_dcsrmultcsr('N',0,8,,4*Lrealdim,4*Rrealdim, &
!               Amat,Acolindex,Arowindex,&
!               Bmatdummy,Bcolindexdummy,Browindexdummy,&
!               buffmat,buffmatcol,buffmatrow,LRoutnelement,info)
!   !   call mkl_dcsrmultd(transA,Anrows,Ancols,Bnrows,Amat,Acolindex,Arowindex,Bmatdummy,Bcolindexdummy,Browindexdummy,Cmat,ldc)
!       deallocate(Browindexdummy,Bcolindexdummy,Bmatdummy)
!   else
!   
!   !   call mkl_dcsrmultd(transA,Anrows,Ancols,Bncols,Amat,Acolindex,Arowindex,Bmatdummy,Bcolindexdummy,Browindexdummy,Cmat,ldc)
!   end if

return

end subroutine SpMMtoDens

!=============================================
!=============================================

subroutine SpMMtoSp(transA,transB,Anrows,Ancols,Bnrows,Bncols,Cnrows,&
    Amat,Acolindex,Arowindex,Bmat,Bcolindex,Browindex,Cmat,Ccolindex,Crowindex,maxnelement)
! this subroutine do sparse matrix matrix product the output is a sparse matrix
    implicit none
    include "mkl_spblas.fi"

    character(len=1) :: transA,transB
    integer(kind=i4) :: Anrows,Ancols,Bnrows,Bncols,maxnelement,Cnrows
    integer(kind=i4) :: Arowindex(Anrows+1),Browindex(Bnrows+1),Crowindex(Cnrows+1)
    integer(kind=i4) :: Acolindex(Arowindex(Anrows+1)-1),Bcolindex(Browindex(Bnrows+1)-1),Ccolindex(maxnelement)
    real(kind=r8) :: Amat(Arowindex(Anrows+1)-1),Bmat(Browindex(Bnrows+1)-1),Cmat(maxnelement)
    
    ! local
    integer(kind=i4),allocatable :: Bcolindexdummy(:),Browindexdummy(:)
    real(kind=r8),allocatable :: Bmatdummy(:)
    integer :: error,info

    if(transB=='T') then
        allocate(Bcolindexdummy(Browindex(Bnrows+1)-1),stat=error)
        if(error/=0) stop
        allocate(Bmatdummy(Browindex(Bnrows+1)-1),stat=error)
        if(error/=0) stop
        allocate(Browindexdummy(Bncols+1),stat=error)
        if(error/=0) stop

        ! copy the Bmat to another work array
        call copy(Bmat,Bmatdummy)
        Bcolindexdummy=Bcolindex

        call CSCtoCSR('RC',Bncols,Bnrows,Bmatdummy,Bcolindexdummy,Browindex,Browindexdummy)
        call mkl_dcsrmultcsr(transA,0,8,Anrows,Ancols,Bnrows, &
                Amat,Acolindex,Arowindex,&
                Bmatdummy,Bcolindexdummy,Browindexdummy,&
                Cmat,Ccolindex,Crowindex,maxnelement,info)
        if(info/=0) then
            write(*,*) "===================="
            write(*,*) "SpMMtoSp 1st failed!",info
            write(*,*) "===================="
            stop
        end if
        deallocate(Browindexdummy,Bcolindexdummy,Bmatdummy)
    else
        call mkl_dcsrmultcsr(transA,0,8,Anrows,Ancols,Bncols, &
                Amat,Acolindex,Arowindex,&
                Bmat,Bcolindex,Browindex,&
                Cmat,Ccolindex,Crowindex,maxnelement,info)
        if(info/=0) then
            write(*,*) "===================="
            write(*,*) "SpMMtoSp 2st failed!",info
            write(*,*) "===================="
            stop
        end if
    end if

return

end subroutine SpMMtoSp

!=============================================
!=============================================

subroutine SpMMtoSptodens(transA,transB,Anrows,Ancols,Bnrows,Bncols,Cnrows,Cncols,&
    Amat,Acolindex,Arowindex,Bmat,Bcolindex,Browindex,Cmat)
    implicit none
    include "mkl_spblas.fi"
    character(len=1) :: transA,transB
    integer(kind=i4) :: Anrows,Ancols,Bnrows,Bncols,Cnrows,Cncols
    integer(kind=i4) :: Arowindex(Anrows+1),Browindex(Bnrows+1)
    integer(kind=i4) :: Acolindex(Arowindex(Anrows+1)-1),Bcolindex(Browindex(Bnrows+1)-1)
    real(kind=r8) :: Amat(Arowindex(Anrows+1)-1),Bmat(Browindex(Bnrows+1)-1),Cmat(Cnrows,Cncols)
    
    ! local
    integer(kind=i4),allocatable :: Ccolindexdummy(:),Crowindexdummy(:)
    real(kind=r8),allocatable :: Cmatdummy(:)
    integer :: maxnelement,job(8)
    integer :: error,info
    
    maxnelement=Cnrows*Cncols

    allocate(Cmatdummy(maxnelement),stat=error)
    if(error/=0) stop
    allocate(Ccolindexdummy(maxnelement),stat=error)
    if(error/=0) stop
    allocate(Crowindexdummy(Cnrows+1),stat=error)
    if(error/=0) stop
    
    ! C is in CSR format
    call SpMMtoSp(transA,transB,Anrows,Ancols,Bnrows,Bncols,Cnrows,&
    Amat,Acolindex,Arowindex,Bmat,Bcolindex,Browindex,cmatdummy,ccolindexdummy,crowindexdummy,maxnelement)

    ! transfer C to the dense format
    job(1)=1
    job(2)=1
    job(3)=1
    job(4)=2
    job(5)=0
    job(6)=1
    call mkl_ddnscsr(job,Cnrows,Cncols,Cmat,Cnrows,cmatdummy,ccolindexdummy,crowindexdummy,info)
    if(info/=0) then
        write(*,*) "======================"
        write(*,*) "SpMMtoSptodens failed!"
        write(*,*) "======================"
        stop
    end if

    deallocate(Cmatdummy,Ccolindexdummy,Crowindexdummy)
return

end subroutine SpMMtoSptodens

!=============================================
!=============================================

subroutine SpMMtrace(trans,Anrows,Ancols,Amat,Acolindex,Arowindex,Bnrows,Bncols,Bmat,Bcolindex,Browindex,output)
! this subroutine do square Matrix Matrix Product trace in CSR format
! if trans=='T' then Bmat need not transpose
! if trans=='N' then Bmat need tranpose
    implicit none
    
    character(len=1),intent(in) ::  trans
    integer,intent(in) :: Anrows,Bnrows,Ancols,Bncols
    integer(kind=i4),intent(in) :: Arowindex(Anrows+1),Browindex(Bnrows+1),&
        Acolindex(Arowindex(Anrows+1)-1),Bcolindex(Browindex(Bnrows+1)-1)
    real(kind=r8) :: Amat(Arowindex(Anrows+1)-1),Bmat(Browindex(Bnrows+1)-1)
    real(kind=r8) :: output

    ! local
    integer(kind=i4),allocatable :: Bcolindexdummy(:),Browindexdummy(:)
    real(kind=r8),allocatable :: Bmatdummy(:)
    integer :: i,j,index1
    real(kind=r8) :: dotout

    if(trans=="T") then
        if(Anrows/=Bnrows .or. Ancols/=Bncols) then
            write(*,*) "Anrows/=Bnrows .or. Ancols/=Bncols",Anrows,Bnrows,Ancols,Bncols
            stop
        end if
    else if(trans=="N") then
        if(Anrows/=Bncols .or. Ancols/=Bnrows) then
            write(*,*) "Anrows/=Bncols .or. Ancols/=Bnrows",Anrows,Bncols,Ancols,Bnrows
            stop
        end if
    end if

    allocate(Browindexdummy(Anrows+1))
    allocate(Bcolindexdummy(Browindex(Bnrows+1)-1))
    allocate(Bmatdummy(Browindex(Bnrows+1)-1))
    
    ! copy the Bmat to another workarray
    call CopySpAtoB(Bnrows,Bmat,Bcolindex,Browindex,Bmatdummy,Bcolindexdummy,Browindexdummy,Browindex(Bnrows+1)-1)
    
    if(trans=='N') then
        call CSCtoCSR('RC',Bncols,Bnrows,Bmatdummy,Bcolindexdummy,Browindex,Browindexdummy)
    end if

    output=0.0D0
    do i=1,Anrows,1
    !    do j=Arowindex(i),Arowindex(i+1)-1,1
    !        call BinarySearch(Browindexdummy(i+1)-Browindexdummy(i),Bcolindexdummy(Browindexdummy(i)),Acolindex(j),index1)
    !        if(index1/=0) then
    !            output=output+Bmatdummy(index1+Browindexdummy(i)-1)*Amat(j)
    !        end if
    !    end do
        call SpDotVec(Arowindex(i+1)-Arowindex(i),Amat(Arowindex(i)),Acolindex(Arowindex(i)),&
                Browindexdummy(i+1)-Browindexdummy(i),Bmatdummy(Browindexdummy(i)),&
                Bcolindexdummy(Browindexdummy(i)),dotout)
        output=dotout+output
    end do
    
    deallocate(Bmatdummy,Bcolindexdummy,Browindexdummy)
return

end subroutine SpMMtrace

!=============================================
!=============================================

subroutine GetSpMat1Vec(dim1,index1,first,last,Amat,Acol,Arow,direction,vecdim,vector)
! this subroutine get the one row or one column of a sparse matrix
! dim1 :: the dimensiton of matrix A
! index1 :: the specific col/row index
! first/last :: the first/last element index of such row and col
! Amat,Acol,Arow
! direction :: row/col
! vecdim,vector :: store vector

    implicit none

    integer :: dim1,vecdim,index1,first,last
    character(len=3) :: direction
    integer(kind=i4) :: Arow(dim1+1),Acol(Arow(dim1+1)-1)
    real(kind=r8) :: Amat(Arow(dim1+1)-1),vector(vecdim)
    
    integer :: i,j
    
    if(vecdim/=last-first+1) then
        write(*,*) "GetSpmat1Vec dim1/=last-first+1",dim1,last,first
        stop
    end if

    vector=0.0D0
    if(direction=="col") then
        do i=first,last,1
            do j=Arow(i),Arow(i+1)-1,1
                if(Acol(j)>index1) then
                    exit
                else if(Acol(j)==index1) then
                    vector(i-first+1)=Amat(j)
                    exit
                end if
            end do
        end do
    else if(direction=="row") then
        do i=Arow(index1),Arow(index1+1)-1,1
            if(Acol(i)>=first .and. Acol(i)<=last) then
                vector(Acol(i)-first+1)=Amat(i)
            end if
        end do
    else
        write(*,*) "direction=",direction,"failed!"
    end if
return

end subroutine GetSpmat1Vec

!=============================================
!=============================================

subroutine CopySpAtoB(nrow,Amat,Amatcol,Amatrow,Bmat,Bmatcol,Bmatrow,Bmaxnelement)
    ! this subroutine do sparse matrix A to B copy
    implicit none
    integer,intent(in) :: nrow,Bmaxnelement
    integer(kind=i4),intent(in) :: Amatrow(nrow+1),Amatcol(Amatrow(nrow+1)-1)
    real(kind=r8),intent(in) :: Amat(Amatrow(nrow+1)-1)
    integer(kind=i4),intent(out) :: Bmatrow(nrow+1),Bmatcol(Bmaxnelement)
    real(kind=r8),intent(out) :: Bmat(Bmaxnelement)

    if(Amatrow(nrow+1)-1>Bmaxnelement) then
        write(*,*) "========================="
        write(*,*) "Bmaxnelement is too small"
        write(*,*) "Bmaxnelement=",Bmaxnelement
        write(*,*) "Amatrow(nrow+1)-1=",Amatrow(nrow+1)-1
        write(*,*) "========================="
        stop
    end if
    
    !Bmatrow(1:nrow+1)=Amatrow(1:nrow+1)
    !Bmatcol(1:Amatrow(nrow+1)-1)=Bmatcol(1:Amatrow(nrow+1)-1)
    call scopy(nrow+1,Amatrow,1,Bmatrow,1)
    call scopy(Amatrow(nrow+1)-1,Amatcol,1,Bmatcol,1)
    call copy(Amat(1:Amatrow(nrow+1)-1),Bmat(1:Amatrow(nrow+1)-1))

return
end subroutine CopySpAtoB

!=============================================
!=============================================

subroutine BinarySearch(num,keyvalue,targetvalue,index1)
    implicit none
    integer,intent(in) :: num,keyvalue(num),targetvalue
    integer,intent(out) :: index1
    ! local
    integer :: firstindex,lastindex,testindex

    index1=0
    if(num<=0 .or. &
        keyvalue(1)>targetvalue .or. &
        keyvalue(num)<targetvalue) return
    
    firstindex=1
    lastindex=num

    if(keyvalue(firstindex)==targetvalue) then
        index1=firstindex
    else if(keyvalue(lastindex)==targetvalue) then
        index1=lastindex
    end if

    testindex=(firstindex+lastindex)/2
    do while(testindex/=firstindex)
        if(keyvalue(testindex)<targetvalue) then
            firstindex=testindex
        else if(keyvalue(testindex)>targetvalue) then
            lastindex=testindex
        else if(keyvalue(testindex)==targetvalue) then
            index1=testindex
            exit
        end if
        testindex=(firstindex+lastindex)/2
    end do
    
    return

end subroutine BinarySearch

!=============================================
!=============================================

subroutine SpDotVec(An,Avec,Acol,Bn,Bvec,Bcol,output)
! calculate sparse vector dot product
    implicit none
    integer,intent(in) :: An,Bn
    real(kind=r8),intent(in) :: Avec(An),Bvec(Bn)
    integer(kind=i4),intent(in) :: Acol(An),Bcol(Bn)
    real(kind=r8),intent(out) :: output

    ! local
    integer :: i,index1

    output=0.0D0

    if(An<=Bn) then
        do i=1,An,1
            call BinarySearch(Bn,Bcol,Acol(i),index1)
            if(index1/=0) then
                output=output+Bvec(index1)*Avec(i)
            end if
        end do
    else
        do i=1,Bn,1
            call BinarySearch(An,Acol,Bcol(i),index1)
            if(index1/=0) then
                output=output+Avec(index1)*Bvec(i)
            end if
        end do
    end if
    
    return

end subroutine SpDotVec

!=============================================
!=============================================
end module MathLib
