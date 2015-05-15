module MathLib
	
	use kinds_mod
	use blas95
	use lapack95
	use f95_precision

	implicit none
	contains
!=============================================================================
!=============================================================================
Subroutine DirectProduct(a,dima,b,dimb,c)
! this subroutine is to do direct product of two matrix
! for example :<sigmaL|<L|O|L>|sigmaL>
! a is the L matrix, b is the sigmaL matrix

implicit none

integer(kind=i4) :: dima,dimb
real(kind=r8) ::a(dima,dima),b(dimb,dimb),c(dima*dimb,dima*dimb)
integer :: iar,ial,ibr,ibl

c=0.0D0
do ibr=1,dimb,1
do iar=1,dima,1
	do ibl=1,dimb,1
	do ial=1,dima,1
		c((ibl-1)*dima+ial,(ibr-1)*dima+iar)=a(ial,iar)*b(ibl,ibr)
	end do
	end do
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
	
	cmat=0.0D0
	crowindex=0
	ccolindex=0

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

!	write(*,*) "nonzero=",nonzero
!	write(*,*) "maxnelement=",maxnelement

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

subroutine spmatrotatebasis(Uncols,Unrows,Umat,Ucolindex,Urowindex,&
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

	if(info==0) then
	!	write(*,*) "maxnelement=",maxnelement
	!	write(*,*) "Arowindex(Anrows+1)-1=",Arowindex(Anrows+1)-1
	else
		write(*,*) "================================================="
		write(*,*) "info/=0 spmatrotatebasis second step failed!",info
		write(*,*) "================================================="
		stop
	end if

	deallocate(bufmat)
	deallocate(bufcolindex)
	deallocate(bufrowindex)
return

end subroutine spmatrotatebasis

!=============================================
!=============================================

subroutine CSCtoCSR(operation,ncols,nrows,Amat,Acolindex,Arowindex,Arowindexnew)
! convert sparse matrix CSC format to CSR format; vice versa
! dcsrcsc only support square matrix
! CSCtoCSR have no such limitation
	implicit none
	include 'mkl_spblas.fi'

	character(len=2) :: operation
	integer(kind=i4) :: ncols,nrows
	integer(kind=i4) :: &
	Arowindex(nrows+1) , &
	Arowindexnew(ncols+1) , &
	Acolindex(Arowindex(nrows+1)-1)
	real(kind=r8) :: Amat(Arowindex(nrows+1)-1)
	! local
	integer :: job(8),dim1,info
	integer(kind=i4),allocatable :: Browindex(:),Crowindex(:),Ccolindex(:)
	real(kind=r8),allocatable :: Cmat(:)
	integer :: error
	integer :: i,j

	if(ncols>nrows) then
		dim1=ncols
	else
		dim1=nrows
	end if

	allocate(Browindex(dim1+1),stat=error)
	if(error/=0) stop
	allocate(Crowindex(dim1+1),stat=error)
	if(error/=0) stop
	allocate(Ccolindex(Arowindex(nrows+1)-1),stat=error)
	if(error/=0) stop
	allocate(Cmat(Arowindex(nrows+1)-1),stat=error)
	if(error/=0) stop

	Browindex(1:nrows+1)=Arowindex(1:nrows+1)
	if(ncols>nrows) then
		Browindex(nrows+2:dim1+1)=Arowindex(nrows+1)
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
	
!	open(unit=11,file="test.tmp",status="replace")
!	do i=1,nrows,1
!		do j=Arowindex(i),Arowindex(i+1)-1,1
!			write(11,*) Amat(j),Acolindex(j),i
!		end do
!	end do
!	do i=1,ncols,1
!		do j=Crowindex(i),Crowindex(i+1)-1,1
!			write(11,*) Cmat(j),Ccolindex(j),i
!		end do
!	end do
!	close(11)


	Amat=Cmat
	Acolindex=Ccolindex
	Arowindexnew=Crowindex(1:ncols+1)

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
	
	integer(kind=i4),allocatable :: Cmatcol(:),Cmatrow(:)
	real(kind=r8),allocatable :: Cmat(:)
	integer :: error,info

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
	Amat=Cmat
	Amatcol=Cmatcol
	Amatrow=Cmatrow

	deallocate(Cmat,Cmatcol,Cmatrow)

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
	logical :: iffind

	iffind=.false.
	do i=Amatrow(irow),Amatrow(irow+1)-1,1
		if(Amatcol(i)==icol) then
			output=Amat(i)
			iffind=.true.
			exit
		end if
		if(Amatcol(i)>icol) exit
	end do
	if(iffind==.false.) then
		output=0.0D0
	end if

	return

end subroutine SpMatIJ

!=============================================
!=============================================
end module MathLib
