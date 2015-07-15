Module MeanField
! this subroutine is to do mean field SCF LCAO-MO calculations
! this is used in PPP model; 
! also can be used in ab initio Hamiltonian, need small change

	use variables
	use communicate
	use exit_mod
	implicit none

	public :: SCFMain

	real(kind=r8),allocatable :: & 
		oneelecH(:,:) , &         ! one electron matrix h
		twoelecG(:,:) , &         ! two electron matrix G
		fockF(:,:)    , &         ! Fock matrix
		coeffC(:,:)   , &         ! coefficient matrix
		energyE(:)    , &         ! orbital energy
		densD(:,:)    , &         ! the densD without *2, bond order matrix needs*2
		densDold(:,:) , &         ! Dij=sum(i,1->occ) C(i,miu) cross c(j,miu)*  ! read Prof.LiJun's pdf
		fockFold(:,:,:) , &       ! store the old Fock matrix  ! used in DIIS acceleration
		errvec(:,:,:)   , &       ! store the DIIS error vector
		doterrvec(:,:)            ! store the error vector inner product(matrix g in DIIS subroutine )
	integer :: nocc
	integer :: diis_subspace_size=10
	integer :: ifDIIS=.true.
	integer :: iscfiter ,&  ! at this step the scf steps index
	           ioldestfock
	contains

!=============================================================

subroutine SCFMain
! the main subrountine of SCF
	
	use mathlib
	use blas95
	use F95_precision
	implicit none
	
	integer :: guessmode,scfmaxiter
	real(kind=r8) :: norm,threshold,HFenergy,nuclrepulsion
	logical :: ifconverged
	integer :: i,j,k
	integer :: error
	real(kind=r8),allocatable :: workarray(:,:)
	
	call master_print_message("enter in SCFMain subroutine")
	
	call SCF_Allocate_Space
	allocate(workarray(norbs,norbs),stat=error)
	if(error/=0) stop
	
	if(mod(realnelecs,2)/=0) then
		call master_print_message(realnelecs,"not closed shell system")
		call SCF_Deallocate_Space
		return
	end if

	nocc=realnelecs/2   ! number of occupied orbitals
	threshold=1.0D-8   ! density matrix difference threshold
	scfmaxiter=100      ! max iterations
	ioldestfock=1

	! contruct one electron matrix
	call OneElecMat

	! get the initial guess coeff
	guessmode=1
	call SCFGuess(guessmode)

	! construct the density matrix
	call gemm(coeffC(:,1:nocc),coeffC(:,1:nocc),densD,'N','T')

	do i=1,scfmaxiter,1
		iscfiter=i
		! the last step density matrix
		densDold=densD

		! construct the two electron matrix G
		call TwoElecMat
		
		! construct the Fock matrix
		call ConstructFockMat
		
		! fockF will be changed using the diagonalization ; can not use again
		call Diagsyev(norbs,fockF,energyE,coeffC)
		
		! construct the new density matrix
		call gemm(coeffC(:,1:nocc),coeffC(:,1:nocc),densD,'N','T')
		
		! check if converged
		ifconverged=.true.
		do j=1,norbs,1
		do k=j,norbs,1  ! densD is a symmetry matrix
			norm=densDold(k,j)-densD(k,j)
			if(abs(norm)>threshold) then
				ifconverged=.false.
				exit
			end if
		end do
		end do
		
		write(*,*) "================================="
		write(*,*) "orbital energy iiter=",i
		write(*,*) energyE
		write(*,*) "================================="

		if(ifconverged==.true.) then
		!	write(*,*) densD
			call master_print_message("The SCF procedure has converged!")
			exit
		else if(i==scfmaxiter) then
			call master_print_message("not converged! The SCF procedure has reached maxiter!")
		end if
	end do

	! write the MO information
	open(unit=150,file="MO.out",status="replace")
	do i=1,norbs,1
		write(150,*) i,energyE(i)
		write(150,*) coeffC(:,i)
	end do
	! AO density matrix 
	write(150,*) densD
	close(150)

	! nuclear repulsion energy
	nuclrepulsion=0.0D0
	do i=1,norbs,1
	do j=i+1,norbs,1
		nuclrepulsion=nuclrepulsion+pppV(j,i)*nuclQ(i)*nuclQ(j)
	end do
	end do

	! Hatree Fock energy  E=sum(ij) [Dji*(2hij+Gij)]
	HFenergy=nuclrepulsion
	fockF=2.0D0*oneelecH+twoelecG                 ! fockF is just a workspace
	call gemm(densD,fockF,workarray,'N','N')     
	do i=1,norbs,1
		HFenergy=HFenergy+workarray(i,i)
	end do
	call master_print_message(nuclrepulsion,"nuclrepulsion=")
	call master_print_message(HFenergy,"HFenergy=")

	call Mean_BondOrd
	
	deallocate(workarray)

return

end subroutine SCFMain

!=============================================================
!=============================================================

subroutine ConstructFockMat
	
	implicit none
	real(kind=r8) :: diis_weight(diis_subspace_size)
	integer :: i

	if(ifDIIS==.true.) then
		call accelerate_DIIS(diis_weight)
		fockF=0.0D0
		do i=1,diis_subspace_size,1
			fockF=diis_weight(i)*fockFold(:,:,i)+fockF
		end do
	else
		fockF=oneelecH+twoelecG
	end if
return
end subroutine ConstructFockMat

!=============================================================
!=============================================================

subroutine accelerate_DIIS(diis_weight)
!  the algorithm follows Q-chem 4.3 DIIS
!  reference : Helgaker's book Moleucular electronic structure theory
!  the  error vector : e(i) = SD(i)F(i)-F(i)D(i)S
!  traces of error vectors' matrix products form a matrix
!  g= ( tr[e(1)^{\dag}e(1)] tr[e(1)^{\dag}e(2)] ... -1  )
!     ( tr[e(2)^{\dag}e(1)] tr[e(2)^{\dag}e(2)] ... -1  )
!     (   ...                  ...              ... ... )
!     (    -1                  -1               ...  0  )
!  the new fock matrix is constructed using the diis coefficient c
!   F(new) = \sum_i c(i) F(i)
!   c=g^{-1} b, where b=(0,0,...,-1)
	USE MKL95_PRECISION
	USE MKL95_LAPACK
	USE MKL95_BLAS
	implicit none
	
	real(kind=r8) :: diis_weight(diis_subspace_size)
	! local
	real(kind=r8) :: midmat(norbs,norbs)
	real(kind=r8),allocatable :: b(:,:),doterrvecdummy(:,:)
	integer,allocatable :: ipiv(:)
	real(kind=r8) :: tmp
	integer :: i,j,dim1
	integer :: info

	diis_weight=0.0D0

	!update the DIIS space fock matrix and errvector matrix
	fockFold(:,:,ioldestfock)=oneelecH+twoelecG
	call gemm(densD,fockFold(:,:,ioldestfock),midmat)
	call gemm(fockFold(:,:,ioldestfock),densD,errvec(:,:,ioldestfock))
	errvec(:,:,ioldestfock)=midmat-errvec(:,:,ioldestfock)
	! <ei|ej>
	do i=1,diis_subspace_size,1
		doterrvec(i,ioldestfock)=0.0D0
		do j=1,norbs,1
			tmp=dot(errvec(:,j,i),errvec(j,:,ioldestfock))
			doterrvec(i,ioldestfock)=doterrvec(i,ioldestfock)+tmp
		end do
		doterrvec(ioldestfock,i)=doterrvec(i,ioldestfock)
	end do

	! AX=B
	dim1=min(iscfiter,diis_subspace_size)

	allocate(doterrvecdummy(dim1+1,dim1+1))
	allocate(b(dim1+1,1))
	allocate(ipiv(dim1+1))
	doterrvecdummy=-1.0D0
	doterrvecdummy(1:dim1,1:dim1)=doterrvec(1:dim1,1:dim1)
	doterrvecdummy(dim1+1,dim1+1)=0.0D0
	b=0.0D0
	b(dim1+1,1)=-1.0D0
	call SYTRF(doterrvecdummy,'U',ipiv,info)
	if(info/=0) then
		write(*,*) "=================="
		write(*,*) "SYTRF info/=0",info
		write(*,*) "=================="
		stop
	end if
	call SYTRS2 (doterrvecdummy,b,ipiv,'U',info)
	if(info/=0) then
		write(*,*) "=================="
		write(*,*) "SYTRS info/=0",info
		write(*,*) "=================="
		stop
	end if
	diis_weight(1:dim1)=b(1:dim1,1)
	deallocate(doterrvecdummy,b,ipiv)

	! update the last Fock matrix index
	ioldestfock=ioldestfock+1
	if(ioldestfock>diis_subspace_size) then
		ioldestfock=ioldestfock-diis_subspace_size
	end if
	
return
end subroutine accelerate_DIIS
	
	
!=============================================================
!=============================================================

subroutine SCF_Allocate_Space
	implicit none
	integer :: error
	
	allocate(oneelecH(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(twoelecG(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(fockF(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(coeffC(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(energyE(norbs),stat=error)
	if(error/=0) stop
	allocate(densD(norbs,norbs),stat=error)
	if(error/=0) stop
	allocate(densDold(norbs,norbs),stat=error)
	if(error/=0) stop
	if(ifDIIS==.true.) then
		allocate(fockFold(norbs,norbs,diis_subspace_size),stat=error)
		if(error/=0) stop
		allocate(errvec(norbs,norbs,diis_subspace_size),stat=error)
		if(error/=0) stop
		allocate(doterrvec(diis_subspace_size+1,diis_subspace_size+1),stat=error)
		if(error/=0) stop
		doterrvec=0.0D0
		doterrvec(diis_subspace_size+1,:)=-1.0D0
		doterrvec(:,diis_subspace_size+1)=-1.0D0
		doterrvec(diis_subspace_size+1,diis_subspace_size+1)=0.0D0
	end if

	return
end subroutine SCF_Allocate_Space

!=============================================================
!=============================================================

subroutine SCF_Deallocate_Space
	
	implicit none
	
	deallocate(oneelecH)
	deallocate(twoelecG)
	deallocate(fockF)
	deallocate(coeffC)
	deallocate(energyE)
	deallocate(densD)
	deallocate(densDold)
	
	if(ifDIIS==.true.) then
		deallocate(fockFold)
		deallocate(errvec)
		deallocate(doterrvec)
	end if
return
end subroutine SCF_Deallocate_Space

!=============================================================
!=============================================================

subroutine OneElecMat
! construct one electron term matrix hij in PPP model
	implicit none

	integer :: icol,irow,j

	do icol=1,norbs,1
	do irow=1,norbs,1
		oneelecH(irow,icol)=t(irow,icol)
		if(irow==icol) then
			do j=1,norbs,1
				if(j/=irow) then
					oneelecH(irow,icol)=oneelecH(irow,icol)-pppV(irow,j)*nuclQ(j)
				end if
			end do
		end if
	end do
	end do

return

end subroutine OneElecMat

!=============================================================
!=============================================================

subroutine TwoElecMat
! construct the two electron G in PPP model
! twoelecG = sum(kl) [Dlk*2*(ij|kl)-(il|kj)]
	implicit none
	integer :: irow,icol,j
	
	twoelecG=0.0D0

	do icol=1,norbs,1
	do irow=1,norbs,1
		if(irow==icol) then
			do j=1,norbs,1
				! coulomb integral without hubbardU
				if(irow/=j) then
					twoelecG(irow,icol)=pppV(irow,j)*2.0D0*densD(j,j)+twoelecG(irow,icol)
				end if
			end do
			! exchange integral with hubbardU 2(ii|ii)-(ii|ii)
			twoelecG(irow,icol)=twoelecG(irow,icol)+hubbardU(irow)*densD(irow,irow)
		else
			twoelecG(irow,icol)=pppV(irow,icol)*(-1.0D0)*densD(irow,icol)
		end if
	end do
	end do

return

end subroutine TwoElecMat

!=============================================================
!=============================================================

subroutine SCFGuess(guessmode)
! the initial Guess of SCF
! guessmode=1 :: the diagonalize oneelecH scheme
	use mathlib
	implicit none
	integer :: guessmode
	
	if(guessmode==1) then
		! the diagonalization oneelecH as guess coeff
		fockF=oneelecH
		call Diagsyev(norbs,fockF,energyE,coeffC)
	end if

return

end subroutine SCFGuess

!=============================================================
!=============================================================

! nouse now
subroutine H2FCI
! this subroutine is to test if the SCF is right
! us the two orbital H2 model 4*4 singlet FCI matrix 
	implicit none

	real(kind=r8) :: fcimat(4,4)

!	fcimat(1,1)=HFenergy-nuclrepulsion
	fcimat(1,2)=0.0D0
	fcimat(1,3)=0.0D0
	fcimat(1,4)=0.0D0
return

end subroutine H2FCI

!=============================================================
!=============================================================

! nouse now
subroutine Motra
! this subroutine is to store the one electron MO integral 
! and two electron MO integral in PPP model
	implicit none
	
	real(kind=r8),allocatable :: AOOneInt(:,:),MOOneInt(:,:),&
	AOTwoInt(:,:,:,:),MOTwoInt(:,:,:,:)
	integer :: reclength

	reclength=2
	
!	open(unit=151,file="AOOneInt.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")
!	open(unit=152,file="MOOneInt.tmp",access="Direct",form="unformatted",recl=reclength,status="replace")

return

end subroutine Motra

!=============================================================
!=============================================================

subroutine Mean_BondOrd
! this subroutine calculate the mean field bond order matrix
	
	use blas95
	use f95_precision
	implicit none

	real(kind=r8) :: mean_bomat(norbs,norbs)
	integer :: i,j

	open(unit=1002,file="mean_bomat.out",status="replace")
	do i=1,norbs,1
	do j=i,norbs,1
		if(bondlink(i,j)/=0) then
			mean_bomat(i,j)=dot(coeffC(i,1:nocc),coeffC(j,1:nocc))
			mean_bomat(i,j)=mean_bomat(i,j)*2.0D0    ! up down spin
			write(1002,*) i,j,mean_bomat(i,j)
		end if
	end do
	end do
	close(1002)
return
end subroutine Mean_BondOrd

!=============================================================
!=============================================================
end module MeanField

