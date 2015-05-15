Module MeanField
! this subroutine is to do mean field SCF LCAO-MO calculations
! this is used in PPP model; 
! also can be used in ab inition Hamiltonian, need small change

	use variables
	use communicate
	use exit_mod
	implicit none
	private
	save

	public :: SCFMain

	real(kind=r8),allocatable :: & 
		oneelecH(:,:) , &         ! one electron matrix h
		twoelecG(:,:) , &         ! two electron matrix G
		fockF(:,:)    , &         ! Fock matrix
		coeffC(:,:)   , &         ! coefficient matrix
		energyE(:)    , &         ! orbital energy
		densD(:,:)    , &         ! the densD without *2, bond order matrix needs*2
		densDold(:,:) , &         ! Dij=sum(i,1->occ) C(i,miu) cross c(j,miu)*  ! read LiJun's pdf
		workarray(:,:)            
	contains

!=============================================================

subroutine SCFMain
! the main subrountine of SCF
	
	use mathlib
	use blas95
	use F95_precision
	implicit none
	
	integer :: guessmode,scfmaxiter,nocc
	real(kind=r8) :: norm,threshold,HFenergy,nuclrepulsion
	logical :: ifconverged
	integer :: i,j,k
	integer :: error
	
	call master_print_message("enter in SCFMain subroutine")

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
	allocate(workarray(norbs,norbs),stat=error)
	if(error/=0) stop
	
	if(mod(realnelecs,2)/=0) then
		call master_print_message(realnelecs,"not closed shell system")
		stop
	end if
	nocc=realnelecs/2   ! number of occupied orbitals
	threshold=1.0D-10   ! density matrix difference threshold
	scfmaxiter=100      ! max iterations

	! contruct one electron matrix
	call OneElecMat

	! get the initial guess coeff
	guessmode=1
	call SCFGuess(guessmode)

	! construct the density matrix
	call gemm(coeffC(:,1:nocc),coeffC(:,1:nocc),densD,'N','T')

	do i=1,scfmaxiter,1
		! the last step density matrix
		densDold=densD

		! construct the two electron matrix G
		call TwoElecMat
		
		! construct the Fock matrix
		fockF=oneelecH+twoelecG
		
		! fockF will be changed using the diagonalization ; can not use again
		call Diagsyev(norbs,fockF,energyE,coeffC)
		
		! construct the new density matrix
		call gemm(coeffC(:,1:nocc),coeffC(:,1:nocc),densD,'N','T')
		
		! check if converged
		ifconverged=.true.
		do j=1,norbs,1
		do k=1,norbs,1
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
			call master_print_message("The SCF procedure has converged!")
			exit
		else if(i==scfmaxiter) then
			call master_print_message("not converged! The SCF procedure has reached maxiter!")
		end if
	end do

	! write the MO information
	open(unit=150,file="MO.tmp",status="replace")
	do i=1,norbs,1
		write(150,*) i,energyE(i)
		write(150,*) coeffC(:,i)
	end do
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

	
	deallocate(oneelecH)
	deallocate(twoelecG)
	deallocate(fockF)
	deallocate(coeffC)
	deallocate(energyE)
	deallocate(densD)
	deallocate(densDold)
	deallocate(workarray)

return

end subroutine SCFMain

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
end module MeanField

