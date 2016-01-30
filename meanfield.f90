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
    
    ! MO integral input
    integer,allocatable :: nactmoa(:)
    integer ::  nmoa
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
            if(ifconverged==.false.) exit
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
    
!   call Motra
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
!   USE F95_PRECISION
!   USE LAPACK95
!   USE BLAS95
    USE F95_PRECISION
    USE LAPACK95
    USE BLAS95
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
        fockFold=0.0D0
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

!   fcimat(1,1)=HFenergy-nuclrepulsion
    fcimat(1,2)=0.0D0
    fcimat(1,3)=0.0D0
    fcimat(1,4)=0.0D0
return

end subroutine H2FCI

!=============================================================
!=============================================================

subroutine Motra
! this subroutine is to store the one electron MO integral 
! and two electron MO integral in PPP model
    implicit none
    
    ! the last letter a indicate the active space
    real(kind=r8),allocatable :: ea(:), xa(:)
    integer(kind=i4) :: nocca,nvira,nblkocca,nblkvira,nblocka, nmo2a, nxa
    integer :: i,j,k,l
    integer :: nstarta,nenda,ii,jj,kk,ll,ij,kl,ijkl
    integer :: ierr
    real(kind=r8),parameter :: autoeV=27.211D0

    
    ! (ij|kl)
    open(unit=56,file='active.inp',status="old")
    read(56,*) nocca,nvira
    nmoa=nocca+nvira
    nmo2a=(nmoa+1)*nmoa/2  ! (ij| or |kl) pair number
    nxa=(nmo2a+1)*nmo2a/2  ! (ij|kl) pair number

    allocate(ea(nmoa), xa(nxa), nactmoa(nmoa), stat=ierr)
    if(ierr/=0) stop

    read(56,*) nblkocca,nblkvira
    nblocka=nblkocca+nblkvira
    
    k=0
    do i=1,nblocka,1
        read(56,*) nstarta,nenda
        do j=nstarta,nenda,1
            k=k+1
            nactmoa(k)=j
        end do
    end do
    close(56)
      
    ! orbital energy
    do i=1,nmoa,1
        ea(i)=energyE(nactmoa(i))/autoeV
    end do
    
    ! two electron integral
    do i=1,nmoa,1
    do j=1,i,1
        ij=(i-1)*i/2+j
        do k=1,nmoa,1
        do l=1,k,1
            kl=(k-1)*k/2+l
            if (ij>=kl) then
                ijkl=(ij-1)*ij/2+kl
                ii=nactmoa(i)
                jj=nactmoa(j)
                kk=nactmoa(k)
                ll=nactmoa(l)
                call IntAOtoMO(xa(ijkl),ii,jj,kk,ll)
                xa(ijkl)=xa(ijkl)/autoeV
            end if
        end do
        end do
    end do
    end do

    open(unit=55,file='moint2.out',form='unformatted',status="replace")
    write(55) nmoa, nxa
    write(*,*) "nmoa=",nmoa,"nxa",nxa
    write(55) (ea(i), i=1,nmoa)
    write(55) (xa(i), i=1,nxa)
    close(55)

    ! calculate the transition dipole
    call transdipol
    deallocate(ea, xa, nactmoa)
return

end subroutine Motra

!=============================================================
!=============================================================

subroutine IntAOtoMO(integral,ii,jj,kk,ll)
! this subroutine do PPP model AO to MO transformation on the fly
    implicit none

    real(kind=r8) :: integral
    integer :: ii,jj,kk,ll
    integer :: l,r
    
    integral=0.0D0
    do l=1,norbs,1
    do r=1,norbs,1
        if(l==r) then
            integral=integral+hubbardU(l)*coeffC(l,ii)*coeffC(l,jj)*coeffC(r,kk)*coeffC(r,ll)
        else
            integral=integral+pppV(l,r)*coeffC(l,ii)*coeffC(l,jj)*coeffC(r,kk)*coeffC(r,ll)
        end if
    end do
    end do
return
end subroutine IntAOtoMO

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

subroutine transdipol
! this subroutine calculate the transition dipole moment between 
! different MO <MO1|r|MO2>, the reference point is the center of mass
! and calculate the HF reference dipole moment <HF|r|HF> and
! the nuclear dipole moment
! in the PPP model only the ni operator contribute the atomic orbital dipole
! moment

    implicit none
    real(kind=r8),allocatable :: trnsdipmo(:,:,:)
    real(kind=r8) :: trnshf(3),trnsnuc(3)
    integer :: i,j,k,l

    ! <MO1|r|MO2>
    allocate(trnsdipmo(norbs,norbs,3))
    trnsdipmo=0.0D0
    ! in the e*angstrom unit
    do l=1,3,1
    do i=1,norbs,1
    do j=1,i,1
        do k=1,norbs,1
            trnsdipmo(j,i,l)=trnsdipmo(j,i,l)-coeffC(k,i)*coeffC(k,j)*(coord(l,k)-cntofmass(l))
            ! the negative sign here represents the negative electron charge
        end do
        trnsdipmo(i,j,l)=trnsdipmo(j,i,l)
    end do
    end do
    end do

    ! HF dipole moment
    trnshf=0.0D0
    do l=1,3,1
    do i=1,norbs,1
        trnshf(l)=trnshf(l)+trnsdipmo(i,i,l)
    end do
    end do
    trnshf=trnshf*2.0D0

    ! nuclear dipole moment
    trnsnuc=0.0D0
    do l=1,3,1
    do i=1,natoms,1
        trnsnuc(l)=trnsnuc(l)+(coord(l,i)-cntofmass(l))*nuclQ(i)
    end do
    end do

    ! MRDCI output
    open(unit=60,file="ed.mo.out",form="unformatted",status="replace")
    do l=1,3,1
        write(60) (trnshf(l)+trnsnuc(l))*eAtoDebye,nmoa,l
        write(60) ((trnsdipmo(nactmoa(i),nactmoa(j),l)*eAtoDebye,i=1,nmoa),j=1,nmoa)
    end do
    close(60)
    
    ! EOM-CCSD output
    open(unit=61,file="dip.mo.out",form="unformatted",status="replace")
    do l=1,3,1
        write(61) trnshf(l)/AutoAngstrom,trnsnuc(l)/AutoAngstrom
        write(61) ((trnsdipmo(nactmoa(i),nactmoa(j),l)/AutoAngstrom,i=1,j),j=1,nmoa)
    end do
    close(61)

    deallocate(trnsdipmo)
return

end subroutine transdipol

!=============================================================
!=============================================================
end module MeanField

