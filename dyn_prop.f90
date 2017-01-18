Module dyn_prop 
!> this module calculate the dynamic property using the 
!! Lanczos fraction representation method + planso solver
!! correction vector method + mkl fgmres solver 
!! dynamic DMRG method + mkl cg solver
!! the cg solver seems more robust than fgmres solver

!! the iterative version is not implemented

    use variables
    use mathlib
    use hamiltonian_mod
    use communicate
    USE MPI
    implicit none
    
    logical :: Ifdyn_prop     ! flag to determine if do dynamic property calculation
    integer(kind=i4) :: maxnlancs, &      ! max # of lanczos vectors
                        dyn_initstate     ! dynamic property initial staet  
    character(len=20) :: dyn_prop_method  ! calculation method
    real(kind=r8) :: normthresh = 1.0D-7 , & ! threshold to neglect the initial state  R|0> ?= 0
                     eta = 1.0D-1            ! Lorentzian broaden parameter
    
contains
!==========================================================
!==========================================================

subroutine dyn_prop_wrapper
    implicit none

    call master_print_message("Enter dynamic property calculation subroutine")
    if(dyn_prop_method == "Lanczos") then
        call master_print_message("Lanczos fraction representation method")
        call dyn_Lanczos_wrapper
    else if(dyn_prop_method == "CV") then
        call master_print_message("correction vector method")
        call dyn_CV_wrapper
    else if(dyn_prop_method == "DDMRG") then
        call master_print_message("dynamic DMRG method")
        call dyn_dDMRG_wrapper
    else
        write(*,*) "dyn_prop method wrong!", dyn_prop_method
    end if
    
    return

end subroutine dyn_prop_wrapper

!==========================================================
!==========================================================

subroutine dyn_Lanczos_wrapper
!> the Lanczos subroutine is 
!! LANSO -- simple LANczos with Selective Orthogonalization (version 1.0)
!! by Professor Beresford Parlett of UC Berkeley
    
    implicit none
    
    real(kind=r8) :: ENDL, ENDR, KAPPA, CONDM, norm(3)
    integer :: LANMAX, worksize, ierr, nlanczos, nconverged, imessage, &
               ijob,ierror,EV,icoord,i
    real(kind=r8),allocatable :: workarray(:), eigenvalue(:), errbnd(:), &
        dummycoeff(:),dummynewcoeff(:), eigenvector(:), dyninitvec(:,:)

    call GetDimSym
    
    LANMAX = MIN(maxnlancs,dimN)
    worksize = 6*dimN+1+4*LANMAX+LANMAX*LANMAX
    
    ! allocate workarray
    if(myid==0) then
        allocate(workarray(worksize))
        allocate(eigenvalue(LANMAX))
        allocate(errbnd(LANMAX))
        allocate(eigenvector(1*dimN))
        allocate(dyninitvec(dimN,3))
    end if

    ! construct the complementary operator matrix
    if(opmethod=="comple") then
        call Complement(operamatbig1,bigcolindex1,bigrowindex1,Lrealdim,Rrealdim,subM,0)
    end if
    
    ! get the initial lanczos vector
    call dyn_init(dimN, dyninitvec, norm)
    
    if(myid==0) then 
        write(*,*) "initial vector norm= ", norm
        ! normalize the initial vector
        do icoord = 1, 3, 1
            dyninitvec(:,icoord) = dyninitvec(:,icoord) / sqrt(norm(icoord))
        end do
    end if

    do icoord = 1, 3, 1 
        ! if the initial vector is too small, neglect it
        if(norm(icoord) < normthresh) cycle
        
        if(myid==0) then
            workarray(1:dimN) = dyninitvec(:,icoord)
            
            ! parameter
            ENDL = -1.0D-30
            ENDR =  1.0D-30
            imessage = 3
            CONDM = 1.0D0
            EV = 99         ! output channel
            KAPPA = 1.0D-6  
            ! KAPPA       ... {DOUBLE PRECISION} relative accuracy of Ritz values
            ! acceptable as eigenvalues; specifically, LANDR
            ! will deliver RITZ(K) along with its error bound
            ! BND(K) and optionally its associated Ritz vector
            ! if BND(K).LE.KAPPA*ABS(RITZ(K)), (KAPPA will be
            ! reset to MAX(KAPPA,macheps**(3/4)) if EV.GT.0,
            ! otherwise KAPPA is untouched)

            
            !CALL LANDR(N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
            !    J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
            call LANDR(dimN, LANMAX, 1, CONDM, ENDL, ENDR, EV, KAPPA, &
                nlanczos, nconverged, eigenvalue, errbnd, workarray, worksize, ierror, imessage, eigenvector)
            
            ! the normal exit signal (IJOB = 0)
            IJOB = 0
            call MPI_BCAST(IJOB,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
            
        else
            do while (.true.)
                call MPI_BCAST(IJOB,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
                if(IJOB == 1) then
                    allocate(dummycoeff(1))
                    allocate(dummynewcoeff(1))
                
                    if(opmethod=="comple") then
                        call op(1,1,dummycoeff,dummynewcoeff,&
                                Lrealdim,Rrealdim,subM,ngoodstates,&
                                operamatbig1,bigcolindex1,bigrowindex1,&
                                Hbig,Hbigcolindex,Hbigrowindex,&
                                quantabigL,quantabigR,goodbasis,goodbasiscol)
                    
                    else if(opmethod=="direct") then
                        call opdirect(1,1,dummycoeff,dummynewcoeff,&
                                Lrealdim,Rrealdim,subM,ngoodstates,&
                                operamatbig1,bigcolindex1,bigrowindex1,&
                                Hbig,Hbigcolindex,Hbigrowindex,&
                                quantabigL,quantabigR,goodbasis,goodbasiscol)
                    end if

                    deallocate(dummycoeff)
                    deallocate(dummynewcoeff)
                else
                    exit
                end if
            end do
        end if

        ! the output 
        if(myid == 0) then
            write(*,*) "nlanczos = ", nlanczos
            write(*,*) "nconverged = ", nconverged
            write(*,*) "index, eigenvalue, errbound"
            do i =1, nlanczos, 1
                write(*,*) i, eigenvalue(i), errbnd(i)
            end do
        end if
    end do
    
    if(opmethod=="comple") then
        call deallocate_Complement
    end if

    if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
        call DestroySymm
    end if

    if(myid==0) then
        deallocate(workarray)
        deallocate(eigenvalue)
        deallocate(errbnd)
        deallocate(eigenvector)
        deallocate(dyninitvec)
    end if

    return

end subroutine dyn_Lanczos_wrapper

!==========================================================
!==========================================================

subroutine dyn_init(vecdim, dyninitvec, norm)
!> master process calculate the initial wfn vector
!! in the absorption spectrum 
!! \hat{r}|psi_{gs}\rangle
!! in the emission spectrum 
!! \hat{r}|psi_{ex}\rangle

    use module_sparse
    use symmetry
    use opc
    implicit none
    integer,intent(in) :: vecdim
    real(kind=r8),intent(out) ::  dyninitvec(:,:), norm(:)
    !real(kind=r8),intent(out) ::  dyninitvec(vecdim,3), norm(3)

    !  local  
    real(kind=r8),allocatable :: dynvec(:,:), dynvec0(:,:,:), dynvecbuf(:,:), dyn1dvec0(:,:)
    integer :: iorb, icoord, iproc, operaindex, ierr, mbasis, ibasis
    character(len=1) :: domain
    ! MPI flag
    integer :: status(MPI_STATUS_SIZE)

    
    if(myid == 0) then
        allocate(dynvec0(4*Lrealdim, 4*Rrealdim, 3))
        dynvec0 = 0.0D0
        allocate(dynvecbuf(4*Lrealdim, 4*Rrealdim))
        allocate(dyn1dvec0(ngoodstates,3))
    else
        allocate(dynvec(4*Lrealdim, 4*Rrealdim))
    end if

    !> in the PPP/hubbard model, the r is simplifiled to n_i
    do iorb = 1, norbs, 1
        if(myid /= 0) then
            iproc = orbid1(iorb,1)
            if(myid == iproc) then
                operaindex = orbid1(iorb,2)*3
                if(iorb <= nleft+1) then
                    domain = "L"
                else
                    domain = "R"
                end if
                
                ! n_iorb * C 
                call SubSpaceOpCdens(dynvec, domain,&
                        Lrealdim,Rrealdim,dyn_initstate,operaindex,&
                        operamatbig1,bigcolindex1,bigrowindex1,&
                        coeffIF,coeffIfcolindex,coeffIFrowindex)

                call MPI_SEND(dynvec,16*Lrealdim*Rrealdim,MPI_real8,0,iorb,MPI_COMM_WORLD,ierr)
            end if
        else
            call MPI_RECV(dynvecbuf,16*Lrealdim*Rrealdim,MPI_real8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
            !   coord(iorb,1:3)*n_iorb*c
            do icoord = 1, 3, 1
                call mataxpy(coord(icoord,status(MPI_TAG)), 4*Lrealdim, 4*Rrealdim, dynvecbuf, dynvec0(:,:, icoord))
            end do
        end if
    end do
    
    if(myid == 0) then

        do icoord = 1, 3, 1
            ! copy to 1d array
            mbasis = 0
            do ibasis = 1, ngoodstates, 1
                mbasis = mbasis +1
                dyn1dvec0(mbasis, icoord) = dynvec0(goodbasis(ibasis,1), goodbasis(ibasis,2), icoord)
            end do

            ! symmetrize
            if(logic_spinreversal/=0) then
                call symmetrizestate(ngoodstates, dyn1dvec0(:,icoord), dyninitvec(:,icoord), 's')
            else
                call copy(dyn1dvec0(:,icoord), dyninitvec(:,icoord))
            end if
            norm(icoord) = dot(dyninitvec(:,icoord), dyninitvec(:,icoord))
            !dyninitvec(:,icoord) = dyninitvec(:,icoord) / sqrt(norm(icoord))
        end do
    end if

    call MPI_BCAST(norm,3,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    ! deallocate memory
    if(myid == 0) then
        deallocate(dynvec0, dynvecbuf, dyn1dvec0)
    else
        deallocate(dynvec)
    end if

    return
end subroutine dyn_init

!==========================================================
!==========================================================

subroutine dyn_CV_wrapper
!> use the correction vector method to calculate the dynamic correlation function
!! use the Intel MKL RCI (P)FGMRES ((Preconditioned) Flexible Generalized Minimal
!! RESidual method) to solve this non-symmetry complex Ax = b equation
!! we transform the complex equation to a real equation but the dim_real = 2* dim_cmplx
!! please see more details in the manual

!> the problem of this CV method is that it can not converge very well
!! especially in the large omega part

!! correction vector method can get a good initial guess like the original DMRG eigen problem
    implicit none
    
    integer :: cvdim
    real(kind=r8),allocatable :: cv(:), rhs(:,:), workarray(:), rhsold(:,:), tmparray(:) ,Hdiag(:)
    integer :: IPAR(128), RCI_REQUEST, itercount, IJOB, maxspace, icoord, ierr, id, i
    real(kind=r8) :: DPAR(128),  norm(3)

    ! omega part
    real(kind=r8) :: omega,omegal,omegar
    integer :: ipoint,npoints
    real(kind=r8),allocatable :: dyn_cor(:,:)

    call GetDimSym
    
    ! transform the complex equation to real equation
    cvdim = 2 * dimN
    
    ! define the omega part
    omegal = 0.0
    omegar = 40.0
    npoints = 400
    if(myid == 0) then
        allocate(dyn_cor(npoints,4))
        dyn_cor = 0.0D0
        do ipoint = 1, npoints, 1
            dyn_cor(ipoint,4) = (omegar-omegal)/npoints*ipoint+omegal
        end do
    end if
    
    !> maximun space to restart calculation
    maxspace = 600
    maxspace = MIN(maxspace, cvdim)

    if(myid == 0) then
        allocate(cv(cvdim))
        cv = 0.0D0
        allocate(rhs(cvdim,3))
        allocate(rhsold(dimN,3))
        allocate(workarray((2*maxspace+1)*cvdim+maxspace*(maxspace+9)/2+1))
        allocate(tmparray(cvdim))
        allocate(HDIAG(dimN))
    else
        allocate(tmparray(2))
        allocate(workarray(2))
        IPAR = 1
        dimN = 1
    end if
    
    ! construct the complementary operator matrix
    if(opmethod=="comple") then
        call Complement(operamatbig1,bigcolindex1,bigrowindex1,Lrealdim,Rrealdim,subM,0)
    end if
    
    ! use the diagonal part to do precondition
    if(logic_spinreversal/=0 .or. &
        (logic_C2/=0 .and. nleft==nright)) then
        call SymmHDiag(HDIAG)
    else 
        call GetHDiag(HDIAG,dimN,&
        operamatbig1,bigcolindex1,bigrowindex1,&
        Hbig,Hbigcolindex,Hbigrowindex,&
        goodbasis,&
        .false.)
    end if

    !> get the rhs of Ax=b
    call dyn_init(dimN, rhsold, norm)
    
    if(myid == 0) then
        open(unit = 222, file="CV.out",status="replace")
    end  if
    
    ! first the icoord then the ipoint
    ! suppose we can use the cv last step as a initial guess
    do icoord = 1,3,1 
    do ipoint =1, npoints, 1 
        if(norm(icoord) < normthresh) cycle
        
        if(myid == 0) then
            ! it seems the zero initial guess can converge fastly.
            ! but it affects a little
            cv = 0.0D0
            !> set the omega and rhs
            omega = dyn_cor(ipoint,4)
            rhs(:, icoord)  = 0.0D0
            rhs(1:dimN,icoord) = rhsold(:,icoord)
            
            DO I = 1, 128, 1
                IPAR(i) = 0
                DPAR(i) = 0.0D0
            END DO
            
            !> initialize the FGMRESsolver
            CALL DFGMRES_INIT(cvdim, cv, rhs(:,icoord), RCI_REQUEST, IPAR, DPAR, workarray)
            if (RCI_REQUEST /= 0) then
                write(*,*) "RCI_REQUEST /= 0, INIT STOP"
                CALL MKL_FREE_BUFFERS
                stop
            end if

            ! set paramters for no-preconditioner calculation
            ! use the default test
            IPAR(9) = 1
            IPAR(10) = 0
            IPAR(11) = 1  ! the diagonal precondition method is generally better than non-preconditioned version
            IPAR(12) = 1
            IPAR(5) = maxspace
            IPAR(15) = maxspace
            
            !> check the parameters
            CALL DFGMRES_CHECK(cvdim, cv, rhs(:,icoord), RCI_REQUEST, IPAR, DPAR, workarray)
            if (RCI_REQUEST /= 0) then
                write(*,*) "RCI_REQUEST /= 0, CHECK STOP"
                CALL MKL_FREE_BUFFERS
                stop
            end if
        end if
        

        do while(.true.)
            if(myid ==0) then
                !> do FGMRES
                CALL DFGMRES(cvdim, cv, rhs(:,icoord), RCI_REQUEST, IPAR, DPAR, workarray)
                if(RCI_REQUEST == 1) then
                    IJOB = 1
                else if (RCI_REQUEST == 0) then
                    IJOB = 0
                else if (RCI_REQUEST == 3) then
                    do id = 1, dimN, 1
                        workarray(IPAR(23)+id-1) = workarray(IPAR(22)+id-1) / (-HDIAG(id)+omega+dmrgenergy(dyn_initstate))
                        workarray(IPAR(23)+id-1+dimN) = workarray(IPAR(22)+id-1+dimN) / (-HDIAG(id)+omega+dmrgenergy(dyn_initstate))
                    end do
                    IJOB = 3
                else
                    write(*,*) "RCI_REQUEST /= 0,/=1, DFGMRES STOP", RCI_REQUEST
                    IJOB = -1
                end if
            end if

            call MPI_BCAST(IJOB,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
            
            if(IJOB==1) then
                !> calculate HC
                if(opmethod=="comple") then
                    call op(dimN,2,workarray(IPAR(22)),tmparray,&
                        Lrealdim,Rrealdim,subM,ngoodstates,&
                        operamatbig1,bigcolindex1,bigrowindex1,&
                        Hbig,Hbigcolindex,Hbigrowindex,&
                        quantabigL,quantabigR,goodbasis,goodbasiscol)
                else if(opmethod=="direct") then
                    call opdirect(dimN,2,workarray(IPAR(22)),tmparray,&
                        Lrealdim,Rrealdim,subM,ngoodstates,&
                        operamatbig1,bigcolindex1,bigrowindex1,&
                        Hbig,Hbigcolindex,Hbigrowindex,&
                        quantabigL,quantabigR,goodbasis,goodbasiscol)
                else
                    stop
                end if
                
                ! A*x+\Gamma*y
                if(myid == 0) then
                    workarray(IPAR(23):IPAR(23)+dimN-1) = -1.0D0*tmparray(1:dimN) + &
                        (omega+dmrgenergy(dyn_initstate))*workarray(IPAR(22):IPAR(22)+dimN-1) - &
                        eta*workarray(IPAR(22)+dimN:IPAR(22)+dimN*2-1)
                    workarray(IPAR(23)+dimN:IPAR(23)+2*dimN-1) = -1.0D0*tmparray(dimN+1:2*dimN) + &
                        (omega+dmrgenergy(dyn_initstate))*workarray(IPAR(22)+dimN:IPAR(22)+2*dimN-1) + &
                        eta*workarray(IPAR(22):IPAR(22)+dimN-1)
                end if

            else if(IJOB == 3) then
                !call master_print_message(IJOB,"preconditioning")
                cycle
            else if(IJOB == 0) then
                call master_print_message(IJOB,"normal exit")
                exit
            else 
                call master_print_message(IJOB,"fatal exit")
                exit
            end if
        end do

        if(myid == 0) then
            CALL DFGMRES_GET(cvdim, cv, rhs(:,icoord), RCI_REQUEST, IPAR, DPAR, workarray, itercount)
            write(*,*) "The system has been solved"
            call master_print_message(itercount, "Number of iterations:")
            CALL MKL_FREE_BUFFERS
            dyn_cor(ipoint,icoord) = dot(cv(dimN+1:cvdim),rhsold(:,icoord))/PI*(-1.0D0)
            write(222,'(2E25.10)') dyn_cor(ipoint,icoord),dyn_cor(ipoint,4)
        end if
    end do
        if(myid==0) then
            write(222,*)
            write(222,*)
        end if
    end do         
    
    if(opmethod=="comple") then
        call deallocate_Complement
    end if

    if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
        call DestroySymm
    end if
    
    if(myid == 0) then
        close(222)
        deallocate(cv)
        deallocate(rhs)
        deallocate(rhsold)
        deallocate(dyn_cor)
        deallocate(HDIAG)
    end if
    deallocate(tmparray)
    deallocate(workarray)

    return
end subroutine dyn_CV_wrapper

!==========================================================
!==========================================================

subroutine dyn_dDMRG_wrapper
!> use the dynamic DMRG method to calculate the correction vector
!> the opt method is based on mkl RCI CG method

    implicit none

    integer,parameter :: nRhs = 3   ! # of right hand side vectors
    real(kind=r8),allocatable :: cv(:,:), rhs(:,:), rhsold(:,:), &
        workarray(:,:), tmparray(:), tmparray2(:)
    integer :: method, IPAR(128 + 2*nRhs), RCI_REQUEST, itercount(3), &
        IJOB, maxspace, icoord, ierr, istep, i, niter
    real(kind=r8) :: DPAR(128 + 2*nRhs), norm(3)

    ! omega part
    real(kind=r8) :: omega,omegal,omegar,thresh, cornew(3)
    integer :: ipoint,npoints
    real(kind=r8),allocatable :: dyn_cor(:,:)
    logical :: converged

    method = 1
    ! relative threshold to determin the correlation function is converged or not?
    thresh = 1.0D-6

    call GetDimSym
    
    omegal = 20.0
    omegar = 40.0
    npoints = 200
    if(myid == 0) then
        allocate(dyn_cor(npoints,4))
        dyn_cor = 0.0D0
        do ipoint = 1, npoints, 1
            dyn_cor(ipoint,4) = (omegar-omegal)/npoints*ipoint+omegal
        end do
    end if
    
    !> maximum space 
    maxspace = 600
    maxspace = MIN(maxspace, dimN)

    if(myid == 0) then
        allocate(cv(dimN,3))
        cv = 0.0D0
        allocate(rhs(dimN,3))
        allocate(rhsold(dimN,3))
        allocate(workarray(dimN,3+nRhs))
        !allocate(workarray(dimN,3))   ! non-preconditioned CG iterations
        allocate(tmparray(dimN))
        allocate(tmparray2(dimN))
    else
        allocate(tmparray(1))
        allocate(tmparray2(1))
        IPAR = 1
        dimN = 1
    end if
    
    ! construct the complementary operator matrix
    if(opmethod=="comple") then
        call Complement(operamatbig1,bigcolindex1,bigrowindex1,Lrealdim,Rrealdim,subM,0)
    end if
    
    !> get the rhs of Ax=b
    call dyn_init(dimN, rhsold, norm)
    if(myid == 0) then
        rhs = rhsold * (-1.0D0) * eta
    end if

    if(myid == 0) then
        open(unit = 222, file="DDMRGCV.out",status="replace")
    end  if
    
    do ipoint =1, npoints, 1
        niter = 0
        if(myid == 0) then
            !> set the omega and rhs
            omega = dyn_cor(ipoint,4)
            
            ! initalize the ipar dpar
            ! it is very essential to get correct result
            ! because the init and check only check the first a few parameters documented in the mkl manual
            DO I = 1, 128 + 2*nRhs
                ipar(i) = 0
                dpar(i) = 0.0D0
            END DO

            !> initialize the FGMRESsolver
            CALL dcgmrhs_init(dimN,cv,nRhs,rhs,method,RCI_request,ipar,dpar,workarray)
            if (RCI_REQUEST /= 0) then
                write(*,*) "RCI_REQUEST /= 0, INIT STOP"
                CALL MKL_FREE_BUFFERS
                stop
            end if
            
            !> set paramters for no-preconditioner calculation
            ipar(5)  = maxspace
            ! residual stopping test
            !ipar(9)  = 1
            !ipar(10) = 0
            
            ! user defined test
            ipar(10) = 1

            !> check the parameters
            CALL dcgmrhs_check(dimN, cv, nRhs, rhs, RCI_REQUEST, IPAR, DPAR, workarray)
            if (RCI_REQUEST /= 0) then
                write(*,*) "RCI_REQUEST /= 0, CHECK STOP"
                CALL MKL_FREE_BUFFERS
                stop
            end if
        end if
        

        do while(.true.)
            if(myid == 0) then
                !> do dcgmrhs
                CALL dcgmrhs(dimN,cv,nRhs,rhs,RCI_request,ipar,dpar,workarray)
                if(RCI_REQUEST == 1) then
                    IJOB = 1
                else if (RCI_REQUEST == 0) then
                    IJOB = 0
                else if (RCI_REQUEST == 2) then
                    IJOB = 2
                else
                    write(*,*) "RCI_REQUEST /= 0,/=1, DFGMRES STOP", RCI_REQUEST
                    IJOB = -1
                end if
            end if

            call MPI_BCAST(IJOB,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
            
            if(IJOB == 1) then
                ! there is an error here
                ! in fact H_eff*H_eff*C /= (H^2)_eff*C
                ! but we use the lhs to approximate the rhs
                niter = niter + 1

                do istep = 1,2,1
                    if(myid ==0 ) then
                        if(istep == 1) then
                            tmparray = workarray(1:dimN,1)
                        else if(istep == 2) then
                            tmparray2 = -1.0D0 * tmparray2  !-H*C
                            CALL AXPY(tmparray, tmparray2, omega+dmrgenergy(dyn_initstate)) ! (-H*C)+(E_0+\omega)*C
                            tmparray = tmparray2
                        end if
                    end if

                    !> calculate HC
                    if(opmethod=="comple") then
                        call op(dimN,1,tmparray,tmparray2,&
                            Lrealdim,Rrealdim,subM,ngoodstates,&
                            operamatbig1,bigcolindex1,bigrowindex1,&
                            Hbig,Hbigcolindex,Hbigrowindex,&
                            quantabigL,quantabigR,goodbasis,goodbasiscol)
                    else if(opmethod=="direct") then
                        call opdirect(dimN,1,tmparray,tmparray2,&
                            Lrealdim,Rrealdim,subM,ngoodstates,&
                            operamatbig1,bigcolindex1,bigrowindex1,&
                            Hbig,Hbigcolindex,Hbigrowindex,&
                            quantabigL,quantabigR,goodbasis,goodbasiscol)
                    else
                        stop
                    end if
                end do
                
                ! update (A^T*A+eta^2)\psi
                if(myid == 0) then
                    workarray(1:dimN,2) = -1.0D0 * tmparray2  !-H*C'
                    CALL AXPY(tmparray, workarray(1:dimN,2), omega+dmrgenergy(dyn_initstate)) ! -H*C'+(E_0+\omega)*C'
                    CALL AXPY(workarray(1:dimN,1), workarray(1:dimN,2), eta**2.0)  ! 
                end if

            else if(IJOB == 2) then
                if(myid == 0) then
                    do icoord = 1,3,1
                        cornew(icoord) = dot(cv(1:dimN,icoord),rhsold(:,icoord))/PI*(-1.0D0)
                    end do
                    converged = .true.
                    do icoord = 1,3,1
                        if(ABS(cornew(icoord)-dyn_cor(ipoint,icoord))/ABS(dyn_cor(ipoint,icoord)) > thresh) then
                            converged = .false.
                            exit
                        end if
                    end do
                    dyn_cor(ipoint,1:3) = cornew
                end if
                call MPI_BCAST(converged,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
                if(converged == .true.) then
                    call master_print_message("user defined test fullfilled")
                    exit
                end if
            else if(IJOB == 0) then
                call master_print_message(IJOB,"normal exit")
                exit
            else 
                call master_print_message(IJOB,"IJOB in DDMRG is wrong, fatal exit")
                exit
            end if
        end do

        if(myid == 0) then
            CALL dcgmrhs_get(dimN,cv,nRhs,rhs,RCI_request,ipar,dpar,workarray,itercount)
            write(*,*) "The system has been solved"
            write(*,*) "Number of iterations:", itercount
            write(*,*) "Number of iterations2:", niter
            CALL MKL_FREE_BUFFERS
            do icoord = 1,3,1
                dyn_cor(ipoint,icoord) = dot(cv(1:dimN,icoord),rhsold(:,icoord))/PI*(-1.0D0)
            end do
            write(222,'(4E25.10)') dyn_cor(ipoint,1:4)
        end if
    end do
    
    if(opmethod=="comple") then
        call deallocate_Complement
    end if

    if(logic_spinreversal/=0 .or. (logic_C2/=0 .and. nleft==nright)) then
        call DestroySymm
    end if
    
    if(myid == 0) then
        close(222)
        deallocate(cv)
        deallocate(rhs)
        deallocate(dyn_cor)
        deallocate(workarray)
    end if
    deallocate(tmparray)
    deallocate(tmparray2)
    
    return

end subroutine dyn_dDMRG_wrapper
!==========================================================
!==========================================================
end Module dyn_prop 
