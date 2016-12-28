Module dyn_prop 
! this module calculate the dyanmic property using the Lanczos fraction representation method
    use variables
    use mathlib
    USE MPI
    implicit none
    
    logical :: Ifdyn_prop
    integer(kind=i4) :: maxnlancs, dyn_initstate 

contains
!==========================================================
!==========================================================

subroutine dyn_prop_wrapper
    use hamiltonian_mod
    implicit none
    
    real(kind=r8) :: ENDL, ENDR, KAPPA, CONDM, norm(3)
    integer :: LANMAX, worksize, ierr, nlanczos, nconverged, imessage
    real(kind=r8),allocatable :: workarray(:), eigenvalue(:), errbnd(:), &
        dummycoeff(:),dummynewcoeff(:), eigenvector(:), dyninitvec(:,:)
    integer :: ijob,ierror,EV,icoord,i

    call master_print_message("Enter dynamic property calculation subroutine")
    
    call GetDimSym
    
    LANMAX = MIN(maxnlancs,dimN)
    worksize = 6*dimN+1+4*LANMAX+LANMAX*LANMAX
    
    ! allocate workarray
    if(myid==0) then
        allocate(workarray(worksize))
        allocate(eigenvalue(LANMAX))
        allocate(errbnd(LANMAX))
        allocate(HDIAG(dimN))
        allocate(eigenvector(1*dimN))
        allocate(dyninitvec(dimN,3))
    end if

    ! construct the complementary operator matrix
    if(opmethod=="comple") then
        call Complement(operamatbig1,bigcolindex1,bigrowindex1,Lrealdim,Rrealdim,subM,0)
    end if
    
    ! Get the diagonal element of hamiltonian
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
    
    call dyn_init(dimN, dyninitvec, norm)
    
    if(myid==0) write(*,*) "dyn_prop norm= ", norm
    do icoord = 1, 3, 1
        if(myid==0) then
            workarray(1:dimN) = dyninitvec(:,icoord)
            !CALL LANDR(N,LANMAX,MAXPRS,CONDM,ENDL,ENDR,EV,KAPPA,
            !    J,NEIG,RITZ,BND,W,NW,IERR,MSGLVL)
            ierror = 0
            
            ! parameter
            ENDL = -1.0D-30
            ENDR =  1.0D-30
            imessage = 3
            CONDM = 1.0D0
            EV = 99
            KAPPA = 1.0D-6
            
            call LANDR(dimN, LANMAX, 1, CONDM, ENDL, ENDR, EV, KAPPA, &
                nlanczos, nconverged, eigenvalue, errbnd, workarray, worksize, ierror, imessage, eigenvector)
            IJOB = 0
            call MPI_BCAST(IJOB,1,MPI_integer,0,MPI_COMM_WORLD,ierr)
            write(*,*) "nlanczos = ", nlanczos
            write(*,*) "nconverged = ", nconverged
            write(*,*) "index, eigenvalue, errbound"
            do i =1, nlanczos, 1
                write(*,*) i, eigenvalue(i), errbnd(i)
            end do
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
                    else
                        stop
                    end if

                    deallocate(dummycoeff)
                    deallocate(dummynewcoeff)
                else
                    exit
                end if
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
        deallocate(HDIAG)
        deallocate(eigenvector)
        deallocate(dyninitvec)
    end if

    return

end subroutine dyn_prop_wrapper

!==========================================================
!==========================================================

subroutine dyn_init(vecdim, dyninitvec, norm)
!> master process calculate the initial wavevector
!! in the absorption spectrum 
!! \hat{r}|psi_o\rangle
    use module_sparse
    use symmetry
    use opc
    implicit none
    integer,intent(in) :: vecdim
    real(kind=r8),intent(out) ::  dyninitvec(vecdim,3), norm(3)

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

                call SubSpaceOpCdens(dynvec, domain,&
                        Lrealdim,Rrealdim,dyn_initstate,operaindex,&
                        operamatbig1,bigcolindex1,bigrowindex1,&
                        coeffIF,coeffIfcolindex,coeffIFrowindex)

                call MPI_SEND(dynvec,16*Lrealdim*Rrealdim,MPI_real8,0,iorb,MPI_COMM_WORLD,ierr)
            end if
        else
            call MPI_RECV(dynvecbuf,16*Lrealdim*Rrealdim,MPI_real8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
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
            dyninitvec(:,icoord) = dyninitvec(:,icoord) / sqrt(norm(icoord))
        end do
    end if

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

end Module dyn_prop 
