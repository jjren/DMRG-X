Module basisindex_mod
    use communicate
    use kinds_mod
    use variables,only : nelecs, totalSz,subM,subMp,logic_perturbation
    implicit none
    ! nmaxgoodbasis and nmaxgoodbasisp is the maximum number of ngoodstates
    integer :: nmaxgoodbasis,nmaxgoodbasisp
    integer(kind=i4),allocatable :: goodbasis(:,:),goodbasisp(:,:),&
     goodbasiscol(:),goodbasiscolp(:)

contains
!=================================================================
!=================================================================

subroutine Allocate_basisindex
    implicit none

    allocate(goodbasis(nmaxgoodbasis,2))
    allocate(goodbasiscol(4*subM+1)) 
    ! every column(R) how many goodquantum number 
    goodbasiscol(1)=1
    
    if(logic_perturbation==1) then
        allocate(goodbasisp(nmaxgoodbasisp,2))
        allocate(goodbasiscolp(4*subMp+1))
        goodbasiscolp(1)=1
    end if

    return
end subroutine allocate_basisindex

subroutine deallocate_basisindex
    implicit none

    deallocate(goodbasis,goodbasiscol)
    if(logic_perturbation==1) then
        deallocate(goodbasisp,goodbasiscolp)
    end if

    return
end subroutine Deallocate_basisindex

!=================================================================
!=================================================================

subroutine basisindex(iLdim,iRdim,cap_quantaL,cap_quantaR,nbasis,basislink,basiscol)
    ! this subroutine do calculate the nonzero state index
    implicit none

    integer(kind=i4),intent(in) :: iRdim,iLdim,&
                         cap_quantaL(iLdim,2),cap_quantaR(iRdim,2)
    integer(kind=i4),intent(out) :: nbasis,basislink(:,:),basiscol(iRdim+1)
    ! local            
    integer :: ir,il
    
    nbasis=0
    do ir=1,iRdim,1
        basiscol(ir+1)=basiscol(ir)
        do il=1,iLdim,1
            if((cap_quantaL(il,1)+cap_quantaR(ir,1)==nelecs) .and. &
                cap_quantaL(il,2)+cap_quantaR(ir,2)==totalSz) then
                basiscol(ir+1)=basiscol(ir+1)+1
                nbasis=nbasis+1
                basislink(nbasis,1)=il
                basislink(nbasis,2)=ir
            end if
        end do
    end do

return
end subroutine basisindex

!=================================================================
!=================================================================

end module basisindex_mod
