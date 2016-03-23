subroutine Free_DMRG
    use module_sparse
    use checkmem_mod
    use basisindex_mod
    implicit none
    
    call Deallocate_basisindex
    call Deallocate_checkmem
    call Deallocate_sparsemat

    return

end subroutine Free_DMRG
