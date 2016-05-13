subroutine Free_DMRG
    use module_sparse
    use checkmem_mod
    use basisindex_mod
    use ppp_term_mod
    use variables
    implicit none
    
    call Deallocate_basisindex
    call Deallocate_checkmem
    call Deallocate_sparsemat
    call Deallocate_loadbalance
    return

end subroutine Free_DMRG


