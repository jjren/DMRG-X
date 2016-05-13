subroutine Free_Program
    use variables
    use PPP_term_mod
    implicit none
    deallocate(pppV,pppVlink)
    if(allocated(pppw)) then
        deallocate(pppw)
    end if
    return
end subroutine Free_Program
