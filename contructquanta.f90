Subroutine System_Constructquanta(domain,dim1,cap_quantabig,cap_quantasma)

! transfer the good quantum number in M basis to 4M basis

    use variables
    use exit_mod
    use communicate

    implicit none
    character(len=1) :: domain
    integer :: dim1
    integer(kind=i4) :: cap_quantabig(4*dim1,2),cap_quantasma(dim1,2)

    call master_print_message("enter system_constructquanta subroutine!")
    
!   write(*,*) "cap_quantasma in contructquanta"
!   write(*,*) "myid==",myid
!   write(*,*) cap_quantasma(:,2)

    if(domain=='R' .and. logic_C2==0) then
    !    R+sigmaR space
        cap_quantabig(1:4*dim1:4,1)=cap_quantasma(1:dim1,1)
        cap_quantabig(2:4*dim1:4,1)=cap_quantasma(1:dim1,1)+1
        cap_quantabig(3:4*dim1:4,1)=cap_quantasma(1:dim1,1)+1
        cap_quantabig(4:4*dim1:4,1)=cap_quantasma(1:dim1,1)+2

        cap_quantabig(1:4*dim1:4,2)=cap_quantasma(1:dim1,2)
        cap_quantabig(2:4*dim1:4,2)=cap_quantasma(1:dim1,2)+1
        cap_quantabig(3:4*dim1:4,2)=cap_quantasma(1:dim1,2)-1
        cap_quantabig(4:4*dim1:4,2)=cap_quantasma(1:dim1,2)

    else if(domain=='L') then
    !    L+sigmaL space
        cap_quantabig(1:dim1,1)=cap_quantasma(1:dim1,1)
        cap_quantabig(dim1+1:2*dim1,1)=cap_quantasma(1:dim1,1)+1
        cap_quantabig(2*dim1+1:3*dim1,1)=cap_quantasma(1:dim1,1)+1
        cap_quantabig(3*dim1+1:4*dim1,1)=cap_quantasma(1:dim1,1)+2

        cap_quantabig(1:dim1,2)=cap_quantasma(1:dim1,2)
        cap_quantabig(dim1+1:2*dim1,2)=cap_quantasma(1:dim1,2)+1
        cap_quantabig(2*dim1+1:3*dim1,2)=cap_quantasma(1:dim1,2)-1
        cap_quantabig(3*dim1+1:4*dim1,2)=cap_quantasma(1:dim1,2)

    else if(domain=='R') then
    !    Rreverse+sigmaR space
        cap_quantabig(1:dim1,1)=cap_quantasma(1:dim1,1)
        cap_quantabig(dim1+1:2*dim1,1)=cap_quantasma(1:dim1,1)+1
        cap_quantabig(2*dim1+1:3*dim1,1)=cap_quantasma(1:dim1,1)+1
        cap_quantabig(3*dim1+1:4*dim1,1)=cap_quantasma(1:dim1,1)+2

        cap_quantabig(1:dim1,2)=cap_quantasma(1:dim1,2)
        cap_quantabig(dim1+1:2*dim1,2)=cap_quantasma(1:dim1,2)+1
        cap_quantabig(2*dim1+1:3*dim1,2)=cap_quantasma(1:dim1,2)-1
        cap_quantabig(3*dim1+1:4*dim1,2)=cap_quantasma(1:dim1,2)
    else
        call exit_DMRG(sigAbort,"domain/=L .and. domain/=R failed!")
    end if

return
end subroutine system_constructquanta

