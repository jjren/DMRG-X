Module opc
!> this module do OP*C operation
    use kinds_mod
    use mathlib
    implicit none

contains
!==========================================================
!==========================================================
subroutine SubSpaceOpCdens(coeffnew,domain,&
            iLrealdim,iRrealdim,stateindex,operaindex,&
            cap_big,cap_bigcol,cap_bigrow,cap_coeff,cap_coeffcol,cap_coeffrow)
!> this subroutine calculate the OP*C of local operator 
!! the output is a dense matrix 
    implicit none
    
    character(len=1),intent(in) :: domain
    integer(kind=i4),intent(in) :: iLrealdim,iRrealdim,stateindex,operaindex
    real(kind=r8),intent(in) :: cap_big(:,:),cap_coeff(:,:)
    integer(kind=i4),intent(in) :: cap_bigcol(:,:),cap_bigrow(:,:),&
                                cap_coeffcol(:,:),cap_coeffrow(:,:)
    
    real(kind=r8),intent(out) :: coeffnew(4*iLrealdim,4*iRrealdim)

    
    if(domain=="L") then
        call SpMMtoSptodens('N','N',4*iLrealdim,4*iLrealdim,4*iLrealdim,4*iRrealdim,4*iLrealdim,4*iRrealdim,&
                cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),&
                cap_coeff(:,stateindex),cap_coeffcol(:,stateindex),cap_coeffrow(:,stateindex),&
                coeffnew)
    else if(domain=="R") then
        call SpMMtoSptodens('N','T',4*iLrealdim,4*iRrealdim,4*iRrealdim,4*iRrealdim,4*iLrealdim,4*iRrealdim,&
                cap_coeff(:,stateindex),cap_coeffcol(:,stateindex),cap_coeffrow(:,stateindex),&
                cap_big(:,operaindex),cap_bigcol(:,operaindex),cap_bigrow(:,operaindex),&
                coeffnew)
    end if

end subroutine SubSpaceOpCdens

!==========================================================
!==========================================================
end Module opc
