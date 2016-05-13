Subroutine ReadIntegral
    use variables
    use kinds_mod

    implicit none
    ! local
    integer :: i
    integer :: link1,link2
    real(kind=r8) :: dummyt

    open(unit=14,file="integral.inp",status="old")

    ! be careful about the structure of the integral formatted
    t=0.0D0
    bondlink=0
    do i=1,nbonds,1
        read(14,*) link1,link2,dummyt
        if(abs(dummyt)>hopthresh) then
            bondlink(link1,link2)=1  ! if linked , bondlink=1
            bondlink(link2,link1)=1
            t(link1,link2)=dummyt
            t(link2,link1)=dummyt
        end if
    end do

    do i=1,norbs,1
        bondlink(i,i)=2    ! bondlink=2 means that it is the site energy
        read(14,*) t(i,i)  ! t(i,i) is the site energy
    end do
    
    hubbardU=0.0D0
    do i=1,norbs,1
        read(14,*) hubbardU(i)  ! read in every hubbard U
    end do
    
!   write(*,*) "--------------------------------------------------"
!   write(*,*) "in the QC-DMRG case readin the FCIDUMP integrals"
!   write(*,*) "--------------------------------------------------"
    
    close(14)
    return

end Subroutine ReadIntegral
