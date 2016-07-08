program Tool_TToverlap
    implicit none

    integer,parameter :: maxdim=8
    real(kind=8) :: onepdm(maxdim,maxdim,2,4),TT,diff
    integer :: i,j

    onepdm=0.0D0

    call Getelement("AO-transOpdm-S.out","S.inp",onepdm(:,:,:,1),maxdim)
    call Getelement("AO-transOpdm-T-1.out","T.inp",onepdm(:,:,:,2),maxdim)
    call Getelement("AO-transOpdm-T0.out","T.inp",onepdm(:,:,:,3),maxdim)
    call Getelement("AO-transOpdm-T1.out","T.inp",onepdm(:,:,:,4),maxdim)
    
   ! diff=0.0D0
   ! do i=1,maxdim,1
   ! do j=i,maxdim,1
   !     TT=(1.0D0/3.0D0*(onepdm(i,j,1,2)+onepdm(i,j,1,3)+onepdm(i,j,1,4)))
   !     write(*,'(3f13.8)') onepdm(i,j,1,1),TT,onepdm(i,j,1,1)/TT
   !     if(abs(i-j)<=1) then
   !         diff=diff+abs((onepdm(i,j,1,1)-TT)/TT)
   !     end if
   ! end do
   ! end do
   ! write(*,*) diff

   ! write(*,*) "============"

    do i=1,maxdim,1
    do j=i,maxdim,1
        write(*,'(6f13.8)') onepdm(i,j,1:2,2:4)
    end do
    end do
end program Tool_TToverlap


subroutine Getelement(fname,inpfname,onepdm,maxdim)
    implicit none

    character(len=*),intent(in) :: fname,inpfname
    integer,intent(in) :: maxdim
    real(kind=8),intent(out) :: onepdm(maxdim,maxdim,2)

    ! local
    integer :: morb,stateindex
    real(kind=8),allocatable :: transDM0(:,:,:,:,:)
    integer :: istate,jstate,norbs,nstates,i,iorb,jorb
    integer(kind=4),allocatable :: orbindex(:)

    open(unit=10,file=inpfname,status="old")
    read(10,*) morb,stateindex
    allocate(orbindex(morb))
    do i=1,morb,1
        read(10,*) orbindex(i)
    end do
    close(10)
     
    open(unit=398,file=fname,form="unformatted",status="old")
    read(398) norbs,nstates
    allocate(transDM0(norbs,norbs,2,nstates,nstates))
    do istate=1,nstates,1
    do jstate=istate,nstates,1
        read(398) transDM0(:,:,:,jstate,istate)
    end do
    end do
    close(398)
    
    do iorb=1,morb,1
    do jorb=iorb,morb,1
        onepdm(iorb,jorb,1:2)=transDM0(orbindex(iorb),orbindex(jorb),1:2,stateindex,stateindex)
    end do
    end do

    deallocate(orbindex)
    deallocate(transDM0)
    return
end subroutine Getelement
