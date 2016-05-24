program main
    use,intrinsic :: iso_c_binding

    implicit none
    include 'fftw3.f03'

    type(C_PTR) :: planf,planb
    ! FFTW plan: planf x to p ; planb p to x
    integer :: ngrids,ncals
    complex((C_DOUBLE_COMPLEX)),allocatable :: fx(:),momentum(:)
    integer :: ical,igrid,idummy
    real(kind=8) :: tmp

    write(*,*) "ngrids"
    read(*,*) ngrids
    write(*,*) "ncals"
    read(*,*) ncals
    
    allocate(fx(ngrids))
    allocate(momentum(ngrids))

    planf=fftw_plan_dft_1d(ngrids,fx,momentum,FFTW_FORWARD,FFTW_MEASURE)
    planb=fftw_plan_dft_1d(ngrids,momentum,fx,FFTW_BACKWARD,FFTW_MEASURE)

    open(unit=100,file="fx.out",status="old")
    open(unit=101,file="momentum.out",status="replace")
    do ical=1,ncals,1
        fx=(0.0D0,0.0D0)
        read(100,*)
        do igrid=1,ngrids,1
            read(100,*) idummy,tmp
            fx(igrid)=CMPLX(tmp,0.0D0)
        end do
        call dft(planf,fx,momentum)
        write(101,*) ical
        write(101,*) "#########################"
        do igrid=1,ngrids,1
            write(101,*) igrid,real(momentum(igrid)),aimag(momentum(igrid))
        end do
        read(100,*)
        read(100,*)
        write(101,*)
        write(101,*)
    end do
    close(100)    
    close(101)
    deallocate(fx,momentum)
    call fftw_destroy_plan(planf)
    call fftw_destroy_plan(planb)

end program

subroutine dft(plan,array1,array2)
! this subroutine do discrete fourier transform between fx and momentum
    use,intrinsic :: iso_c_binding

    implicit none
    include 'fftw3.f03'

    type(C_PTR) :: plan
    complex(kind=8) :: array1,array2

    call fftw_execute_dft(plan,array1,array2)
return
end subroutine
