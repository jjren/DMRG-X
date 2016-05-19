program main
    use,intrinsic :: iso_c_binding

    implicit none
    include 'fftw3.f03'

    type(C_PTR) :: planf,planb
    ! FFTW plan: planf x to p ; planb p to x
    integer :: ngrids,ncals
    complex(kind=8),allocatable :: fx(:),momentum(:)
    integer,parameter :: nhalf=75,nboundary=50
    integer :: ical,i,idummy
    real(kind=8) :: tmp

    write(*,*) "ngrids"
    read(*,*) ngrids
    write(*,*) "ncals"
    read(*,*) ncals
    ngrids=ngrids+nboundary*2
    allocate(fx(ngrids))
    allocate(momentum(ngrids))

    planf=fftw_plan_dft_1d(ngrids,fx,momentum,FFTW_FORWARD,FFTW_MEASURE)
    planb=fftw_plan_dft_1d(ngrids,momentum,fx,FFTW_BACKWARD,FFTW_MEASURE)

    open(unit=100,file="ddfx.out",status="old")
    open(unit=101,file="ddmomentum.out",status="replace")
    do ical=1,ncals,1
        read(100,*)
        do i=1,nhalf,1
            read(100,*)
        end do
        fx=(0.0D0,0.0D0)
        do i=1,nhalf,1
            read(100,*) idummy,tmp
            if(i<=ngrids-2*nboundary) then
                fx(i+nboundary)=CMPLX(tmp,0.0D0)
            end if
        end do
        !write(*,*) fx
        momentum=(0.0D0,0.0D0)
        call dft(planf,fx,momentum)
        write(101,*) ical
        write(101,*) "#########################"
        do i=1,ngrids,1
            write(101,*) i,real(momentum(i)),aimag(momentum(i))
        end do
        write(101,*)
        write(101,*)
        read(100,*)
        read(100,*)
    end do
    close(100)    
    close(101)

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
