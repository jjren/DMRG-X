subroutine checkinfo(info)

    implicit none

    integer :: info

    if(info/=0) then
        write(*,*) "===================="
        write(*,*) "info/=0 failed!",info
        write(*,*) "===================="
        stop
    end if
return

end subroutine checkinfo

    
