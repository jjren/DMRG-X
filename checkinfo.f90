subroutine checkinfo(info,char1)

    implicit none

    integer,intent(in) :: info
    character(len=6),intent(in) :: char1

    if(info/=0) then
        write(*,*) "===================="
        write(*,*) "info/=0 failed!",info
        write(*,*) "char=",char1
        write(*,*) "===================="
        stop
    end if
return

end subroutine checkinfo

    
