! Daniel Fletcher 2012

! module for increasing array sizes dynamically
! currently new indices must be larger than old

! TODO : extend to handle array size reduction

MODULE reallocate

    IMPLICIT NONE
    PUBLIC :: reallocate_2d

CONTAINS

    SUBROUTINE reallocate_1d(a,ni_new)

        INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT) :: a
        INTEGER,DIMENSION(:),ALLOCATABLE :: temp
        INTEGER,INTENT(IN) :: ni_new
        INTEGER :: ni_old

        ni_old = SIZE(a)

        ALLOCATE(temp(ni_new))

        temp(1:ni_old) = a

        CALL MOVE_ALLOC(temp,a)

    END SUBROUTINE reallocate_1d

    SUBROUTINE reallocate_2d(a,ni_new,nj_new)

        INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) :: a
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: temp
        INTEGER,INTENT(IN) :: ni_new,nj_new
        INTEGER :: ni_old,nj_old

        ni_old = UBOUND(a,1)
        nj_old = UBOUND(a,2)

        ALLOCATE(temp(ni_new,nj_new))

        temp(1:ni_old,1:nj_old) = a

        CALL MOVE_ALLOC(temp,a)

    END SUBROUTINE reallocate_2d

    SUBROUTINE reallocate_3d(a,ni_new,nj_new,nk_new)

        INTEGER,DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT) :: a
        INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: temp
        INTEGER,INTENT(IN) :: ni_new,nj_new,nk_new
        INTEGER :: ni_old,nj_old,nk_old

        ni_old = UBOUND(a,1)
        nj_old = UBOUND(a,2)
        nk_old = UBOUND(a,3)

        ALLOCATE(temp(ni_new,nj_new,nk_new))

        temp(1:ni_old,1:nj_old,1:nk_old) = a

        CALL MOVE_ALLOC(temp,a)

    END SUBROUTINE reallocate_3d

    SUBROUTINE reallocate_4d(a,ni_new,nj_new,nk_new)

        INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE,INTENT(INOUT) :: a
        INTEGER,DIMENSION(:,:,:,:),ALLOCATABLE :: temp
        INTEGER,INTENT(IN) :: ni_new,nj_new,nk_new,nl_new
        INTEGER :: ni_old,nj_old,nk_old,nl_old

        ni_old = UBOUND(a,1)
        nj_old = UBOUND(a,2)
        nk_old = UBOUND(a,3)
        nl_old = UBOUND(a,4)

        ALLOCATE(temp(ni_new,nj_new,nk_new,nl_new)

        temp(1:ni_old,1:nj_old,1:nk_old,1:nl_old) = a

        CALL MOVE_ALLOC(temp,a)

    END SUBROUTINE reallocate_4d

    
END MODULE reallocate
