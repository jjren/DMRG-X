subroutine countnonzero
! this subroutine is to count the nonzero element of operamatbig and Hbig

    use variables
    use communicate
    use module_sparse
    use mpi

    implicit none

    integer :: operaindex
    integer :: i,j
    integer :: ierr

    call master_print_message("print the sparse matrix infomation")

    if(myid==0) then
        write(*,*) "bigdim1=",bigdim1
        write(*,*) "smadim1=",smadim1
        write(*,*) "bigdim2=",bigdim1
        write(*,*) "smadim2=",smadim1
        write(*,*) "Hbigdim=",Hbigdim
        write(*,*) "Hsmadim=",Hsmadim
        write(*,*) "coeffIFdim",coeffIFdim
        write(*,*) "---------------------------"
        write(*,*) "HLbig",Hbigrowindex(4*subM+1,1)-1
        write(*,*) "HLsma",Hsmarowindex(subM+1,1)-1
        write(*,*) "HRbig",Hbigrowindex(4*subM+1,2)-1
        write(*,*) "HRsma",Hsmarowindex(subM+1,2)-1
        write(*,*) "coeffIF",coeffIFrowindex(4*subM+1,:)-1
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    do i=1,norbs,1
        if(myid==orbid1(i,1)) then
            operaindex=orbid1(i,2)
            write(*,'(1I4,6I10)') i,bigrowindex1(4*subM+1,operaindex*3-2:operaindex*3)-1,&
                smarowindex1(subM+1,operaindex*3-2:operaindex*3)-1
        end if
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if(logic_bondorder==1) then
        do i=1,(norbs+1)/2,1
        do j=i,(norbs+1)/2,1
            if(bondlink(i,j)/=0) then
                if(myid==orbid2(i,j,1)) then
                    operaindex=orbid2(i,j,2)
                    write(*,'(2I4,4I10)') i,j,bigrowindex2(4*subM+1,operaindex*2-1:operaindex*2)-1,&
                        smarowindex2(subM+1,operaindex*2-1:operaindex*2)-1
                end if
            end if
        end do
        end do

        do i=(norbs+1)/2+1,norbs,1
        do j=i,norbs,1
            if(bondlink(i,j)/=0) then
                if(myid==orbid2(i,j,1)) then
                    operaindex=orbid2(i,j,2)
                    write(*,'(2I4,4I10)') i,j,bigrowindex2(4*subM+1,operaindex*2-1:operaindex*2)-1,&
                        smarowindex2(subM+1,operaindex*2-1:operaindex*2)-1
                end if
            end if
        end do
        end do
    end if

return

end subroutine countnonzero




