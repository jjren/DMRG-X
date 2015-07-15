subroutine C2_copy(direction)
! this subroutine is to copy the L space operator and R space operator in
! operamatbig1,operamatbig2
! the C2 like symmetry condition
! direction is 'i' 'l' 'r'

	use mpi
	use variables
	use communicate
	use module_sparse
	use blas95
	use f95_precision

	implicit none
	character(len=1) :: direction
	integer :: operaindex1,operaindex2
	integer :: Hindex,copyHindex,orbstart,orbend
	integer :: i,j,k,nonzero
	integer :: status(MPI_STATUS_SIZE)
	character(len=1),allocatable :: packbuf(:)
	integer :: packsize,position1
	integer :: ierr,error
	

	call master_print_message("enter in C2_copy subroutine")
	
	if(Lrealdim/=Rrealdim) then
		write(*,*) "=================================="
		write(*,*) "C2 copy Lrealdim/=Rrealdim failed!"
		write(*,*) "=================================="
		stop
	end if
	
	if(direction=='l' .or. direction=='i') then
		orbstart=1
		orbend=nleft+1
		Hindex=1
		copyHindex=2
	else if(direction=='r') then
		orbstart=norbs-nright
		orbend=norbs
		Hindex=2
		copyHindex=1
	end if

	packsize=(bigdim1*12+4*(4*Lrealdim+1))*3+1000

	if(myid/=0) then
		allocate(packbuf(packsize),stat=error)
		if(error/=0) stop
	end if

	do i=orbstart,orbend,1
		if(myid==orbid1(i,1)) then
			operaindex1=orbid1(i,2)
			
			! in the same process
			if(myid==orbid1(norbs-i+1,1)) then
				operaindex2=orbid1(norbs-i+1,2)
				
				do j=1,3,1
					bigrowindex1(:,operaindex2*3-3+j)=bigrowindex1(:,operaindex1*3-3+j)
					nonzero=bigrowindex1(4*Lrealdim+1,operaindex1*3-3+j)-1
					bigcolindex1(1:nonzero,operaindex2*3-3+j)=bigcolindex1(1:nonzero,operaindex1*3-3+j)
					call copy(operamatbig1(1:nonzero,operaindex1*3-3+j),operamatbig1(1:nonzero,operaindex2*3-3+j))
				end do
			else
			! in different process
				position1=0
				do j=1,3,1
					call MPI_PACK(bigrowindex1(1,operaindex1*3-3+j),(4*Lrealdim+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
					nonzero=bigrowindex1(4*Lrealdim+1,operaindex1*3-3+j)-1
					call MPI_PACK(bigcolindex1(1,operaindex1*3-3+j),nonzero,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
					call MPI_PACK(operamatbig1(1,operaindex1*3-3+j),nonzero,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
				end do
				call MPI_SEND(packbuf,position1,MPI_PACKED,orbid1(norbs-i+1,1),i,MPI_COMM_WORLD,ierr)
			end if

		else if(myid==orbid1(norbs-i+1,1)) then
			operaindex2=orbid1(norbs-i+1,2)
			
			call MPI_RECV(packbuf,packsize,MPI_PACKED,orbid1(i,1),i,MPI_COMM_WORLD,status,ierr)
			position1=0
			do j=1,3,1
				call MPI_UNPACK(packbuf,packsize,position1,bigrowindex1(1,operaindex2*3-3+j),(4*Lrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
				nonzero=bigrowindex1(4*Lrealdim+1,operaindex2*3-3+j)-1
				call MPI_UNPACK(packbuf,packsize,position1,bigcolindex1(1,operaindex2*3-3+j),nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
				call MPI_UNPACK(packbuf,packsize,position1,operamatbig1(1,operaindex2*3-3+j),nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
			end do
		end if
	end do

	if(myid==0) then
		Hbigrowindex(:,copyHindex)=Hbigrowindex(:,Hindex)
		nonzero=Hbigrowindex(4*Lrealdim+1,Hindex)-1
		Hbigcolindex(1:nonzero,copyHindex)=Hbigcolindex(1:nonzero,Hindex)
		Hbig(1:nonzero,copyHindex)=Hbig(1:nonzero,Hindex)
		
		if(logic_spinreversal/=0) then
			symmlinkbig(:,1,copyHindex)=symmlinkbig(:,1,Hindex)
		end if
	end if
	
	if(Hindex==1) then
		quantabigR=quantabigL
	else
		quantabigL=quantabigR
	end if
	
	! bond order operator
	! only be copied L to R in the last l direction sweep
	if(logic_bondorder/=0 .and. direction=='l') then
		do i=orbstart,orbend,1
		do j=i,orbend,1
			if(bondlink(i,j)/=0 .or. logic_bondorder==2) then
				! check the corresponds bond
				if(bondlink(norbs-i+1,norbs-j+1)==0) then
					write(*,*) "=========================================="
					write(*,*) "bondlink failed!",i,j,norbs-i+1,norbs-j+1
					write(*,*) "=========================================="
					stop
				end if

				if(myid==orbid2(i,j,1)) then
					operaindex1=orbid2(i,j,2)

					if(myid==orbid2(norbs-i+1,norbs-j+1,1)) then
						operaindex2=orbid2(norbs-i+1,norbs-j+1,2)
					
						do k=1,2,1
							bigrowindex2(:,operaindex2*2-2+k)=bigrowindex2(:,operaindex1*2-2+k)
							nonzero=bigrowindex2(4*Lrealdim+1,operaindex1*2-2+k)-1
							bigcolindex2(1:nonzero,operaindex2*2-2+k)=bigcolindex2(1:nonzero,operaindex1*2-2+k)
							call copy(operamatbig2(1:nonzero,operaindex1*2-2+k),operamatbig2(1:nonzero,operaindex2*2-2+k))
						end do
					else
						! in different process
						position1=0
						do k=1,2,1
							call MPI_PACK(bigrowindex2(1,operaindex1*2-2+k),(4*Lrealdim+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
							nonzero=bigrowindex2(4*Lrealdim+1,operaindex1*2-2+k)-1
							call MPI_PACK(bigcolindex2(1,operaindex1*2-2+k),nonzero,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
							call MPI_PACK(operamatbig2(1,operaindex1*2-2+k),nonzero,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
						end do
						call MPI_SEND(packbuf,position1,MPI_PACKED,orbid2(norbs-i+1,norbs-j+1,1),i,MPI_COMM_WORLD,ierr)
					end if

				else if(myid==orbid2(norbs-i+1,norbs-j+1,1)) then
					operaindex2=orbid2(norbs-i+1,norbs-j+1,2)
					call MPI_RECV(packbuf,packsize,MPI_PACKED,orbid2(i,j,1),i,MPI_COMM_WORLD,status,ierr)
					
					position1=0
					do k=1,2,1
						call MPI_UNPACK(packbuf,packsize,position1,bigrowindex2(1,operaindex2*2-2+k),(4*Lrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
						nonzero=bigrowindex2(4*Lrealdim+1,operaindex2*2-2+k)-1
						call MPI_UNPACK(packbuf,packsize,position1,bigcolindex2(1,operaindex2*2-2+k),nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
						call MPI_UNPACK(packbuf,packsize,position1,operamatbig2(1,operaindex2*2-2+k),nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
					end do
				end if
			end if
		end do
		end do
	end if

	! local spin operator
	! only be copied L to R in the last l direction sweep
	if(logic_localspin==1 .and. direction=='l') then
		do i=orbstart,orbend,1
		do j=i,orbend,1
			if(myid==orbid3(i,j,1)) then
				operaindex1=orbid3(i,j,2)

				if(myid==orbid3(norbs-i+1,norbs-j+1,1)) then
					operaindex2=orbid3(norbs-i+1,norbs-j+1,2)

					do k=1,2,1
						bigrowindex3(:,operaindex2-2+k)=bigrowindex3(:,operaindex1-2+k)
						nonzero=bigrowindex3(4*Lrealdim+1,operaindex1-2+k)-1
						bigcolindex3(1:nonzero,operaindex2-2+k)=bigcolindex3(1:nonzero,operaindex1-2+k)
						call copy(operamatbig3(1:nonzero,operaindex1-2+k),operamatbig3(1:nonzero,operaindex2-2+k))
					end do
				else
					! in different process
					position1=0
					do k=1,2,1
						call MPI_PACK(bigrowindex3(1,operaindex1-2+k),(4*Lrealdim+1),MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
						nonzero=bigrowindex3(4*Lrealdim+1,operaindex1-2+k)-1
						call MPI_PACK(bigcolindex3(1,operaindex1-2+k),nonzero,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
						call MPI_PACK(operamatbig3(1,operaindex1-2+k),nonzero,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
					end do
					call MPI_SEND(packbuf,position1,MPI_PACKED,orbid3(norbs-i+1,norbs-j+1,1),i,MPI_COMM_WORLD,ierr)
				end if
			else if(myid==orbid3(norbs-i+1,norbs-j+1,1)) then
				operaindex2=orbid3(norbs-i+1,norbs-j+1,2)
				call MPI_RECV(packbuf,packsize,MPI_PACKED,orbid3(i,j,1),i,MPI_COMM_WORLD,status,ierr)
				
				position1=0
				do k=1,2,1
					call MPI_UNPACK(packbuf,packsize,position1,bigrowindex3(1,operaindex2-2+k),(4*Lrealdim+1),MPI_integer4,MPI_COMM_WORLD,ierr)
					nonzero=bigrowindex3(4*Lrealdim+1,operaindex2-2+k)-1
					call MPI_UNPACK(packbuf,packsize,position1,bigcolindex3(1,operaindex2-2+k),nonzero,MPI_integer4,MPI_COMM_WORLD,ierr)
					call MPI_UNPACK(packbuf,packsize,position1,operamatbig3(1,operaindex2-2+k),nonzero,MPI_real8,MPI_COMM_WORLD,ierr)
				end do
			end if
		end do
		end do
	end if

	if(myid/=0) then
		deallocate(packbuf)
	end if

return
end subroutine C2_copy


