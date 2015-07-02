subroutine Store_Operator(domain)
! this subroutine is to store the operator of every site L/R space
! in fact only the 4M basis operator matrix need to be store
	use variables
	use mpi
	use communicate
	use exit_mod
	use module_sparse

	implicit none

	character(len=1) :: domain
	! local
	integer :: i,j,k
	integer :: reclength,operaindex,Hindex,recindex
	integer :: orbstart,orbend,orbref,orbnow,nsuborbs
	integer :: count1,nonzero,dim1
	! orbref is the reference orbindex 1 or norbs
	! orbnow is nleft+1 or norbs-nright
	character(len=50) :: filename,Hfilename,Hcolfilename,Hrowfilename,symmfilename,quantafilename
	logical :: alive

	integer :: thefile , status(MPI_STATUS_SIZE) ! MPI_flag
	integer(kind=MPI_OFFSET_KIND) :: offset
	integer :: ierr

	! in fortran io, reclength is the length of direct io
	! ifort use 4byte as 1 by default
	
	call master_print_message("enter in store_operator subroutine")
	
	! L+sigmaL/R+sigmaR space operator
	! MPI parallel io

	if(domain=='L') then
		orbnow=nleft+1
		write(filename,'(i5.5,a8)') orbnow,'left.tmp'
		orbstart=1
		orbend=nleft+1
		orbref=1
		nsuborbs=nleft+1
		dim1=Lrealdim
	else if (domain=='R') then
		orbnow=norbs-nright
		write(filename,'(i5.5,a9)') orbnow,'right.tmp'
		orbstart=norbs-nright
		orbend=norbs
		orbref=norbs
		nsuborbs=nright+1
		dim1=Rrealdim
	else
		call exit_DMRG(sigAbort,"domain/=L .and. domain/='R' failed!")
	end if

	call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,thefile,ierr)
	if(myid/=0) then
		do i=orbstart,orbend,1
			if(myid==orbid1(i,1)) then
				operaindex=orbid1(i,2)
				! store mat-> colindex->rowindex in CSR format
				offset=abs(i-orbref)*(bigdim1*12+(4*subM+1)*4)*3
				do j=1,3,1
					call MPI_FILE_WRITE_AT(thefile,offset,bigrowindex1(1,3*operaindex-3+j),(4*subM+1),mpi_integer4,status,ierr)
					offset=offset+(4*subM+1)*4
					nonzero=bigrowindex1(4*dim1+1,3*operaindex-3+j)-1
					call MPI_FILE_WRITE_AT(thefile,offset,operamatbig1(1,3*operaindex-3+j),nonzero,mpi_real8,status,ierr)
					offset=offset+nonzero*8
					call MPI_FILE_WRITE_AT(thefile,offset,bigcolindex1(1,3*operaindex-3+j),nonzero,mpi_integer4,status,ierr)
					offset=offset+nonzero*4
				end do
			end if
		end do
	else if(myid==0) then
		if(domain=='L') then
			Hfilename="0-left.tmp"
			Hcolfilename="0-leftcol.tmp"
			Hrowfilename="0-leftrow.tmp"
			symmfilename="symmlink-left.tmp"
			quantafilename="quantabigL.tmp"
			Hindex=1
		else if(domain=='R') then
			Hfilename="0-right.tmp"
			Hcolfilename="0-rightcol.tmp"
			Hrowfilename="0-rightrow.tmp"
			symmfilename="symmlink-right.tmp"
			quantafilename="quantabigR.tmp"
			Hindex=2
		end if
		recindex=abs(orbnow-orbref)+1
!----------------open a binary file-------------
		
		! rowindex
		reclength=subM*4+1
		inquire(file=trim(Hrowfilename),exist=alive)
		if(alive) then
			open(unit=111,file=trim(Hrowfilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=111,file=trim(Hrowfilename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		! sigmaR index =norbs is 1; norbs-1 is 2
		write(111,rec=recindex) Hbigrowindex(:,Hindex)
		close(111)

		! Hmat
		nonzero=Hbigrowindex(4*dim1+1,Hindex)-1   ! nonzero element numbers in Hbig
	
		reclength=2*Hbigdim
		inquire(file=trim(Hfilename),exist=alive)
		if(alive) then
			open(unit=101,file=trim(Hfilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=101,file=trim(Hfilename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		! sigmaR index =norbs is 1; norbs-1 is 2
		write(101,rec=recindex) Hbig(1:nonzero,Hindex)
		close(101)
		
		! colindex
		reclength=Hbigdim
		inquire(file=trim(Hcolfilename),exist=alive)
		if(alive) then
			open(unit=109,file=trim(Hcolfilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=109,file=trim(Hcolfilename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		! sigmaR index =norbs is 1; norbs-1 is 2
		write(109,rec=recindex) Hbigcolindex(1:nonzero,Hindex)
		close(109) 

!------------------------------write the symmetry matrix---
		if(logic_spinreversal/=0) then
			reclength=2*subM
			! we use integer kind=2 in this symmlink matrix
			inquire(file=trim(symmfilename),exist=alive)
			if(alive) then
				open(unit=103,file=trim(symmfilename),access="Direct",form="unformatted",recl=reclength,status="old")
			else
				open(unit=103,file=trim(symmfilename),access="Direct",form="unformatted",recl=reclength,status="replace")
			end if
			write(103,rec=recindex) symmlinkbig(:,1,Hindex)
			close(103)
		end if
!------------------------------write the quantabigL/R(4*subM,2)---
		reclength=4*subM*2
		! quantabigL is  integer(kind=4) so without *2
		inquire(file=trim(quantafilename),exist=alive)
		if(alive) then
			open(unit=107,file=trim(quantafilename),access="Direct",form="unformatted",recl=reclength,status="old")
		else
			open(unit=107,file=trim(quantafilename),access="Direct",form="unformatted",recl=reclength,status="replace")
		end if
		if(domain=='L') then
			write(107,rec=recindex) quantabigL
		else if(domain=='R') then
			write(107,rec=recindex) quantabigR
		end if
		close(107)
	end if
	call MPI_FILE_CLOSE(thefile,ierr)

!============================================================
	! store the bondorder matrix
	if(logic_bondorder==1 .and. nsuborbs==exactsite+1) then
		write(filename,'(a1,a6)') domain,'bo.tmp'
		call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(filename),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,thefile,ierr)
		if(myid/=0) then
			count1=0
			do i=orbstart,orbend,1
			do j=i,orbend,1
				if(bondlink(i,j)/=0) then
					count1=count1+1
					if(myid==orbid2(i,j,1)) then
						operaindex=orbid2(i,j,2)
						offset=(count1-1)*(bigdim2*12+(4*subM+1)*4)*2
						do k=1,2,1
							! store mat-> colindex->rowindex in CSR format
							call MPI_FILE_WRITE_AT(thefile,offset,bigrowindex2(1,2*operaindex-2+k),(4*subM+1),mpi_integer4,status,ierr)
							nonzero=bigrowindex2(4*dim1+1,2*operaindex-2+k)-1
							offset=offset+(4*subM+1)*4
							call MPI_FILE_WRITE_AT(thefile,offset,operamatbig2(1,2*operaindex-2+k),nonzero,mpi_real8,status,ierr)
							offset=offset+nonzero*8
							call MPI_FILE_WRITE_AT(thefile,offset,bigcolindex2(1,2*operaindex-2+k),nonzero,mpi_integer4,status,ierr)
							offset=offset+nonzero*4
						end do
					end if
				end if
			end do
			end do
		end if
		call MPI_FILE_CLOSE(thefile,ierr)
	end if

!=====================================================
! print operamatbig
! only used in test 
!
!	inquire(file="exact.tmp",exist=alive)
!	if(alive) then
!		open(unit=107,file="exact.tmp",status="old",position="append")
!	else
!		open(unit=107,file="exact.tmp",status="replace")
!	end if
!	do i=orbstart,orbend,1
!		if(myid==orbid1(i,1)) then
!			write(107,*) "==========================",i
!			do j=1,3,1
!				operaindex=orbid1(i,2)*3-3+j
!				do k=1,4*subM,1
!				do l=bigrowindex1(k,operaindex),bigrowindex1(k+1,operaindex)-1,1
!					if(abs(operamatbig1(l,operaindex))>relazero) then
!						write(107,*) operamatbig1(l,operaindex),k,bigcolindex1(l,operaindex)
!					end if
!				end do
!				end do
!			end do
!		end if
!		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!	end do
!	close(107)
!	if(myid==0)  then
!		inquire(file="exact0.tmp",exist=alive)
!		if(alive) then
!			open(unit=108,file="exact0.tmp",status="old",position="append")
!		else
!			open(unit=108,file="exact0.tmp",status="replace")
!		end if
!		write(108,*) "=================",Hindex,"===================="
!		do k=1,4*subM,1
!		do l=Hbigrowindex(k,Hindex),Hbigrowindex(k+1,Hindex)-1,1
!			if(abs(Hbig(l,Hindex))>relazero) then
!				write(108,*) Hbig(l,Hindex),k,Hbigcolindex(l,Hindex)
!			end if
!		end do
!		end do
!		close(108)
!	end if
!=============================================================

return

end subroutine Store_Operator


