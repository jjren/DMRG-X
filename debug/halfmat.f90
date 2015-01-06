subroutine halfmat
! this subroutine is used to print the half space operator matrix and HL
! and HR.operamatsma,operamatbig
! only used in debug mode
	use mpi
	use variables
	implicit none

	integer :: operaindex
	integer :: i,j,k,l

	if(myid==0) then
	write(*,*) "HbigL"
	do i=1,4*Lrealdim,1
	do j=1,4*Lrealdim,1
		if(abs(Hbig(j,i,1))>1.0D-2) then
			write(*,*) Hbig(j,i,1),j,i
		end if
	end do
	end do

	write(*,*) "HbigR"
	do i=1,4*Rrealdim,1
	do j=1,4*Rrealdim,1
		if(abs(Hbig(j,i,2))>1.0D-2) then
			write(*,*) Hbig(j,i,2),j,i
		end if
	end do
	end do

	write(*,*) "HsmaL"
	do i=1,Lrealdim,1
	do j=1,Lrealdim,1
		if(abs(Hsma(j,i,1))>1.0D-2) then
			write(*,*) Hsma(j,i,1),j,i
		end if
	end do
	end do

	write(*,*) "HsmaR"
	do i=1,Rrealdim,1
	do j=1,Rrealdim,1
		if(abs(Hsma(j,i,2))>1.0D-2) then
			write(*,*) Hsma(j,i,2),j,i
		end if
	end do
	end do
	end if


	do i=1,nleft+1,1
		if(myid==orbid(i)) then
			write(*,*) "orb=",i
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			do l=1,3,1
				write(*,*) "opera=",l
				do j=1,4*Lrealdim,1
				do k=1,4*Lrealdim,1
					if(abs(operamatbig(k,j,(operaindex-1)*3+l))>1.0D-2) then
						write(*,*) operamatbig(k,j,(operaindex-1)*3+l),k,j
					end if
				end do
				end do
			end do
		end if
	end do

	do i=norbs,norbs-nright,-1
		if(myid==orbid(i)) then
			write(*,*) "orb=",i
			if(mod(i,nprocs-1)==0) then
				operaindex=i/(nprocs-1)
			else
				operaindex=i/(nprocs-1)+1
			end if
			do l=1,3,1
				write(*,*) "opera=",l
				do j=1,4*Rrealdim,1
				do k=1,4*Rrealdim,1
					if(abs(operamatbig(k,j,(operaindex-1)*3+l))>1.0D-2) then
						write(*,*) operamatbig(k,j,(operaindex-1)*3+l),k,j
					end if
				end do
				end do
			end do
		end if
	end do

return

end subroutine
