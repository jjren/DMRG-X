subroutine countnonzero
! this subroutine is to count the nonzero element of operamatbig and Hbig

	use variables
	use communicate

	implicit none

	integer(i4) :: nnonzero,nzero
	integer :: operaindex
	integer :: i,j,k,l
	

	if(myid==0) then
		do k=1,2,1
			nnonzero=0
			nzero=0
			do j=1,4*subM,1
			do i=1,4*subM,1
				if(abs(Hbig(i,j,k))<relazero) then
					nzero=nzero+1
				else
					nnonzero=nnonzero+1
				end if
			end do
			end do
			write(*,*) "Hbig",k,"nzero=",nzero,"nnonzero=",nnonzero
		end do
	else
		do l=1,norbs,1
			if(myid==orbid(l)) then
				if(mod(l,nprocs-1)==0) then
					operaindex=l/(nprocs-1)
				else
					operaindex=l/(nprocs-1)+1
				end if
				do k=1,3,1 
					nnonzero=0
					nzero=0
					do j=1,4*subM,1
					do i=1,4*subM,1
						if(abs(operamatbig(i,j,3*operaindex-3+k))<relazero) then
							nzero=nzero+1
						else
							nnonzero=nnonzero+1
						end if
					end do
					end do
					write(*,*) "operator",l,"matrix",k,"nzero=",nzero,"nnonzero=",nnonzero
				end do
			end if
		end do
	end if

return

end subroutine countnonzero




