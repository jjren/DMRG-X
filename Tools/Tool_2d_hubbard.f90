module var
	implicit none
	integer :: nsite
	integer,allocatable :: indexxy(:,:),bondlink(:,:,:,:)
	real(kind=8) :: bondlength
	real(kind=8),parameter :: hopt=-2.4D0,hubbardU=11.26D0
	
end module var

program Tool_2d_hubbard
! this program prepare the 2d_hubbard integral for DMRG-X
	use var
	implicit none

	write(*,*) "input nsite:"
	read(*,*) nsite
	write(*,*) "bond length:"
	read(*,*) bondlength

	allocate(indexxy(nsite,nsite))
	allocate(bondlink(nsite,nsite,nsite,nsite))
	
	call bondnet
	call multichain
	call output
	
	deallocate(indexxy)
	deallocate(bondlink)
end program Tool_2d_hubbard

subroutine FCIDUMP
	use var
	implicit none

end subroutine FCIDUMP

subroutine output
	use var
	implicit none
	integer :: icol,irow,jcol,jrow,i
	logical :: iffind

	open(unit=10,file="coord.xyz",status="replace")
	write(10,*) nsite*nsite
	write(10,*) 
	do i=1,nsite*nsite,1
		iffind=.false.
		do icol=1,nsite,1
		do irow=1,nsite,1
			if(indexxy(irow,icol)==i) then
				write(10,'(1A,4F8.5)') 'C',DBLE(icol-1)*bondlength,DBLE(irow-1)*bondlength,0.0D0,1.0D0
				iffind=.true.
				exit
			end if
			if(iffind==.true.) exit
		end do
		end do
	end do
	close(10)
	
	open(unit=11,file="integral.inp",status="replace")
	open(unit=13,file="FCIDUMP",status="replace")
	do icol=1,nsite,1
	do irow=1,nsite,1
	do jcol=1,nsite,1
	do jrow=1,nsite,1
		if(bondlink(jrow,jcol,irow,icol)==1) then
			write(11,*) indexxy(irow,icol),indexxy(jrow,jcol),hopt
			write(13,*) hopt,indexxy(irow,icol),indexxy(jrow,jcol),0,0
		end if
	end do
	end do
	end do
	end do
	do i=1,nsite*nsite,1
		write(11,*) 0.0D0
	end do
	do i=1,nsite*nsite,1
		write(11,*) hubbardU
		write(13,*) hubbardU,i,i,i,i
	end do
	close(11)
	close(13)

	return

return
end subroutine output

subroutine multichain
	use var
	implicit none
	integer :: i,j

	do j=1,nsite,1
	do i=1,nsite,1
		if(mod(j,2)==1) then
			indexxy(i,j)=i+(j-1)*nsite
		else
			indexxy(i,j)=(nsite-i)+1+(j-1)*nsite
		end if
	end do
	end do
	return
end subroutine multichain

subroutine bondnet
	use var
	implicit none
	integer :: icol,irow,jcol,jrow
	
	bondlink=0
	do icol=1,nsite,1
	do irow=1,nsite,1
	do jcol=1,nsite,1
	do jrow=1,nsite,1
		if(jrow+jcol*nsite>irow+icol*nsite) then
			if ((icol-jcol)**2+(irow-jrow)**2==1) then
				bondlink(jrow,jcol,irow,icol)=1
			end if
		end if
	end do
	end do
	end do
	end do
	return
end subroutine bondnet




















