Subroutine spincorrect(DavidWORK)
	use mpi
	use variables
	use blas95
	use f95_precision

	implicit none
	real(kind=8) :: DavidWORK(ngoodstates*nstate)
	real(kind=8) :: norm
	integer :: i,j
	logical :: done


	if(myid==0 .and. logic_spinreversal/=0) then
		do i=1,ngoodstates,1
			if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(i,1) .and. &
			abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(i,2)) then
				if(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2))&
				/=logic_spinreversal*((-1)**mod(nelecs/2,2))) then
					DavidWORK(i:nstate*ngoodstates:ngoodstates)=0.0D0
				end if
			else if(abs(symmlinkbig(symmlinkgood(i,1),1,1))/=symmlinkgood(i,1) .or. &
			abs(symmlinkbig(symmlinkgood(i,2),1,2))/=symmlinkgood(i,2)) then
				done=.false.
				do j=1,ngoodstates,1
					if(abs(symmlinkbig(symmlinkgood(i,1),1,1))==symmlinkgood(j,1) .and. &
					abs(symmlinkbig(symmlinkgood(i,2),1,2))==symmlinkgood(j,2)) then
						DavidWORK(j:nstate*ngoodstates:ngoodstates)=DavidWORK(i:nstate*ngoodstates:ngoodstates)&
						*DBLE(sign(1,symmlinkbig(symmlinkgood(i,1),1,1))*sign(1,symmlinkbig(symmlinkgood(i,2),1,2)))&
						*DBLE(logic_spinreversal)*((-1.0D0)**mod(nelecs/2,2))
						done=.true.
						exit
					end if
				end do
				if(done==.false.) then
					write(*,*) "-----------------------------------------------------------------------------"
					write(*,*) "in spincorrect did't find the state",i,symmlinkgood(i,1),symmlinkgood(i,2),"corrsponds",&
					symmlinkbig(symmlinkgood(i,1),1,1),symmlinkbig(symmlinkgood(i,2),1,2)
					write(*,*) "-----------------------------------------------------------------------------"
					stop
				end if
			end if
		end do

		norm=dot(DavidWORK(1:ngoodstates),DavidWORK(1:ngoodstates))
		write(*,*) "spincorrect state1 norm=",norm
		if(norm<1.0D-10) then
			write(*,*) "--------------------------"
			write(*,*) "in op norm is < 1.0D-10,caution!"
			write(*,*) "--------------------------"
		end if
		DavidWORK(1:ngoodstates)=DavidWORK(1:ngoodstates)/sqrt(norm)

! Gram-Schmit Orthogonalization
		if(nstate >= 2) then
		do i=2,nstate,1
			do j=1,i-1,1
				norm=dot(DavidWORK((i-1)*ngoodstates+1:i*ngoodstates),&
				DavidWORK((j-1)*ngoodstates+1:j*16*ngoodstates))
				DavidWORK((i-1)*ngoodstates+1:i*ngoodstates)=&
					DavidWORK((i-1)*ngoodstates+1:i*ngoodstates)-&
					norm*DavidWORK((j-1)*ngoodstates+1:j*ngoodstates)
			end do
				norm=dot(DavidWORK((i-1)*ngoodstates+1:i*ngoodstates),&
				DavidWORK((i-1)*ngoodstates+1:i*ngoodstates))
				write(*,*) "spincorrect state",i,"norm=",norm
				if(norm<1.0D-10) then
					write(*,*) "norm is < 1.0D-10,caution!"
				end if
				DavidWORK((i-1)*ngoodstates+1:i*ngoodstates)=&
					DavidWORK((i-1)*ngoodstates+1:i*ngoodstates)/sqrt(norm)
		end do
		end if
	end if

return
end subroutine

