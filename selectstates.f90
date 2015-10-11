subroutine selectstates(valuework,dim1,valueindex,singularvalue,&
		subspacenum,syssite,szzero,szl0)
	
	use variables
	use kinds_mod

	implicit none
	integer :: dim1
	integer :: szzero,szl0
	real(kind=r8) :: valuework(dim1),singularvalue(subM)
	integer :: valueindex(subM)
	real(kind=r8) :: percent
	integer :: i,j,m,k
	integer :: directly,syssite
	! directly is the num of states we have selected
	integer :: subspacenum((syssite*2+1)*(syssite*2+3)+1)
	logical :: noequal,done,ifexist,iffind

	write(*,*) "enter in selectstates subroutine!"
	
	singularvalue=-1.0D0
	valueindex=0

	if(logic_spinreversal==0) then
		
		do i=1,dim1,1
			do j=1,subM,1
				if(valuework(i)>singularvalue(j)) then
					valueindex(j+1:subM)=valueindex(j:subM-1)
					singularvalue(j+1:subM)=singularvalue(j:subM-1)
					valueindex(j)=i
					singularvalue(j)=valuework(i)
					exit
				end if
			end do
		end do
		
	! use the garnet chan proposed select states rule
	!	if(isweep==0 .and. nelecs/=realnelecs) then
	!		percent=0.5
	!	else if(isweep==0) then
	!		percent=0.7
	!	else if(isweep==1) then
	!		percent=0.9
	!	else
	!		percent=1.1
	!	end if
		percent=2.0+DBLE(isweep)*0.1
		
		if(percent<1.0D0) then
			directly=INT(DBLE(subM)*percent)
			write(*,*) "Garnet Chan's select rule,directly=",directly
			do while(directly<subM)
			do i=2,subspacenum(1)+1,1
				do j=sum(subspacenum(2:i)),sum(subspacenum(2:i-1))+1,-1
					do m=1,directly,1
						noequal=.true.
						if(j==valueindex(m)) then
							noequal=.false.
							exit
						end if
					end do
					if(noequal==.true.) then
						valueindex(directly+1)=j
						singularvalue(directly+1)=valuework(j)
						directly=directly+1
						exit
					end if
				end do
				if(directly==subM) then
					exit
				end if
			end do
			end do
		end if
	
		! check if every valueindex is not 0
		do i=1,subM,1
			if(valueindex(i)==0) then
				write(*,*) "valueindex(i)==0",i
				stop
			end if
		end do
	else
		do i=1,szzero+szl0,1
			do j=1,subM,1
				if(valuework(i)>singularvalue(j)) then
					if(i<=szl0 .and. j<subM) then
						valueindex(j+2:subM)=valueindex(j:subM-2)
						valueindex(j)=i
						valueindex(j+1)=szl0+szzero+i
						singularvalue(j+2:subM)=singularvalue(j:subM-2)
						singularvalue(j)=valuework(i)
						singularvalue(j+1)=valuework(i)
						exit
					else
						valueindex(j+1:subM)=valueindex(j:subM-1)
						singularvalue(j+1:subM)=singularvalue(j:subM-1)
						valueindex(j)=i
						singularvalue(j)=valuework(i)
						exit
					end if
				end if
			end do
		end do
		
		! check if every valueindex is not 0
		do i=1,subM-1,1
			if(valueindex(i)==0) then
				write(*,*) "valueindex(i)==0",i
				stop
			end if
		end do
		
		! check last index
		if(valueindex(subM)<=szl0) then
			iffind=.false.
			do i=szl0+szzero,szl0+1,-1
				do j=1,subM-1,1
					if(valueindex(j)==i) then
						done=.false.
						exit
					else 
						done=.true.
					end if
				end do
				if(done==.true.) then
					valueindex(subM)=i
					singularvalue(subM)=valuework(i)
					iffind=.true.
					exit
				end if
			end do
			if(iffind==.false.) then
				write(*,*) "-----------------------------------------------------"
				write(*,*) "selectstates did not find the last index valueindex=0"
				write(*,*) "-----------------------------------------------------"
				stop
			end if
		end if
	end if

return

end subroutine selectstates
