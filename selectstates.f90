subroutine selectstates(valuework,dim1,valueindex,singularvalue,&
		subspacenum,syssite,szzero,szl0)
	
	use variables
	use kinds_mod
	use communicate

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
	
	if(dim1<subM) then
		call master_print_message(dim1,"dim1<subM")
		stop
	end if

	singularvalue=-1.0D0
	valueindex=0

	! use the garnet chan proposed select states rule
	if(isweep==0 .and. nelecs/=realnelecs) then
		percent=0.4
	else if(isweep==0) then
		percent=0.6
	else if(isweep==1) then
		percent=0.7
	else if(isweep==2) then
		percent=0.9
	else
		percent=1.1
	end if
	!percent=2.0+DBLE(isweep)*0.1
		
	if(percent<1.0D0) then
		directly=INT(DBLE(subM)*percent)
	else
		directly=subM
	end if
	write(*,*) "Garnet Chan's select rule,directly=",directly

	if(logic_spinreversal==0) then
		
		do i=1,dim1,1
			do j=1,directly,1
				if(valuework(i)>singularvalue(j)) then
					valueindex(j+1:directly)=valueindex(j:directly-1)
					singularvalue(j+1:directly)=singularvalue(j:directly-1)
					valueindex(j)=i
					singularvalue(j)=valuework(i)
					exit
				end if
			end do
		end do
		
	else
		do i=1,szzero+szl0,1
			do j=1,directly,1
				if(valuework(i)>singularvalue(j)) then
					if(i<=szl0 .and. j<directly) then
						valueindex(j+2:directly)=valueindex(j:directly-2)
						valueindex(j)=i
						valueindex(j+1)=szl0+szzero+i
						singularvalue(j+2:directly)=singularvalue(j:directly-2)
						singularvalue(j)=valuework(i)
						singularvalue(j+1)=valuework(i)
						exit
					else
						valueindex(j+1:directly)=valueindex(j:directly-1)
						singularvalue(j+1:directly)=singularvalue(j:directly-1)
						valueindex(j)=i
						singularvalue(j)=valuework(i)
						exit
					end if
				end if
			end do
		end do
		
		! left the last valueindex(directly) to be paired
		do while(valueindex(directly)<=szl0)
			directly=directly-1
		end do
	end if

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
				if(logic_spinreversal==0) then
					valueindex(directly+1)=j
					singularvalue(directly+1)=valuework(j)
					directly=directly+1
				else
					if(j<=szl0 .and. directly+1<subM) then
						valueindex(directly+1)=j
						valueindex(directly+2)=szl0+szzero+j
						singularvalue(directly+1)=valuework(j)
						singularvalue(directly+2)=valuework(j)
						directly=directly+2
					else
						valueindex(directly+1)=j
						singularvalue(directly+1)=valuework(j)
						directly=directly+1
					end if
				end if
				exit
			end if
		end do
		if(directly==subM) then
			exit
		end if
	end do
	end do
	
	if(logic_spinreversal/=0) then
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
				do i=subM-1,1,-1
					if(valueindex(i)>szl0 .and. valueindex(i)<=szl0+szzero) then
						valueindex(i:subM-2)=valueindex(i+1:subM-1)
						singularvalue(i:subM-2)=singularvalue(i+1:subM-1)
						exit
					end if
				end do
				do i=1,szl0,1
					done=.true.
					do j=1,subM-2,1
						if(valueindex(j)==i) then
							done=.false.
							exit
						end if
					end do
					if(done==.true.) then
						valueindex(subM-1)=i
						valueindex(subM)=i+szl0+szzero
						singularvalue(subM-1)=valuework(i)
						singularvalue(subM)=valuework(i)
						exit
					end if
				end do 
			end if
		end if
	end if

	! check if every valueindex is not 0
	do i=1,subM,1
		if(valueindex(i)==0) then
			write(*,*) "valueindex(i)==0",i
			stop
		end if
	end do
return

end subroutine selectstates
