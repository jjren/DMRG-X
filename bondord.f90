
	! two electron operator
	do i=orbstart,orbend,1
	do j=i,orbend,1
	if(bondlink(i,j)/=0 .and. myid==orbid2(i,j,1)) then
		! send the PPP operator
		if(bondlink(i,j)==2) then
			operaindex=orbid2(i,j,2)*2-1
			position1=0
			call MPI_PACK(smarowindex2(1,operaindex),subM+1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(smacolindex2(1,operaindex),smadim2,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			call MPI_PACK(operamatsma2(1,operaindex),smadim2,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
			call MPI_SEND(packbuf,packsize,MPI_PACKED,0,i+1000,MPI_COMM_WORLD,ierr)
		end if

		if(domain=='R' .and. logic_C2==0 ) then
			do j=1,2,1
				operaindex=(orbid2(i,j,2)-1)*2+j
				call SparseDirectProduct(4,4,II,IIcolindex,IIrowindex,&
								dim1,dim1,operamatsma2(:,operaindex),smacolindex2(:,operaindex),smarowindex2(:,operaindex),&
								operamatbig2(:,operaindex),bigcolindex2(:,operaindex),bigrowindex2(:,operaindex),bigdim2)
			end do
		else
			do j=1,2,1
				operaindex=(orbid2(i,j,2)-1)*2+j
				call SparseDirectProduct(dim1,dim1,operamatsma2(:,operaindex),smacolindex2(:,operaindex),smarowindex2(:,operaindex),&
								4,4,II,IIcolindex,IIrowindex,&
								operamatbig2(:,operaindex),bigcolindex2(:,operaindex),bigrowindex2(:,operaindex),bigdim2)
			end do
		end if
	end if
	end do
	end do


	if(myid==orbid2(orbadd,orbadd,1))
		if(domain=='R' .and. logic_C2==0) then
			do j=1,2,1
				operaindex=(orbid2(orbadd,orbadd,2)-1)*2+j
				call SparseDirectProduct(4,4,onesitemat(:,3*j),osmcolindex(:,j),osmrowindex(:,j),&
								dim1,dim1,IM,IMcolindex,IMrowindex,&
								operamatbig2(:,operaindex),bigcolindex2(:,operaindex),bigrowindex2(:,operaindex),bigdim2)
			end do
		else 
			do j=1,2,1
				operaindex=(orbid2(orbadd,orbadd,2)-1)*2+j
				call SparseDirectProduct(dim1,dim1,IM,IMcolindex,IMrowindex,&
								4,4,onesitemat(:,j*3),osmcolindex(:,j),osmrowindex1(:,j),&
								operamatbig2(:,operaindex),bigcolindex2(:,operaindex),bigrowindex2(:,operaindex),bigdim2)
			end do
		end if
	end if
