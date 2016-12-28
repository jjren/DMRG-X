Program FCIDUMPsymmadaped
! this program is to created the FCIDUMP using the homemade MO-twoelectron.out created by molcas to the molpro structure FCIDUMP
! and this is the second edition to include symmetry adapted 

  implicit none

!-----------------------------------------------------
  real,parameter :: error=1.0D-17
  integer,parameter :: Ms2=0,isym=1
! Ms2 is the totalspin*2 ; 0 for singlet; 1 for doublet; 2 for triplet
! error is the if abs(integral)< error, then integral=error
!-----------------------------------------------------
  character(len=5) Group
! Group is the name of the Group D2h C2V
  integer :: nirrep,testnirrep,totalnorbs,nreadcheck,irrepprodu,noneelec,nelec
! nirrep is the number of irreps, testnirrep is check if nirrep is right
! totalnorbs is the total number of orbitals
! nreadcheck is the total number of two electron integral elements in every notzero irrep space
! irrepprdu is two irrep's product
! noneelec is the total number of one electron integral element
! nelec is the total number of electrons
  integer,allocatable :: norb(:),orbirrep(:),irrepproduexist(:)
! norb(nirrep) is in every irrep how many orbitals
! orbirrep(totalorbitals) is what the irrep of every orbital
! irrepproduexist(:) is the have exsit irrep product
  real(kind=8),allocatable :: twoelecterm(:,:,:,:),oneelecterm(:,:)
  integer :: ierror,i,j,k,quasik,l,testtotalnorbs,p,q,r,s,quasip,quasis,tmp1,tmp2,index1,ileft,leftnum,testnoneelec,zero
! leftnum is the dummy buffer contain not real integral elements
  logical :: haveexist,notzero
! haveexist is if the irrep product have exsit
! notzero is the four irrep product is not zero
  real(kind=8) :: coreenergy
! orbirrep(totalnorbs) is the irrep index of every orb

!-----------test part--------------------------------
!real(kind=8) :: testterm(2)
!integer :: p1(2),q1(2),r1(2),s1(2)
!----------------------------------------------------
  
  
  
  Group='C1'
  
  open(unit=10,file="info.inp",status="old")
  read(10,*)  Group
  read(10,*)  nirrep
  read(10,*)  totalnorbs
  read(10,*) nelec
  
  call WhatGroup(Group,testnirrep)
  
  if(testnirrep/=nirrep) then
      write(*,*) "! the number of irreps is wrong!"
  end if

  testtotalnorbs=0
  
  allocate(norb(nirrep),stat=ierror)
  if(ierror/=0) stop
  allocate(orbirrep(totalnorbs),stat=ierror)
  if(ierror/=0) stop
  allocate(twoelecterm(totalnorbs,totalnorbs,totalnorbs,totalnorbs),stat=ierror)
  if(ierror/=0) stop
  
  do i=1,nirrep,1
    read(10,*) norb(i)
    if(norb(i)/=0) then
    do j=testtotalnorbs+1,testtotalnorbs+norb(i),1
        orbirrep(j)=i
    end do
    end if
    testtotalnorbs=testtotalnorbs+norb(i)
  end do
  
  close(10)

  if(totalnorbs/=testtotalnorbs) then
    write(*,*) "----------------------------------------"
    write(*,*) "! the total number of orbitals is wrong!"
    write(*,*) "----------------------------------------"
  end if

write(*,*) "read two electron integral begain!"

!-read the two electron integral---------------------------
  open(unit=11,file="MO-twoelectron.out",status="old")
!  open(unit=12,file="test.out",status="replace")

twoelecterm=88.88D0

!-readin the nouse number-------------------------
!  do i=1,1027200,1
!     read(11,*) 
!  end do
!-------------------------------------------------

  do i=1,nirrep,1
      do j=1,i,1
          do k=1,i,1              
             if(k<i) then
                quasik=k
             else
                quasik=j
             end if
                  do l=1,quasik,1
                      
                      nreadcheck=0

                      call Ifnotzero(Group,nirrep,i,j,k,l,notzero)
                      if(notzero == .true.) then

                          do p=1,totalnorbs,1
                              if(orbirrep(p)/=k) cycle
                              do q=1,p,1
                                  if(orbirrep(q)/=l) cycle
                                !  if(i/=k) then
                                !     quasip=1
                                !  else
                                !     quasip=p
                                !  end if
                                  do r=p,totalnorbs,1
                                      if(orbirrep(r)/=i) cycle
                                        if(r==p) then
                                          do s=q,r,1
                                              if(orbirrep(s)/=j) cycle
                                              read(11,*) twoelecterm(p,q,r,s)
                                              if(abs(twoelecterm(p,q,r,s))<error) then
                                                  twoelecterm(p,q,r,s)=error
                                              end if
                                            !  write(12,*) twoelecterm(p,q,r,s),p,q,r,s
                                              nreadcheck=nreadcheck+1
                                          end do
                                       else 
                                           do s=1,r,1
                                              if(orbirrep(s)/=j) cycle
                                              read(11,*) twoelecterm(p,q,r,s)
                                              if(abs(twoelecterm(p,q,r,s))<error) then
                                                  twoelecterm(p,q,r,s)=error
                                              end if
                                            !  write(12,*) twoelecterm(p,q,r,s),p,q,r,s
                                              nreadcheck=nreadcheck+1
                                          end do
                                       end if
                                  end do
                              end do
                          end do

                        write(*,*) i,j,k,l,nreadcheck

! this part is to read the left dummy number of the integral in the buffer                      
                        if(nreadcheck/=0) then
                          leftnum=9600-Mod(nreadcheck,9600)
                          do ileft=1,leftnum,1
                              read(11,*)
                          end do
                        end if
!-------------------------------------------------------------------------           
                      end if
                  end do
    end do
    end do
    end do

close(11)
!close(12)
!------the twoelecterm read part is over-------------------------------------



open(unit=14,file="FCIDUMP",status="replace")

!! the header part
write(14,'(1X,''&FCI NORB='',I3,'',NELEC='',I2,'',MS2='',I2,'','')') totalnorbs,nelec,Ms2
write(14,'(2X,''ORBSYM='',30(I1,'',''),6(/1X,40(I1,'','')))') (orbirrep(i),i=1,totalnorbs)
write(14,'(2X,''ISYM='',I1)') isym
write(14,'(1X,''&END'')')



!------the twoelecterm write part begain----------------------------------


  allocate(irrepproduexist(nirrep),stat=ierror)
  if(ierror/=0) stop

index1=0
irrepproduexist=0

do i=1,nirrep,1
    do j=1,nirrep,1
      call irrepproduct(group,i,j,irrepprodu,nirrep)
        haveexist= .false.
        do k=1,nirrep,1
          if(irrepprodu==irrepproduexist(k)) then
            haveexist= .true.
            exit
          end if
        end do
        if(haveexist == .true.) cycle
        index1=index1+1
        irrepproduexist(index1)=irrepprodu
        
     do p=1,totalnorbs,1
       do q=1,p,1
         call irrepproduct(group,orbirrep(p),orbirrep(q),tmp1,nirrep)
         if (tmp1/=irrepprodu) cycle
         do r=1,p,1
           if(p==r) then 
             quasis=q
           else
             quasis=r
           end if
           do s=1,quasis,1
             call irrepproduct(group,orbirrep(r),orbirrep(s),tmp2,nirrep)
             if (tmp2/=irrepprodu) cycle
             if(abs(twoelecterm(r,s,p,q)-88.88D0)>1.0D-5) then 
               write(14,'(E28.20,4I4)') twoelecterm(r,s,p,q),p,q,r,s
             else
               write(14,'(E28.20,4I4)') twoelecterm(p,q,r,s),p,q,r,s
             !  write(*,*) "hello"
             end if

           end do
         end do
        end do
    end do
end do
end do

!twoelecterm write end!------------------------------------

! oneelecterm read begin

open(unit=13,file="MO-active-FOCK.out",status="old")


allocate(oneelecterm(totalnorbs,totalnorbs),stat=ierror)
if(ierror/=0) stop


read(13,*) noneelec
testnoneelec=0
zero=0
do i=1,totalnorbs,1
    do j=1,i,1
        call Ifnotzero(Group,nirrep,orbirrep(i),orbirrep(j),1,1,notzero)
        if(notzero== .true.) then
            testnoneelec=testnoneelec+1
            read(13,*) oneelecterm(i,j)
            if(abs(oneelecterm(i,j))<error) then
                oneelecterm(i,j)=error
            end if
            write(14,'(E28.20,4I4)') oneelecterm(i,j),i,j,zero,zero
        end if
    end do
end do

if(noneelec/=testnoneelec) then
    write(*,*) "-----------------------------------"
    write(*,*) "!the number of onelecterm is wrong!"
    write(*,*) "-----------------------------------"
end if
close(13)

open(unit=15,file="MO-totalcoreenergy.out",status="old")
read(15,*) coreenergy
!write(*,*) "hello",coreenergy 
write(14,'(E28.20,4I4)') coreenergy,zero,zero,zero,zero
close(15)

close(14)







!test section------------------------

!open(unit=100,file="FCIDUMPmolpro",status="old")
!open(unit=101,file="FCIDUMP",status="old")
!do i=1,6,1
!read(100,*)
!read(101,*)
!end do

!do i=1,1916996,1
!   read(100,*) testterm(1),p1(1),q1(1),r1(1),s1(1)
!   read(101,*) testterm(2),p1(2),q1(2),r1(2),s1(2)
!    if(abs(abs(testterm(1))-abs(testterm(2)))>1.0D-4 .or. p1(1)/=p1(2) .or. q1(1)/=q1(2) .or. r1(1)/=r1(2) .or. s1(1)/=s1(2) )  then
!         write(*,*) abs(testterm(1))-abs(testterm(2)),p1(1),q1(1),r1(1),s1(1),p1(2),q1(2),r1(2),s1(2)
!   end if
!end do

!close(101)
!close(100)



write(*,*) norb
write(*,*) orbirrep

!------------------------------------


write(*,*) "happy landing!"

end program
                                              
                                      




  
  

  




Subroutine WhatGroup(Group,testnirrep)
  implicit none
  character(len=5) Group
  integer :: testnirrep
  if(Group=='D2h') then
    testnirrep=8
  else if(Group=='C2v') then
      testnirrep=4
  else if(Group=="C1") then
      testnirrep=1
  else if(Group=="C2h") then
      testnirrep=4
  end if
end Subroutine

Subroutine groupbook(Group,chara,nirrep)
! the groupbook is follow the block(garnet chan) irrep order
  implicit none
  character(len=5) Group
  integer :: nirrep
  integer :: chara(nirrep,nirrep)
  if(Group=="D2h") then
      chara(1,1:8)=1

      chara(2,1)=1
      chara(2,2)=-1
      chara(2,3)=-1
      chara(2,4)=1
      chara(2,5)=-1
      chara(2,6)=1
      chara(2,7)=1
      chara(2,8)=-1

      chara(3,1)=1
      chara(3,2)=-1
      chara(3,3)=1
      chara(3,4)=-1
        chara(3,5)=-1
      chara(3,6)=1
      chara(3,7)=-1
      chara(3,8)=1

      chara(4,1)=1
      chara(4,2)=1
      chara(4,3)=-1
      chara(4,4)=-1
      chara(4,5)=1
      chara(4,6)=1
      chara(4,7)=-1
      chara(4,8)=-1

      chara(5,1)=1
      chara(5,2)=1
      chara(5,3)=-1
      chara(5,4)=-1
      chara(5,5)=-1
      chara(5,6)=-1
      chara(5,7)=1
      chara(5,8)=1

      chara(6,1)=1
      chara(6,2)=-1
      chara(6,3)=1
      chara(6,4)=-1
      chara(6,5)=1
      chara(6,6)=-1
      chara(6,7)=1
      chara(6,8)=-1

      chara(7,1)=1
      chara(7,2)=-1
      chara(7,3)=-1
      chara(7,4)=1
      chara(7,5)=1
      chara(7,6)=-1
      chara(7,7)=-1
      chara(7,8)=1

      chara(8,1)=1
      chara(8,2)=1
      chara(8,3)=1
      chara(8,4)=1
      chara(8,5)=-1
      chara(8,6)=-1
      chara(8,7)=-1
      chara(8,8)=-1
    else if(Group=="C2v") then
      chara(1,1:4)=1
      chara(2,1)=1
      chara(2,2)=-1
      chara(2,3)=1
      chara(2,4)=-1
      chara(3,1)=1
      chara(3,2)=-1
      chara(3,3)=-1
      chara(3,4)=1
      chara(4,1)=1
      chara(4,2)=1
      chara(4,3)=-1
      chara(4,4)=-1
    else if(Group=="C1") then
        chara(1,1)=1
    else if(Group=="C2h") then
        chara(1,:)=1
        chara(2,1)=1
        chara(2,2)=1
        chara(2,3)=-1
        chara(2,4)=-1
        chara(3,1)=1
        chara(3,2)=-1
        chara(3,3)=-1
        chara(3,4)=1
        chara(4,1)=1
        chara(4,2)=-1
        chara(4,3)=1
        chara(4,4)=-1
    end if

end subroutine


Subroutine Ifnotzero(Group,nirrep,i,j,k,l,notzero)
  implicit none
  character(len=5) Group
  integer :: nirrep
  integer :: chara(nirrep,nirrep)
  integer :: i,j,k,l,p
  integer :: sum1
  logical :: notzero
  
  call groupbook(Group,chara,nirrep)
    
    notzero=.true.
    do p=1,nirrep,1
        sum1=chara(i,p)*chara(j,p)*chara(k,p)*chara(l,p)
        if(sum1/=1) then
            notzero=.false.
            exit
        end if
    end do
end Subroutine

Subroutine irrepproduct(group,i,j,irrepprodu,nirrep)
  implicit none
  character(len=5) :: group
  integer :: nirrep
  integer :: i,j,irrepprodu,k,l,tmp(nirrep)
  integer :: chara(nirrep,nirrep)
  logical :: ifsearch
  
  call groupbook(group,chara,nirrep)
  
  do k=1,nirrep,1
    tmp(k)=chara(i,k)*chara(j,k)
  end do
  
  do k=1,nirrep,1
    do l=1,nirrep,1
      if(chara(k,l)/=tmp(l)) then
        ifsearch=.false.
        exit
      else
        ifsearch= .true.
      end if
    end do
    if (ifsearch == .true.) then
      irrepprodu=k
      exit
    end if
  end do
end subroutine  
  

  
    
  
  
