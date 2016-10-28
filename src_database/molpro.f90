program ch3ohpes

implicit double precision (a-h,o-z)
implicit integer (i-n)

  integer Ndata
  character fname*12,fileindex*3, line1*99,line2*45, order*2, symb
  dimension cart(8,3)
  character(len=10) ndat

  call getarg(1,ndat)
  read(ndat,*) Ndata

  open(201,file='fort.1001.dat',status='old')
  open(202,file='fort.1002.dat',status='unknown')

  do i=1,Ndata

     read(201,*)
     read(201,*)
     do j=1,8
        read(201,*) symb,cart(j,1),cart(j,2),cart(j,3)
     end do
    
     do j=1,8
        write(202,'(3F15.8)') cart(j,:)
     end do

     open(27,file='mp1.inp',status='unknown')
     write(27,*)'*** CH4H2O DMS, using cartesian coor'
     write(27,*)'memory,100,m;'
     write(27,*)'basis,C=avtz,O=avtz,H=vtz'
     write(27,*)
     write(27,*)'geomtyp=xyz'
     write(27,*)'nosym'
     write(27,*)'geometry={'
     write(27,*)'8'
     write(27,*)'CH4 H2O'
120  format(I8,A15)
121  format(A3,F14.9,F14.9,F14.9)
     write(27,121)' H ',cart(1,1),cart(1,2),cart(1,3)
     write(27,121)' H ',cart(2,1),cart(2,2),cart(2,3)
     write(27,121)' H ',cart(3,1),cart(3,2),cart(3,3)
     write(27,121)' H ',cart(4,1),cart(4,2),cart(4,3)
     write(27,121)' H ',cart(5,1),cart(5,2),cart(5,3)
     write(27,121)' H ',cart(6,1),cart(6,2),cart(6,3)
     write(27,121)' C ',cart(7,1),cart(7,2),cart(7,3)
     write(27,121)' O ',cart(8,1),cart(8,2),cart(8,3)
     write(27,121)'}'
     write(27,*)
     write(27,*)'{rhf}'
     write(27,*)'{mp2;expec,dm}'
     write(27,*)
     write(27,*)'---'
     close(27)

     open(29,file='mp2.inp',status='unknown')
     write(29,*)'*** CH4 DMS, using cartesian coor'
     write(29,*)'memory,100,m;'
     write(29,*)'basis,C=avtz,H=vtz'
     write(29,*)
     write(29,*)'geomtyp=xyz'
     write(29,*)'nosym'
     write(29,*)'geometry={'
     write(29,*)'5'
     write(29,*)'CH4'
     write(29,121)' H ',cart(1,1),cart(1,2),cart(1,3)
     write(29,121)' H ',cart(2,1),cart(2,2),cart(2,3)
     write(29,121)' H ',cart(3,1),cart(3,2),cart(3,3)
     write(29,121)' H ',cart(4,1),cart(4,2),cart(4,3)
     write(29,121)' C ',cart(7,1),cart(7,2),cart(7,3)
     write(29,121)'}'
     write(29,*)
     write(29,*)'{rhf}'
     write(29,*)'{mp2;expec,dm}'
     write(29,*)
     write(29,*)'---'
     close(29)

     open(25,file='mp3.inp',status='unknown')
     write(25,*)'*** H2O DMS, using cartesian coor'
     write(25,*)'memory,100,m;'
     write(25,*)'basis,O=avtz,H=vtz'
     write(25,*)
     write(25,*)'geomtyp=xyz'
     write(25,*)'nosym'
     write(25,*)'geometry={'
     write(25,*)'3'
     write(25,*)'H2O'
     write(25,121)' H ',cart(5,1),cart(5,2),cart(5,3)
     write(25,121)' H ',cart(6,1),cart(6,2),cart(6,3)
     write(25,121)' O ',cart(8,1),cart(8,2),cart(8,3)
     write(25,121)'}'
     write(25,*)
     write(25,*)'{rhf}'
     write(25,*)'{mp2;expec,dm}'
     write(25,*)
     write(25,*)'---'
     close(25)

     call system('/usr/local/bin/molprop_2010_1_Linux_x86_64_i8 -n 8 mp1.inp')
     call system('mv mp1.inp mp1.inp.old')
     call system('/usr/local/bin/molprop_2010_1_Linux_x86_64_i8 -n 8 mp2.inp')
     call system('mv mp2.inp mp2.inp.old')
     call system('/usr/local/bin/molprop_2010_1_Linux_x86_64_i8 -n 8 mp3.inp')
     call system('mv mp3.inp mp3.inp.old')

     open(28,file='mp1.out',status='old')
     MP2success=0
300  read(28,'(A99)')line1
     if (line1(2:8).eq.'? Error') then
        line2=line1(1:45)
        write(202,'(i8,A45)')i,line2
        goto 301
     else
        if (line1(2:20).eq.'!MP2 dipole moments') then
           MP2success=1
           line2=line1(36:80)
           write(202,'(i8,A45)')i,line2
           goto 301
        else
           goto 300
        end if
     end if
301  continue
     close(28)

     open(26,file='mp2.out',status='old')
     MP2success=0
400  read(26,'(A99)')line1
     if (line1(2:8).eq.'? Error') then
        line2=line1(1:45)
        write(202,'(i8,A45)')i,line2
        goto 401
     else
        if (line1(2:20).eq.'!MP2 dipole moments') then
           MP2success=1
           line2=line1(36:80)
           write(202,'(i8,A45)')i,line2
           goto 401
        else
           goto 400
        end if
     end if
401  continue
     close(26)

     open(24,file='mp3.out',status='old')
     MP2success=0
500  read(24,'(A99)')line1
     if (line1(2:8).eq.'? Error') then
        line2=line1(1:45)
        write(202,'(i8,A45)')i,line2
        goto 501
     else
        if (line1(2:20).eq.'!MP2 dipole moments') then
           MP2success=1
           line2=line1(36:80)
           write(202,'(i8,A45)')i,line2
           goto 501
        else
           goto 500
        end if
     end if
501  continue
     close(24)

  end do

  close(201)
  close(202)

end
