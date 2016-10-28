program co2h2o


  integer Ndata
  character fname*12,fileindex*3, line1*99,line2*45, order*2, symb(6)
  dimension cart(6,3)
  character(len=10) ndat

Ndata=2
ndat='sample.66'

open(unit=11,file=ndat)
  do n=1,Ndata
     read(11,*) natm
     do k=1,natm
        read(11,*) symb(k),cart(k,:)
     end do
!     cart=cart*0.5291772083d0
    

     open(27,file='mp.inp',status='unknown')
     write(27,*)'*** co2h2o'
     write(27,*)'memory,100,m;'
     write(27,*)'basis=avtz'
     write(27,*)
     write(27,*)'geomtyp=xyz'
     write(27,*)'geometry={'
     write(27,*)'6'
     write(27,*)'CO2 H2O'
120  format(I8,A15)
121  format(A3,F14.9,F14.9,F14.9)
     write(27,121)' C ',cart(1,1),cart(1,2),cart(1,3)
     write(27,121)' O ',cart(2,1),cart(2,2),cart(2,3)
     write(27,121)' O ',cart(3,1),cart(3,2),cart(3,3)
     write(27,121)' O ',cart(4,1),cart(4,2),cart(4,3)
     write(27,121)' H ',cart(5,1),cart(5,2),cart(5,3)
     write(27,121)' H ',cart(6,1),cart(6,2),cart(6,3)
     write(27,121)'}'
     write(27,*)
     write(27,*)'{rhf}'
     write(27,*)''
!     write(27,*)'complexenrg = energy'
!     write(27,*)'table,complexenry'
     write(27,*)'---'
     close(27)

     call system('./sub-regnew mp.inp 1 4')
     open(28,file='mp.ut',status='old')

300  read(28,'(A99)')line1

     if (line1(2:8).eq.'? Error') then
        line2=line1(1:45)
        write(21,'(2I5,A45)')n,line2
        goto 301
     else
        if (line1(9:29).eq.'COMPLEXENRY') then
           read(28,*) COMPLEXENRY
           write(21,'(I5)') natm
           write(21,'(F15.8)') COMPLEXENRY
           do k = 1,natm
              write(21,'(A,3F15.8)') symb(k),cart(k,:)
           end do
           goto 301
        else
           goto 300
        end if
     end if

301  continue
     close(28)
  end do

  close(11)
  close(21)

end
