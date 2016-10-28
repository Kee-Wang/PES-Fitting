program main
implicit none
integer::i,j,natm
real::shift(3),R0,R,diagI(3)
character(len=25)::fname
double precision, parameter::auang=0.5291772083d0
real,dimension(:),allocatable::mass
real,dimension(:,:),allocatable::x,y
character,dimension(:),allocatable::sym

  call getarg(1,fname)
  open(1,file=trim(fname))
  open(10,file="molden.inp",status='unknown')
  open(11,file="H2CO.out",status='unknown')
  read(1,*) natm

  allocate(x(3,natm))
  allocate(y(3,natm))
  allocate(mass(natm))
  allocate(sym(natm))

  read(1,*)
  do i=1,natm
     read(1,*) sym(i),x(:,i)
     selectcase(sym(i))
        case("H")
           mass(i)=1837.152582
        case("C")
           mass(i)=21874.66176
        case("O")
           mass(i)=29156.94556
     endselect
  enddo

  write(11,'(I4)') natm
  write(11,*) "Input geometry"
  do i=1,natm
     write(11,'(A,3F15.8)') sym(i),x(:,i)
  end do

  call com_coor(natm,x,mass)
  write(11,'(I4)') natm
  write(11,*) "Move to center-of-mass"
  do i=1,natm
     write(11,'(A,3F15.8)') sym(i),x(:,i)
  end do

  call prin_axis_coor(natm,x,mass,diagI)

  do i=1,natm
     y(3,i) = x(1,i)   
     y(1,i) = x(2,i)
     y(2,i) = x(3,i)
  end do

  write(11,'(I4)') natm
  write(11,*) "Rotate to principle axis"
  do i=1,natm
     write(11,'(A,3F15.8)') sym(i),y(:,i)!/0.5291772083
  enddo

end
