program main
use pes_shell
  implicit none

  real(kind=wp)::auang=0.5291772083
  real(kind=wp)::aucm=219474.63
  real(kind=wp)::xx(3,5) ! change
  character(len=2)::symb
  character(len=32)::filename
  integer::i,j,natm,ierr
  real(kind=wp)::pot,pot0

  call getarg(1,filename)
  open(21,status='old',file=trim(filename))

  call pes_init()

  do
     read(21,*,iostat=ierr) natm
     if (ierr < 0) exit
     read(21,*) pot0
     do i=1,natm
        read(21,*) symb,xx(:,i)
     end do
     xx=xx/auang
     
     pot=f(xx) ! Calculate the potential

     write(*,*) pot0,pot
  end do

end program main
