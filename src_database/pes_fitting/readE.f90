program main
use pes_shell
  implicit none


  real::auang=0.5291772083
  real::aucm=219474.63
  real::xx(3,6)
  character::symb(6)
  character(len=32)::filename
  integer::i,j,natm,ierr,n
  real::pot
  real::rms

  call getarg(1,filename)

open(21,status='old',file=trim(filename))
i=index(filename,'.',.true.)
filename=filename(1:i-1)
open(22,status='unknown',file=trim(filename)//".eng")

call pes_init()

n = 0
do 
    read(21,*,iostat=ierr) natm
    if (ierr<0) exit
    n = n+1
    read(21,*) 
    do i=1,natm
      read(21,*) symb(i),xx(:,i)
    end do
    xx=xx/auang

   pot=f(xx)
  write(*,*) pot

end do

end program main













    
