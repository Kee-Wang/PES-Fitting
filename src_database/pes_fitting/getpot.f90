program main
use pes_shell
  implicit none

!  integer,parameter::wp=selected_real_kind(12,300)
  real::auang=0.5291772083
  real::aucm=219474.63
  real::xx(3,6),pot0,diff
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

  rms=0.0
  n = 0
  do
     read(21,*,iostat=ierr) natm
     if (ierr < 0) exit
     n = n + 1
     read(21,*) pot0
     do i=1,natm
        read(21,*) symb(i),xx(:,i)
     end do
     xx=xx/auang

     pot=f(xx) ! Calculate the potential
     diff = abs(pot-pot0)*aucm
     write(22,'(2F15.8,F13.2)') pot0,pot,diff

     rms=rms+(pot0-pot)**2

  end do

  rms=sqrt(rms/real(n))*219474.63
  write(*,*) rms

end program main
