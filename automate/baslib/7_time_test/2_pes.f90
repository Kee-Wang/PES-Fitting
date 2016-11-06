program main
use pes_shell
implicit none

  character (len=1) :: symb
  integer :: i, j, k
  real :: pot,start,finish,avgtime,xx(3,6)
  real :: t(10)
  character(len=32)::fname
  real,parameter::auang=0.5291772083d0


fname='40000.inp'
call pes_init()


write(*,*) 'Program begin'
do k=1,10
call cpu_time(start)

open(12,file=trim(fname),status='old')
do j=1,40000
  read (12,*) 
  read (12,*) 

  do i=1,6
     read (12,*)symb, xx(:,i)
  end do
  xx=xx/auang
  pot=f(xx)
  end do

close (12)
call cpu_time(finish)
t(k)=finish-start

write(*,*) 'Finished cycle:',k,'Time elpased:',t(k)
end do

avgtime=0
do i=1,10
avgtime=avgtime+t(i)
end do
write(*,*) 'Averaged time:',avgtime/10.0
write(*,*) 'Time test finished!'
end program
