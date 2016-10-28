program sorting
implicit none


real,allocatable::ab(:),pot(:),rms(:)
integer::N=23339
character::filename
integer:: i,j
real::switch1=0,switch2=0
real::a


allocate(ab(N))
allocate(pot(N))
allocate(rms(N))




!call getarg(1,filename)
open(11,status='old',file='timetest.eng')

do i=1,N
read(11,*) ab(i),pot(i),rms(i)
end do

do i=1,N
 do j=i,N
   if ( ab(i) - ab(j) > 0) then
   switch1=ab(i)
   switch2=rms(i)
   ab(i)=ab(j)
   rms(i)=rms(j)
   ab(j)=switch1
   rms(j)=switch2
end if
end do
end do

a=0
open(12,status='unknown',file='order')
do i=1,N
a=a+rms(i)**2
write(12,*) ab(i)*219474.63,sqrt(a/real(i)),rms(i)
end do
end

