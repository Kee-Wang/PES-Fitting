program sorting
implicit none


real,allocatable::ab(:),pot(:),diff(:)
integer::N=15000
character::filename
integer:: i,j
real::switch1=0,switch2=0
real::rms


allocate(ab(N))
allocate(pot(N))




open(11,status='old',file='1d_cut.energy')

do i=1,N
read(11,*) pot(i),ab(i)
end do

do i=1,N
 do j=i,N
   if ( ab(i) - ab(j) > 0) then
   switch1=ab(i)
   switch2=pot(i)
   ab(i)=ab(j)
   pot(i)=pot(j)
   ab(j)=switch1
   pot(j)=switch2
end if
end do
end do

rms=0.0
open(12,status='unknown',file='order')
do i=1,N
write(12,*) pot(i),ab(i)
end do
end

