Program FixEnerg
integer,parameter::wp=selected_real_kind(12,300)
character(100):: line
character:: F*10
integer::N=15144    !number of configurations
integer:: nl=8    !number of lines perconfiguraiton
real(kind=wp):: co2h2o
real(kind=wp):: h2o
real(kind=wp),dimension(15144):: co2
F='(F14.10)'


open(unit=11,file='co2h2o.dmd',status='old')

do i=1,N
  read(11,*) 
  read(11,F) co2h2o(i)
  read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
end do
close(unit=11)

open(unit=12,file='h2o.dmd',status='old')
do i=1,N
  read(12,*)
  read(12,F) h2o(i)
read(12,*)
read(12,*)
read(12,*)
end do
open(unit=13,file='co2.dmd',status='old')

do i=1,N
  read(13,*)
  read(13,F) co2(i)
read(13,*)
read(13,*)
read(13,*)
end do 

open(unit=11,file='co2h2o.dmd',status='old')
open(unit=14,file='relativeE',status='unknown')

do i=1,N,1 
read(11,*)
read(11,*)
write(14,*)  '6'
write(14,F) co2h2o(i)-co2(i)-h2o(i)

do j=1,nl-2,1
read(11,'(A)') line
write(14,'(A80)')  line
end do 
end do
!write(*,*) co2h2o

end
