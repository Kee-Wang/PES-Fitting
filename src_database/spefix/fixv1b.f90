Program FixEnerg
integer,parameter::wp=selected_real_kind(12,300)
character(100):: line
character:: F*10
integer::N=15144    !number of configurations
integer:: nl=8    !number of lines perconfiguraiton
real(kind=wp),dimension(15144):: co2h2o
real(kind=wp):: h2o
real(kind=wp):: co2

open(unit=11,file='co2h2o.dmd',status='old')
do i=1,N
  read(11,*) 
  read(11,'(F15.8)') co2h2o(i)
  read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
read(11,*)
end do
close(unit=11)

open(unit=12, file='monomer.txt',status='old')
read(12,*)
read(12,'(F15.8)') co2
read(12,*)
read(12,'(F15.8)') h2o
close(unit=12)

open(unit=11,file='co2h2o.dmd',status='old')
open(unit=14,file='v1b',status='unknown')

do i=1,N,1 
read(11,*)
read(11,*)
write(14,*)  '6'
write(14,'(2F15.8)') co2h2o(i)-co2-h2o

do j=1,nl-2,1
read(11,'(A)') line
write(14,'(A80)')  line
end do 
end do
!write(*,*) co2h2o

end
