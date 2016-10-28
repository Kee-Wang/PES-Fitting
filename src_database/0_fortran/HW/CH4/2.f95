program area1
implicit none

real,parameter ::  pi=3.1415926
real r

write(*,*) "Please enter the radius of the circle:"
read (*,*) r
write(*,*) "The radius you entered is:,",r,".","The area is:",pi*r*r

stop
end
