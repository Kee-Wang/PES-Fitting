program a1
implicit none

real x
real tax

read(*,*) x

if (x<1000) then
	tax=0.03*x
	
else if(x<5000) then
	tax=0.03*1000+(x-1000)*0.1

else 
	tax=0.03*1000+(5000-1000)*0.1+(x-5000)*0.15
end if

write(*,"(f10.2)") tax

stop
end
