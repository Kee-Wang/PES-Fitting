program ex0410

real(kind=4)::  x
real(kind=4)::  y

read(*,*) x,y 

if (x>0) then
        if (y>0) then
                write(*,*) "1"
        else if (y<0) then
                write(*,*) "4"
        else
                write(*,*) "0"
        end if
else if (x<0) then
        write(*,*) "x<0"
end if

write(*,*) "Helloo"
stop
end
