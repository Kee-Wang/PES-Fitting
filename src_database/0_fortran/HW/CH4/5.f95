program unitconversion
implicit none


type::distance
        real::m
        real::cm
        real::inch
end type distance

type(distance) ::a

write(*,*) "Please enter the distance in meter:"
read(*,*) a%m

a%cm = a%m*100
a%inch = a%m*39.370

write(*,*) a%m,"m=",a%cm,"cm=",a%inch,"inch"

stop
end
