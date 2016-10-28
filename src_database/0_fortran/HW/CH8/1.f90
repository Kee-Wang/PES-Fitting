program ch8
implicit none

call area()
end

subroutine area ()
implicit none
real :: pi = 3.1415926
real :: r
real :: s
character :: name(10)


write(*,*) "Please input the radius r="
read(*,*) r
write(*,*) "Please input the aliase of this area:"
read(*,'(A)') name
s=pi*r**2

write(*,"(A10,A,F13.8)") name,'=',s

return
end 


