Program a2

integer :: f(0:10)

f(0)=0
f(1)=1

do i = 2,10
  f(i)=f(i-1)+f(i-2)
end do
write(*,*) f

end
