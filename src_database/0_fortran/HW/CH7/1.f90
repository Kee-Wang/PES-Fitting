Program a1

integer :: a(10) = (/(2*i, i=1,10)/)
integer :: b
integer :: i


b = 0

do i=1,10
b = b + a(i)
end do
write(*,*) real(b/10)

end

