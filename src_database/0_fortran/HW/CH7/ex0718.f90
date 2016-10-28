Program ee

implicit none
integer :: a(10) = (/ 5,3,6,4,8,7,1,9,2,10/)
integer :: i,j
integer :: t

do  i=1,9
 do  j=i+1,10
   if (a(i)<a(j)) then
        t = a(j)
        a(j) = a(i)
        a(i) = t
        end if
 end do 
end do 
write(*,*) a
end
