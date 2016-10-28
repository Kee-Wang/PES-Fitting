program a4
implicit none


real::ans=0
real::prod=1
integer::x=10
integer i
integer j

do i=1,10
        do j=1,i
                prod=prod*j
                write(*,*) "i:",i,"j:",j
        end do
        ans=1/prod+ans
end do

write (*,*) "Final answer:",ans

stop
end
