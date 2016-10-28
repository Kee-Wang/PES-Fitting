program a2
implicit none

integer::i,j
integer::ans=0
j=1

write(*,*) "1"

do i=1,99,2
        ans=ans+i
        j=j+2
        write(*,*) i       
end do

write (*,*) "=",ans

stop
end
