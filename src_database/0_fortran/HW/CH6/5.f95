program a5
implicit none

character(len=20)::string
character(len=20)::new
integer::strlen
integer::i=1
integer::j=1



write(*,*) "Please enter a phrase:"
read(*,"(A)") string
strlen=len_trim(string)
!write (*,*) string, strlen
new="          "
do i=1,strlen

        if(ichar(string(i:i))/=32) then
                new(j:j)=string(i:i)
!                write(*,*) j,i,trim(new),"|||||",string
                j=j+1
        end if

!write (*,*) "i=",i

end do
write(*,"(A10)") trim(new)
stop 
end

