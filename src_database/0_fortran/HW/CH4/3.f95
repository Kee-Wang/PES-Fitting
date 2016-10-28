program score1
implicit none

real a

write (*,*) "Please enter the student's original grade:"
read (*,*) a
write (*,*) "The curved grade is:",a**0.5*10

stop
end

