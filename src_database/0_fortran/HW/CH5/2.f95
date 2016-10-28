program a2
implicit none

integer week

write(*,*) "Please insert which day of the week"
read (*,*) week


select case(week)

case(1,4)
        write(*,*) "News"
case(2,5)
        write(*,*) "TV"
case(3,6)
        write(*,*) "Cartoon"
case(7)
        write(*,*) "Movie"
case default
        write(*,*) "None exist, alien!"
end select      
stop
end
