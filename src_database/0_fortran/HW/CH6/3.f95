program a3
implicit none

real weight
real::key=45.0
real::e=0.001
integer i



outter:  do i=1,5
        write(*,*) "Enter body weight:"
        read(*,*) weight
        if ((abs(weight-key))<=e) then
                write(*,*) "Yes, that's right."
                exit outter
        else
        Write(*,*) "Trial left:",5-i
        end if
end do outter

stop
end
