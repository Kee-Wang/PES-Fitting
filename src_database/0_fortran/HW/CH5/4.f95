program a4
implicit none

integer y

read (*,*) y

if (Mod(y,4)==0) then
        if (Mod(y,100)==0) then
                if(Mod(y,400)==0) then
                        write(*,*) "闰年"
                else 
                        write(*,*) "non闰年"
                end if
        else
                write(*,*) "闰年"
        end if
else
write(*,*) "non闰年"
end if

stop
end
