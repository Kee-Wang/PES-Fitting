program a3
implicit none


integer age
real income
real tax
real x

write (*,*) "Please enter your age"
read (*,*) age
write (*,*) "Please enter your income"
read (*,*) income
x=income

if (age<50) then
        if (income<1000) then
                tax=x*0.03
        else if (income<5000) then
                tax=1000*0.03+(x-1000)*0.10
        else 
                tax=1000*0.03+(5000-1000)*0.10+(x-5000)*0.15
        end if
else if (income<1000) then
                tax=x*0.05
else if (income<5000) then
                tax=1000*0.03+(x-1000)*0.07
else
                tax=1000*0.03+(5000-1000)*0.07+(x-5000)*0.10
 end if

write (*,*) "Tax is:",tax

stop
end

        
