program test

use pes_shell
use xxralpha

implicit none
double precision::energy
call pes_init()


files(1)='input'
units(1)=11
open(unit=units(1),file=files(1),status='old')

call readcart(units(1),xx)

xx = xx/auang
energy = f(xx) * hatreecm

write(*,*) 'Energy is:', energy, 'cm-1'

end program
