#!/bin/csh -f                                                                                                                  
set nmax = 24
set i = 1
while ( $i <= $nmax )
  set tail = $i
  if ( $i < "10" ) then
    set tail = 0$i
  endif

  set k = `expr $i - 0`
  set j = `expr $k \* 4417`

  head -n $j fort.11 | tail -4417 > fort.11.$tail

  set  i = `expr $i + 1`
end

#rm fort.11
