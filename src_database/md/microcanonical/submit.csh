#!/bin/csh -f
set nmax=9
set i=0
while ( $i <= $nmax )
  set tail=$i
  if ( $i < "10" ) then
    set tail=0$i
  endif
  set dir="traj"$tail
 
  mkdir $dir
  cp 'ch3oh.int' ./$dir
  cp 'subtraj' ./$dir
  cp 'dyn.x' ./$dir
  cp 'coeff.dat' ./$dir

  cd $dir
  ./subtraj ch3oh.int 1 50000 10 2
  cd .. 

  set  i=`expr $i + 1 `
end

