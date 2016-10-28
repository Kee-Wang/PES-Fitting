#!/bin/csh -f
set nmax=34
set i=30
while ( $i <= $nmax )
  set tail=$i
  if ( $i < "10" ) then
    set tail=0$i
  endif
  set dir="traj_min."$tail
 
  mkdir $dir
  cp 'minimum.int' ./$dir
  cp 'subtraj' ./$dir
  cp 'ext_mpeg.pl' ./$dir
  cp 'mp_eg.tmp' ./$dir
  cp 'dmd.x' ./$dir

  cd $dir
  ./subtraj minimum.int 1 3500 10 168
  cd .. 

  set  i=`expr $i + 1 `
end

