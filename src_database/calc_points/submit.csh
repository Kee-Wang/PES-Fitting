#!/bin/csh -f
set nmax=24
set i=1
while ( $i <= $nmax )
  set tail=$i
  if ( $i < "10" ) then
    set tail=0$i
  endif
  set dir="pts"$tail
  set inpfile="fort.11."$tail
 
  mkdir $dir
  mv 'fort.11.'$tail ./$dir
  cp 'que.pts' ./$dir
  cp 'generate.x' ./$dir

  cd $dir
  mv 'fort.11.'$tail 'fort.11'
  qsub que.pts
  cd ..

  set  i=`expr $i + 1 `
end

