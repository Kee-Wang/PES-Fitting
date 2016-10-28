#!/bin/csh -f
set nmax=19
set i=0
while ( $i <= $nmax )
  set tail=$i
  if ( $i < "10" ) then
    set tail=0$i
  endif
  set dir="traj"$tail
#  set inpfile="fort.1001.dat."$tail
#  set  inpfile=$tail.int
#  echo $dir
#  echo $i
 
  mkdir $dir
  mv 'ch3oh.int.'$tail ./$dir
  cp 'subtraj' ./$dir
  cp 'dynam.x' ./$dir
#  cp 'copy.csh' ./$dir
  cp 'coeff.34369.dat' ./$dir

  cd $dir
  mv 'ch3oh.int.'$tail 'ch3oh.int'
  ./subtraj ch3oh.int 1 5000 5 2
  cd .. 

#  cp $root/$origindir $root/$dir -r
#  rm $root/$dir -r
#  mv $root/$dir/"CD3CHO-5Mprp61.int" $root/$dir/"CD3CHO-5Mprp"$inpfile
  set  i=`expr $i + 1 `
end
