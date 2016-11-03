#!/bin/csh -f                                                                                                                  
#This file split a file into nfile files, and creat standard molpro input file for each configuration.
#n=number of lines per configuration
#nfile=number of file after spliting

set n=8
set nfile = 2

set i = 1
while ( $i <= $nfile )
  set tail = $i
  if ( $i < "10" ) then
    set tail = 0$i
  endif

  set k = `expr $i - 0`
  set j = `expr $k \* $n`

echo '*** co2h2o' 	 >> sample.66.temp.head
echo 'memory,100,m' 	 >> sample.66.temp.head
echo 'basis avtz' 	 >> sample.66.temp.head
echo 'geomtype = xyz'    >> sample.66.temp.head
echo 'geometry = {' 	 >> sample.66.temp.head
echo '6'		 >> sample.66.temp.head

head -n $j sample.66 | tail -$n > sample.66.temp.config
cat sample.66.temp.head  sample.66.temp.config >> sample.66.$tail
rm sample.66.temp.head

echo ''  		 >> sample.66.$tail
echo '}'		 >> sample.66.$tail
echo 'hf'		 >> sample.66.$tail
echo 'ccsd(t)-f12'	 >> sample.66.$tail
echo '---'		 >> sample.66.$tail

  set  i = `expr $i + 1`
end