#!/bin/csh -f                                                                                                                  
set nmax = 21
set i = 2
while ( $i <= $nmax )
rsh node$i rm -r  /scratch/kee/*
set i = `expr $i + 1`

end

#rm fort.11
