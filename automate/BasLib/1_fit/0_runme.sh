#!/bin/bash
path='../pes_shell/coef'

make
./fit.x $path *.fitting
rm fit.x

cd $path
for x in *.out
do
  base=${x%.out}
  mv $x $base.dat
done

echo 'Fitting PIP is ready to call!'
