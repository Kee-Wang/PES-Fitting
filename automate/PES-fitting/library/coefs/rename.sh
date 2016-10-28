#!/bin/bash

for x in *.out
do
  base=${x%.out}
  mv $x $base.dat
done
