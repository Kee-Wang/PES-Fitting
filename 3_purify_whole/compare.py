#!/usr/bin/env python

f = open('coeff.dat')
f2 = open('compare.dat','w')
for coef in f:
    coef=float(coef.strip())
    if abs(coef)<=10**(-7):
        coef = 0
    f2.write(str(coef)+'\n')
f.close()
f2.close()
