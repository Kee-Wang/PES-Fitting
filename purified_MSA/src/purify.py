#!/usr/bin/env python
import subprocess
import os
import shlex
def cl(command):
    #ip::string, command line as string input
    #op::string, return value is the output of command line
    #Notice, each time when change directly, cl starts from currect directory.
    #Use three \' if you want to input multiple line
    arg = shlex.split(command)
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    #print(output)
    return output

cl('cp basis.f90 basis_copy.f90')

f = open('purified_coeff.dat')
count = 0
purify = list()
for i in f:
    i=float(i.strip())
    if abs(i-0) >= 1*10**(-10):
        print (i)
        purify.append(count)
    count += 1
f.close()

def replace(number,file='basis_copy.f90'):
    import re

    f = open(file) #Original file
    f_new = open('temp','w') #A temporary file
    for line in f:
        line = line.strip()
        #print(line)
        poly = 'p('+str(number)+')'
        doly = 'd('+str(number)+')'
        key = '.*p\('+str(number)+'\).*'
        #if re.search(key,line):
        if re.search(key, line): #If find that polynomial
            if re.search('p\(0\)',line):
                line='p=0d0'
            else:
                line=line.replace(poly,doly) #Replace poly by doly
            #print (line)  #Print out the line that has just been replaced
        f_new.write(line+'\n') #Write new lines into temporary file
    f_new.close()
    f.close()
    cl('mv temp '+file) #Replace original file by temporary file





#print(purify)

j = 0
for i in purify[::-1]:
    #print (i)
    replace(i)
    j=j+1
    print('Count: {:d}, replacing p({:d}) by d({:d}),{:d}<--list element'.format(j,i,i,purify[-j]))
print(purify)
#print ('Count: {:d}'.format(len(purify)))
#replace(509)