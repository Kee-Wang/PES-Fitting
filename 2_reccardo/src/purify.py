#!/usr/bin/env python

#class purify():

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

#cl('cp basis.f90 basis_copy.f90')

"""First generate purify.f90 so that non-zero terms will be printed out"""
def purify_f90(natom=None,ncoeff=None,inter=None):

    if inter is None:
         inter = list()
         while True:
            a = input('Input r_n that is intermolecular label:')
            a=a.strip()
            if a is 'n':
                break
            else:
                inter.append(int(a))

    natom= int(natom)
    ncoeff = int(ncoeff)
    inter_mol = int(natom * (natom - int(1)) / int(2))
    prime = [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
        61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
        149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
        229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
        313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
        409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491,
        499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
        601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,
        691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
        809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887,
        907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997]

    print('Hit n as the ending:')
    f1 = open('purify_test.f90','w')
    f1.write("""program purify
    use basis

    real,dimension(1:"""+str(inter_mol)+""")::x
    real,dimension(0:"""+str(ncoeff-1)+""")::p
    integer::i
    \n""")

    j=0
    for i in range(1,inter_mol+1):
        if j<=len(inter)-1 and inter[j] is i:
            f1.write('x({:d})={:d}\n'.format(i,prime[j]))
            j += 1
        else:
            f1.write('x({:d})={:d}\n'.format(i,0))

    f1.write('''
    call bemsav(x,p)
    open(file='purified_coeff.dat',status='unknown',unit=11)

    do i=1,'''+str(ncoeff)+'''
    write(11,'(F20.10)') p(i-1)
    end do

    end
    ''')
    f1.close()
    return None


"""Run the purify.90 to print out all non-zero terms and write into file 'purified_coeff.dat'"""


"""Read the 'purified_coeff.dat' file and then modify the 'basis.f90'"""
def purified_coeff():
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

    j = 0
    for i in purify[::-1]:
        #print (i)
        replace(i,file='basis.f90')
        j=j+1
        print('Count: {:d}, replacing p({:d}) by d({:d}),{:d}<--list element'.format(j,i,i,purify[-j]))
    print(purify)

interd=[1,5,9,10,11,13]
natom = 6
ncoeff=5835
purify_f90(natom,ncoeff,interd)
cl('gfortran basis purify.f90')
cl('./a.out')
purified_coeff()