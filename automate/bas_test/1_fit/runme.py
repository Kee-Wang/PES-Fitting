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
    print(output)
    return output

#order = raw_input('Please input the order of the polynomial:')
#symmetry = raw_input('Please input the symmetry of the polynomial(for example x3y2z1):')
#train_x = raw_input('Please input the data file name: ')
order = '2'
symmetry = 'x7'
train_x = 'x7.dat'
arg = order +' '+ symmetry

f = open(train_x)
nol = 0
for line in f:
  nol=nol+1
  if nol==1:
    natom=int(line)
nconfig = nol/(natom+2)
f.close()
print('1. Input file info:')
print('Number of atom is : ' + str (natom))
print('Number of line is : ' + str(nol))
print('Number of configuration is: '+str(nconfig)+'\n')


repeat = ' '
i =int(order)
while i > 2:
  repeat = repeat + ' cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), & \n'
  i = i-1

f = open('fit.f90','w')
f.write('''PROGRAM fit
use pes, wp=>pes_wp
use px
implicit none
character (len=255) :: iargv, dname, fname
real (kind=wp), allocatable :: coef(:)
call getarg (1, dname)  ! coeff dir
call getarg (2, fname)  ! data file

pes2_dwt=(/0.02_wp,0.2_wp/)

call pes0_init (dir=dname)
call pes1_init (pes_'''+symmetry+'''_sysold)
call pes2_init (fname, pes_'''+symmetry+'''_nk, pes_'''+symmetry+'''_pot)
px_pcv(2:'''+order+''') = (/ &
'''
+repeat+
'''  cx_t('''+order+''', 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr) /)
allocate (coef(0:px_'''+symmetry+'''_nbase()-1))

call px_lsq (px_'''+symmetry+'''_base, coef)
call px_errf (px_'''+symmetry+'''_base, coef)
call cx_getcf (pes_'''+symmetry+'''_getcf, pes_'''+symmetry+'''_nki, pes_'''+symmetry+'''_sysnew, px_pcv, coef)
call pes2_reinit (pes_'''+symmetry+'''_pot)
call pes2_errf (pes_'''+symmetry+'''_nk, 1)

END PROGRAM fit
''')
f.close()

f=open('Makefile','w')
f.write('''FC =ifort -r8  -O # gfortran -fdefault-real-8 -fdefault-double-8 #ifort -r8 -O
FFLAGS  = -I ../pes_shell/lib/mod
LIBS    = -L../pes_shell/lib/pes-xyz -lpes -lpx -mkl=sequential


SRC1=fit.f90
OBJ1=$(SRC1:.f90=.o)


1_fit.x : $(OBJ1)
\t$(FC) $(FFLAGS) -o fit.x $(OBJ1) $(LIBS)
\trm -rf *.o *.mod
%.o : %.f90
\t$(FC) -c $(FFLAGS) $(LIBS) $<
clean:
\trm -rf *.o *.mod
''')


f.close()


cl('make')
print('Fitting...(This might take time...)')
cl('''./fit.x ../pes_shell/coef '''+train_x+'''
rm fit.x
cd ../pes_shell/coef
./readme_transform''')

print('Fitting is finishd')







