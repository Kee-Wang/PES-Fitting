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
    print output
    return output

order = raw_input('Please input the order of the polynomial:')
symmetry = raw_input('Please input the symmetry of the polynomial:')
train_x = raw_input('Please input the data file name: ')
#order = '5 '
#symmetry = ' 2 2 1 '
arg = order +' '+ symmetry
#train_x = 'points.dat'


cl('''
cd emsa
make
cp msa ../
cd ../
./msa '''+ arg +  '''
./postemsa.pl ''' + arg
)


f = open(train_x)
nol = 0
for line in f:
  nol=nol+1
  if nol==1:
    natom=int(line)
nconfig = nol/(natom+2)
f.close()
print('Number of atom is : ' + str (natom))
print('Number of line is : ' + str(nol))
print('Number of configuration is: '+str(nconfig))

f = open('basis.f90')
nol=1 #Num of lines in file
for line in f:
  if nol==8:
    ncoeff = int(line.split(':')[1].split(')')[0])
    ncoeff = ncoeff + 1
    break
  nol=nol+1
f.close()
print ('Number of coefficient is: ' + str(ncoeff))



a = raw_input('Do you want to continue? y/n \n')
if a == 'n':
    print('Fitting program terminated')
    quit()

g = open('fit.f90','w')
g.write('''program fit
use basis
implicit none

  external dgelss
  integer, parameter:: dp=kind(0.d0)
  real(dp) :: rmse, mue, mse, a0
  integer :: data_size
  real(dp),allocatable::xyz(:,:,:),x(:)
  real(dp),allocatable::v(:),b(:),p(:)
  real(dp),allocatable::yij(:,:), A(:,:)
  real(dp),allocatable::coeff(:),v_out(:),s(:)
  real(dp) :: work(150000), dr(3)
  integer :: ncoeff, natm, ndis
  integer :: i, j, k, m, info, rank
  character(len=32) :: data_file
  character :: symb

  natm=''' + str(natom) + '''         ! change to the number of atoms
  ncoeff=''' + str(ncoeff)+ '''       ! change to the number of coeff. (size of c in bemsa.f90)
  data_size=''' + str(nconfig)+ '''    ! change to the number of data points in pts.dat
  a0=2.0_dp

  ndis=natm*(natm-1)/2

  open(10,file='coeff.dat',status='unknown')
  open(11,FILE='points.eng',status='unknown')
  open(12,file='points.dat',status='old')

  allocate(x(ndis))
  allocate(xyz(data_size,natm,3))
  allocate(v(data_size),v_out(data_size),b(data_size),coeff(ncoeff),s(ncoeff))

  do i=1,data_size
     read(12,*)
     read(12,*) v(i)
     do j=1,natm
        read(12,*) symb,xyz(i,j,:)
     end do
  end do

  allocate(yij(data_size,ndis))

  do m=1,data_size
     k = 1
     do i=1,natm-1
        do j=i+1,natm

           yij(m,k)=0.0_dp
           dr=xyz(m,i,:)-xyz(m,j,:)
           yij(m,k)=sqrt(dot_product(dr,dr))
           yij(m,k)=yij(m,k)/0.5291772083_dp
           yij(m,k)=exp(-yij(m,k)/a0)

           k=k+1
        end do
     end do
  end do

  deallocate(xyz)
  allocate(p(ncoeff))
  allocate(A(data_size,ncoeff))

  do m=1,data_size
     x=yij(m,:)
     call bemsav(x,p)
     A(m,:)=p
  end do
  b=v

  call dgelss(data_size,ncoeff,1,A,data_size,b,data_size,s,1.0d-8,rank,work,150000,info)


  coeff(:)=b(1:ncoeff)

  do i=1,ncoeff
     write(10,*) coeff(i)
  end do

  mse=0.0_dp
  rmse=0.0_dp
  mue=0.0_dp
  do i=1,data_size
     v_out(i)=emsav(yij(i,:),coeff)
     write (11,*) v(i),v_out(i),abs(v(i)-v_out(i))*219474.63
     mse=mse+abs(v(i)-v_out(i))
     rmse=rmse+(v(i)-v_out(i))**2
     mue=mue+sqrt((v(i)-v_out(i))**2)
  end do

  mse=mse/dble(data_size)
  rmse=sqrt(rmse/dble(data_size))
  mue=mue/dble(data_size)

  print*, 'MSE  = ', mse , ' Hartree'
  print*, 'RMSE = ', rmse , ' Hartree'
  print*, 'MUE  = ', mue , ' Hartree'

  close (10)
  close (11)
  close (12)

end program
''')
g.close() #Must close the file handle if you want to compile this file.

print("Fitting...(Press Ctrl+C to interupt anytime) \n")
cl('make')
cl('./fit.x')

print ('End of program')
