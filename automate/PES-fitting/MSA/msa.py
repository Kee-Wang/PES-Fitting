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

#order = raw_input('Please input the order of the polynomial:')
#symmetry = raw_input('Please input the symmetry of the polynomial:')
#train_x = raw_input('Please input the data file name: ')
order = '5 '
symmetry = ' 2 2 1 '
arg = order +' '+ symmetry
train_x = 'points.dat'


cl('''
cd src
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
print('1. Input file info:')
print('Number of atom is : ' + str (natom))
print('Number of line is : ' + str(nol))
print('Number of configuration is: '+str(nconfig)+'\n')

f = open('./src/basis.f90')
nol=1 #Num of lines in file
for line in f:
  if nol==8:
    ncoeff = int(line.split(':')[1].split(')')[0])
    ncoeff = ncoeff + 1
  if nol==24:
    nmonomial = int(line.split(':')[1].split(')')[0])
    nmonomial = nmonomial + 1
    break
  nol=nol+1
f.close()

print('2. Polynomial info:')
print('Given order: '+ order)
print('Given symmetry: '+ symmetry)
print ('Number of coefficient is: ' + str(ncoeff) +'\n')



a = raw_input('Do you want to continue? y/n \n')
#if a == 'n':
#    print('Fitting program terminated')
#    quit()

g = open('./src/fit.f90','w')
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
  
  print*, '3. Fitting is finished: '
  print*, 'MSE  = ', mse , ' Hartree'
  print*, 'RMSE = ', rmse , ' Hartree'
  print*, 'MUE  = ', mue , ' Hartree'

  close (10)
  close (11)
  close (12)

end program
''')
g.close() #Must close the file handle if you want to compile this file.

cl('''cd src
make''')
print("Fitting... (This might take time) \n")

cl('''cp ./src/fit.x ./
./fit.x '''+train_x+'''
rm fit.x''')


cl('cp ./src/basis.f90 ./')
cl('cp ./src/gradient.f90 ./')
cl('cp ./src/Makefile ./')
cl('cp ./src/getpot.f90 ./')
g = open('pes_shell.f90','w')
g.write('''module pes_shell
use basis
use gradient
implicit none

  real::coeff(1:'''+str(ncoeff)+''') ! change to number of coefficients
                     ! (size of c in bemsa.f90)
  save coeff

contains
  !==================================!
  ! read the coefficients of the PES !
  !==================================!
  subroutine pes_init()
    !::::::::::::::::::
    integer::i

    open(10,file='coeff.dat',status='old')

    do i=1,size(coeff)
       read (10,*) coeff(i)
    end do

    return
    close (10)
  end subroutine pes_init

  !====================================!
  ! Function to evaluate the potential !
  !====================================!
  function f(xyz)
    real,dimension(:,:),intent(in)::xyz
    real::f
    !::::::::::::::::::::::::::::::
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2)::x
    real,dimension(3)::dr
    real::a0  ! the same as the fitting code
    integer::i,j,k

    a0 = 2.d0

    k = 1
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr = xyz(:,i) - xyz(:,j)
          x(k) = sqrt(dot_product(dr,dr))
          k = k+1
       end do
    end do

    do i=1,size(x)
       x(i)=exp(-x(i)/a0)
    end do

    f=emsav(x,coeff)

    return
  end function f

  !===========================!
  ! function to calculate the !
  ! analytical gradient       !
  !===========================!
  function g(xyz)
    real,dimension(:,:),intent(in)::xyz
    real,dimension(size(xyz,2)*3)::g
    !::::::::::::::::::::
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2)::r,x
    real,dimension(3,size(xyz,2)*(size(xyz,2)-1)/2)::dr
    real,dimension(size(xyz,2)*3,size(xyz,2)*(size(xyz,2)-1)/2)::drdx
    real,dimension(1:''' + str(ncoeff) + ''')::p   ! change to number of popynomials
                               ! (size of p in bemsa.f90)
    real,dimension(1:'''+str(nmonomial)+''')::m  ! change to number of monomials
                              ! (size of m in bemsa.f90
    real::a0
    integer::i,j,k

    a0 = 2.d0  ! the same as in fitting code

    k = 1
    drdx = 0.d0
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr(:,k) = xyz(:,i) - xyz(:,j)
          r(k) = sqrt(dot_product(dr(:,k),dr(:,k)))

          drdx(3*i-2:3*i,k) = dr(:,k)/r(k)
          drdx(3*j-2:3*j,k) = -drdx(3*i-2:3*i,k)
          k = k+1
       end do
    end do

    do i=1,size(x)
       x(i)=exp(-r(i)/a0)
    end do

    call evmono(x,m)
    call evpoly(m,p)

    do i=1,3*size(xyz,2)
       g(i) = demsav(drdx,coeff,m,p,i)
    end do

    return
  end function g

end module pes_shell
''')
g.close()

cl('''mv '''+train_x+''' ./src''')
cl('''make getpot.x''')
cl('cp ./src/test.xyz ./')
print ('4. In order to run the test program, use command:')
print ('./getpot.x test.xyz \n')
print ('End of program \n')





