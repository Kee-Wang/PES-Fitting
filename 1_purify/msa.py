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


"""First generate purify.f90 so that non-zero terms will be printed out"""
def purify_f90(natom=None,ncoeff=None,inter=None):

    if inter is None:
         inter = list()
         while True:
            a = raw_input('Input r_n that is intermolecular label, end with "n":')
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

    #print('Hit n as the ending:')
    f1 = open('./src/purify.f90','w')
    f1.write("""program purify
    use basis

    real,dimension(1:"""+str(inter_mol)+""")::x
    real,dimension(0:"""+str(ncoeff-1)+""")::p
    integer::i
    \n""")

    j=0
    for i in range(1,inter_mol+1):
        if j<=len(inter)-1 and inter[j] is i:
            f1.write('x({0:d})={1:d}\n'.format(i,prime[j]))
            j += 1
        else:
            f1.write('x({0:d})={1:d}\n'.format(i,0))

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
    f = open('./src/purified_coeff.dat')
    count = 0
    purify = list()
    for i in f:
	try:
       	    i=float(i.strip())
            if abs(i-0) >= 1*10**(-10):
                print (i)
                purify.append(count)
	except:
	    purify.append(count)
        count += 1
    f.close()

    def replace(number,file='basis.f90'):
        import re

        f = open(file) #Original file
        f_new = open('./src/temp','w') #A temporary file
        for line in f:
            line = line.strip()
            #print(line)
            poly = 'p('+str(number)+')'
            doly = 'ddd('+str(number)+')'
            key = '.*p\('+str(number)+'\).*'
            #if re.search(key,line):
            if re.search(key, line): #If find that polynomial
                if re.search('p\(0\)',line):
                    line='''real,dimension(0:'''+str(ncoeff-1)+''')::ddd
				p=0d0'''
                else:
                    line=line.replace(poly,doly) #Replace poly by doly
                #print (line)  #Print out the line that has just been replaced
            f_new.write(line+'\n') #Write new lines into temporary file
        f_new.close()
        f.close()
        cl('''cd src
		mv temp basis.f90''') #Replace original file by temporary file

    j = 0
    for i in purify[::-1]:
        #print (i)
        replace(i,file='./src/basis.f90')
        j=j+1
        print('Count: {0:d}, replacing p({1:d}) by d({2:d}),{3:d}<--list element'.format(j,i,i,purify[-j]))
    print(purify)




order = raw_input('Please input the maximum order of the polynomial: ')
symmetry = raw_input('Please input the permutation symmetry of the molecule: ')
train_x = raw_input('Please input the name of the data file: ')
arg = order +' '+ symmetry

print("")
print("Generating the fitting bases... (This might take time) \n")

cl('''
cd src
cd emsa
make
cp msa ../
cd ../
./msa '''+ arg +  '''
./postemsa.pl ''' + arg + '''
./derivative.pl ''' + arg
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
print('Number of atoms is : ' + str (natom))
print('Number of configurations is: '+str(nconfig)+'\n')

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
print('Given polynomial order: '+ order)
print('Given symmetry: '+ symmetry)
print ('Number of coefficients is: ' + str(ncoeff) +'\n')


ans = raw_input('Would you like to continue? y/n \n')
ans ='y'
if ans == 'n':
    print('Fitting program terminated')
    quit()

ans = raw_input('Do you want to use purified MSA? y/n \n')
ans='y'
if ans == 'y':
    #interd = [1, 8, 9, 10]
    purify_f90(natom, ncoeff) #Generate purify.f90 that can give 'purified_coeff.dat'
    cl('''cd src
            gfortran basis.f90 purify.f90''')
    cl('''cd src
    ./a.out''')  #Printed zero terms in 'purified_coeff.dat'
    purified_coeff() # Read zero terms and modify 'basis.f90'

a0='2.5'
wt = '1000000'
g = open('./src/fit.f90','w')
g.write('''program fit
use basis
implicit none

  external dgelss
  real :: rmse, wrmse, a0, dwt, vmin
  integer :: data_size
  real,allocatable::xyz(:,:,:),x(:)
  real,allocatable::v(:),b(:),p(:),wt(:)
  real,allocatable::yij(:,:), A(:,:)
  real,allocatable::coeff(:),v_out(:),s(:)
  real :: work(150000), dr(3)
  integer :: ncoeff, natm, ndis
  integer :: i, j, k, m, info, rank
  character(len=32) :: data_file
  character :: symb

  natm=''' + str(natom) + '''         ! change to the number of atoms
  ncoeff=''' + str(ncoeff)+ '''       ! change to the number of coeff. (size of c in bemsa.f90)
  data_size=''' + str(nconfig)+ '''    ! change to the number of data points in pts.dat
  a0=''' + a0 + '''
  dwt=''' + wt + '''

  ndis=natm*(natm-1)/2

  open(10,file='coeff.dat',status='unknown')
  open(11,FILE='points.eng',status='unknown')
  open(12,file="'''+train_x+'''",status='old')

  allocate(x(ndis))
  allocate(xyz(data_size,natm,3))
  allocate(v(data_size),v_out(data_size),b(data_size),coeff(ncoeff),s(ncoeff))
  allocate(yij(data_size,ndis))
  allocate(p(ncoeff))
  allocate(A(data_size,ncoeff))
  allocate(wt(data_size))

  do i=1,data_size
     read(12,*)
     read(12,*) v(i)
     do j=1,natm
        read(12,*) symb,xyz(i,j,:)
     end do
  end do
  vmin = minval(v)

  do m=1,data_size
     k = 1
     do i=1,natm-1
        do j=i+1,natm

           yij(m,k)=0.0
           dr=xyz(m,i,:)-xyz(m,j,:)
           yij(m,k)=sqrt(dot_product(dr,dr))
           yij(m,k)=yij(m,k)/0.5291772083
           yij(m,k)=exp(-yij(m,k)/a0)

           k=k+1
        end do
     end do
  end do

  do i=1,data_size
     wt(i)=dwt/(dwt+v(i)-vmin)
     x=yij(i,:)
     call bemsav(x,p)
     A(i,:)=p*wt(i)
     b(i)=v(i)*wt(i)
  end do

  call dgelss(data_size,ncoeff,1,A,data_size,b,data_size,s,1.0d-8,rank,work,150000,info)
  coeff(:)=b(1:ncoeff)

  do i=1,ncoeff
     write(10,*) coeff(i)
  end do

  rmse=0.0
  wrmse=0.0
  do i=1,data_size
     v_out(i)=emsav(yij(i,:),coeff)
     write (11,*) v(i),v_out(i),abs(v(i)-v_out(i))*219474.63
     rmse=rmse+(v(i)-v_out(i))**2
     wrmse=wrmse+(wt(i)*(v(i)-v_out(i)))**2
  end do

  rmse=sqrt(rmse/dble(data_size))
  wrmse=sqrt(wrmse/dble(data_size))
  write(*,'(A)') '3. Fitting is finished: '
  write(*,'(A,F15.10,A)') 'Overall Root-mean-square fitting error: ', rmse , ' Hartree'
  write(*,'(A,F15.10,A)') 'Weighted Root-mean-square fitting error: ', wrmse , ' Hartree'

  deallocate(xyz,x,v,b,p,yij,A,coeff,v_out,s,wt)
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
rm fit.x
mv ./src/basis.f90 ./
mv ./src/gradient.f90 ./
cp ./src/Makefile ./
cp ./src/getpot.f90 ./ '''
)

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

    a0 = ''' + a0 + '''

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

    a0 = ''' + a0 + '''  ! the same as in fitting code

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

cl('''make getpot.x
cp ./src/test.xyz ./
cp ./src/expected.out ./'''
)

print ('4. In order to run the test program, use command:')
print ('./getpot.x test.xyz \n')
print ('End of program')
