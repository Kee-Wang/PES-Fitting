program fit
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

  natm=5         ! change to the number of atoms
  ncoeff=2304       ! change to the number of coeff. (size of c in bemsa.f90)
  data_size=44623    ! change to the number of data points in pts.dat
  a0=2.5
  dwt=1000000

  ndis=natm*(natm-1)/2

  open(10,file='coeff.dat',status='unknown')
  open(11,FILE='points.eng',status='unknown')
  open(12,file="points.dat",status='old')

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
