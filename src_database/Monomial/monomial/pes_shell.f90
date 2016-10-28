module pes_shell
use bemsa
use dbemsa
implicit none

  real::coeff(1:904) ! change to number of coefficients
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
    real,dimension(1:904)::p   ! change to number of p
    real,dimension(1:64)::m  ! change to number of m
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

  !==============================!
  ! Calculate the gradient using !
  ! finite difference            !
  !==============================!
  function grad(xyz) result(gd)
    real,dimension(:,:),intent(in)::xyz
    !::::::::::::::::::
    real::eps
    real,dimension(3,size(xyz,2))::xt
    real,dimension(3*size(xyz,2))::gd
    real::fa,fb
    integer::i,j,k

    eps = 1.d-3

    do i=1,size(xyz,2)
       do j=1,3
          xt=xyz; xt(j,i)=xt(j,i)-eps
          fa=f(xt)
          xt=xyz; xt(j,i)=xt(j,i)+eps
          fb=f(xt)
          k=3*(i-1)+j
          gd(k)=0.5d0*(fb-fa)/eps
       end do
    end do

  end function grad

end module pes_shell
