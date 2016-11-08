module nma_proc
  use pes_shell
  implicit none

  ! Define the mass of different atoms
  real,parameter::c_mass = 12.0000000 !21874.66
  real,parameter::h_mass =  1.0078250  !1837.15
  real,parameter::d_mass =  2.0141018
  real,parameter::o_mass = 15.9949146  !29156.95

  ! Define constants
  real,parameter::emass=1822.88848
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.6313710
  
  ! Define globe variables
  real,dimension(:),allocatable::nma_x     ! nma geom in bohr
  real,dimension(:),allocatable::nma_mass  ! atom mass
  character(len=2),dimension(:),allocatable :: nma_symb
  real::nma_det                            ! default 5.0E-3 bohr
  
contains
  !==========================================
  ! calculate the mass weighted hessian at p        
  !==========================================
  subroutine mw_hessian(p,H)
    real,dimension(:),intent(in)::p
    real,dimension(:,:),intent(inout)::H
    real,dimension(3,size(p)/3)::xx
    real::f_ff,f_fb,f_bf,f_bb,fp
    real::rmass
    real,dimension(1:size(p))::tp
    integer::dim,i,j,ii,natm

    natm=size(p)/3

    tp=p
    do ii=1,natm
       xx(:,ii)=p(3*ii-2:3*ii)
    end do
!    call calcpot(fp,xx)
    fp=f(xx)

    dim=size(p)

    do i=1,dim-1
       do j=i+1,dim
          rmass=sqrt(nma_mass(ceiling(i/3.0)))*sqrt(nma_mass(ceiling(j/3.0)))

          call pt(tp,i,j, 1, 1)
          do ii=1,natm
             xx(:,ii)=tp(3*ii-2:3*ii)
          end do
          f_ff=f(xx)
          tp=p
!         f( 1, 1)

          call pt(tp,i,j, 1,-1)
          do ii=1,natm
             xx(:,ii)=tp(3*ii-2:3*ii)
          end do
          f_fb=f(xx)
          tp=p 
!         f( 1,-1)

          call pt(tp,i,j,-1, 1)
          do ii=1,natm
             xx(:,ii)=tp(3*ii-2:3*ii)
          end do
          f_bf=f(xx)
          tp=p 
!         f(-1, 1)

          call pt(tp,i,j,-1,-1)
          do ii=1,natm
             xx(:,ii)=tp(3*ii-2:3*ii)
          end do
          f_bb=f(xx)
          tp=p
!         f(-1,-1)

          H(i,j)=0.25*(f_ff-f_fb-f_bf+f_bb)/nma_det/rmass/nma_det

          H(j,i)=H(i,j)
         
       end do
    end do

    do i=1,dim
       rmass=nma_mass(ceiling(i/3.0))
       call pt(tp,i,i, 0, 1)
       do ii=1,natm
          xx(:,ii)=tp(3*ii-2:3*ii)
       end do
       f_ff=f(xx)
       tp=p
!      f( 1)

       call pt(tp,i,i, 0,-1)
       do ii=1,natm
          xx(:,ii)=tp(3*ii-2:3*ii)
       end do
       f_bb=f(xx)
       tp=p
!      f(-1)
       
       H(i,i)=(f_ff-2*fp+f_bb)/nma_det/nma_det/rmass

    end do

    return 
  end subroutine mw_hessian

  !==================================================!
  ! move the point p in i,j direction m and n        !
  ! steps(step length is equal to opt_det            !
  !==================================================!
  subroutine pt(p,i,j,m,n)
    real,dimension(:),intent(inout)::p        !the original point
    integer::i,j,m,n
    ! ::::::::::::::::::::
    p(i)=p(i)+m*nma_det
    p(j)=p(j)+n*nma_det

    return
  end subroutine pt  

  !==================================================
  ! diagonalize the hessian matrix and return        
  ! the eigen value and eigenvectors.         
  ! The original Hessian matrix  will be destroied
  !==================================================
  subroutine diag_hessian(H,w)
    real,dimension(:,:),intent(inout)::H
    real,dimension(:),intent(out)::w
    ! ::::::::::::::::::::
    real,dimension(:),allocatable::work
    integer::dim,lwork,info,i,j
    
    dim=size(H,1)
    lwork=dim*dim*10;
    allocate(work(1:lwork))
    
    call dsyev('v','u',dim,H,dim,w,work,lwork,info) 
    
    do i=1,dim
       w(i)=sign(sqrt(abs(w(i)))*aucm,w(i))
    end do
    
    return
  end subroutine diag_hessian

  !==================================================
  ! print vectors that can be visualized by xmakemole
  !==================================================
  subroutine prtxvec(x,q,symbs,mode,f)
    real,dimension(:),intent(in)::x
    real,dimension(:),intent(in)::q
    character(len=2),dimension(:),intent(in)::symbs
    integer,intent(in)::mode
    integer,intent(in)::f
    !:::::::::::::::::::::::::::
    real,dimension(1:size(symbs))::mass
    integer::i,natm,dim

    natm=size(symbs,1)

    do i=1,natm
       mass(i)=nma_mass(i)/emass/h_mass
    end do

    write(f,'(I2)') natm
    write(f,'(A,I8)') "Mode",mode 

    do i=1,natm
       write(f,'(A,3F13.8,3F10.5)') symbs(i),x(3*i-2:3*i),&
            q(3*i-2:3*i)/mass(i)*2.0
    end do

    return
  end subroutine prtxvec
  
end module nma_proc
