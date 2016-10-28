module nma_proc
  use pes_shell
  implicit none

  ! Define the mass of different atoms
  real,parameter::c_mass = 12.0000000 !21874.66
  real,parameter::h_mass =  1.0078250  !1837.15
  real,parameter::d_mass =  2.0141018
  real,parameter::o_mass = 15.9949146  !29156.95

  ! Define constants
  real,parameter::emass=1822.88848 !mass of electron
  real,parameter::auang=0.5291772083
  real,parameter::aucm=219474.6313710
  
  ! Define globe variables
  real,dimension(:),allocatable::nma_x     ! nma geom in bohr
  real,dimension(:),allocatable::nma_mass  ! atom mass
  character(len=2),dimension(:),allocatable :: nma_symb !atom name
  real::nma_det                            ! default 5.0E-3 bohr
  
contains
  !==========================================
  ! calculate the mass weighted hessian at p        
  !==========================================
  subroutine mw_hessian(p,H)                    !mass weighted hessian
    real,dimension(:),intent(in)::p             !p is one dimensional matrix,
!intent(in) means the variable can enter but can't change values
!p is the xb, 1D cooridnate in Bohr
    real,dimension(:,:),intent(inout)::H        !H is two dimensional matrix,
!intent(inout) means the variable can enter a value and return a value 
    real,dimension(3,size(p)/3)::xx   !xx is 2D matrix, 3*(p/3)=p elements(all
!coordinates
    real::f_ff,f_fb,f_bf,f_bb,fp
    real::rmass         
    real,dimension(1:size(p))::tp !transpose of p
    integer::dim,i,j,ii,natm

    natm=size(p)/3 !number of atoms is 1/3 of all coordinates

    tp=p   !turing column p into array p
    do ii=1,natm                 !count all atoms.Used to reform the matrix
       xx(:,ii)=p(3*ii-2:3*ii)          !(x1y1z1);(x2y2z2);...(xnynzn)
    end do
!    call calcpot(fp,xx)
    fp=f(xx) !used already fitted pes f(x) to calcualte energy fp

    dim=size(p)

    do i=1,dim-1
       do j=i+1,dim
          rmass=sqrt(nma_mass(ceiling(i/3.0)))*sqrt(nma_mass(ceiling(j/3.0)))

          call pt(tp,i,j, 1, 1)
          do ii=1,natm
             xx(:,ii)=tp(3*ii-2:3*ii) !put incresed tp BACK to xx
          end do        !now xx is the cooridnate matrix increased det_step
          f_ff=f(xx) !new energy with little increased displacement
          tp=p          !turning the tp value BACK! Restore the tp matrix
!         f( 1, 1)

          call pt(tp,i,j, 1,-1) ! now first cooridnate increase, second decrease
          do ii=1,natm 
             xx(:,ii)=tp(3*ii-2:3*ii)
          end do
          f_fb=f(xx) !calculate new energy with new xx
          tp=p          !restore value of tp
!         f( 1,-1)

          call pt(tp,i,j,-1, 1) !now frist decrese, second increase
          do ii=1,natm
             xx(:,ii)=tp(3*ii-2:3*ii)
          end do
          f_bf=f(xx)
          tp=p 
!         f(-1, 1)

          call pt(tp,i,j,-1,-1) !now first decrease, second decrease
          do ii=1,natm
             xx(:,ii)=tp(3*ii-2:3*ii)
          end do
          f_bb=f(xx)
          tp=p
!         f(-1,-1)

          H(i,j)=0.25*(f_ff-f_fb-f_bf+f_bb)/nma_det/rmass/nma_det
! this is the numerical expression for hessian, off-diagonal element
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
! this is the diagonal element
! it is more convenient to separate diagnal and off-dia, but they're same
! eqution
    end do

    return 
  end subroutine mw_hessian

  !==================================================!
  ! move the point p in i,j direction m and n        !
  ! steps(step length is equal to opt_det            !
  !==================================================!
  subroutine pt(p,i,j,m,n)
    real,dimension(:),intent(inout)::p        !the original point matrix,1D
    integer::i,j,m,n
    ! ::::::::::::::::::::
    p(i)=p(i)+m*nma_det !i could be x or y or z
    p(j)=p(j)+n*nma_det !j could be x or y or z

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
    lwork=dim*dim*10; !????????why do this
    allocate(work(1:lwork))
    
    call dsyev('v','u',dim,H,dim,w,work,lwork,info)!from package
    
    do i=1,dim
       w(i)=sign(sqrt(abs(w(i)))*aucm,w(i))
    end do
    
    return
  end subroutine diag_hessian

  !==================================================
  ! print vectors that can be visualized by xmakemole
  !==================================================
  subroutine prtxvec(x,q,symbs,mode,f)
    real,dimension(:),intent(in)::x ! 1D coordinate of config
    real,dimension(:),intent(in)::q ! ith column eigv
    character(len=2),dimension(:),intent(in)::symbs
    integer,intent(in)::mode !the (modeth) ith mode
    integer,intent(in)::f !the unit of file
    !:::::::::::::::::::::::::::
    real,dimension(1:size(symbs))::mass !corresponding mass for each coordinate
    integer::i,natm,dim

    natm=size(symbs,1) !number of atoms

    do i=1,natm
       mass(i)=nma_mass(i)/emass/h_mass
    write(f,'(I2)') natm
    write(f,'(A,I8)') "Mode",mode 

    do i=1,natm
       write(f,'(A,3F13.8,3F10.5)') symbs(i),x(3*i-2:3*i),&
            q(3*i-2:3*i)/mass(i)*2.0
    end do

    return
  end subroutine prtxvec
  
end module nma_proc
