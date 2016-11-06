!=============================================================================!
! this code is used to do geometry optimization using Newton's Method.        !
! USAGE: xopt geom.xyz 0|1                                                    !
! Zhen Xie@Bowman Group                                                       !
! zxie3@emory.edu                                                             !
!                                                                             !
! Dec. 16 2004                                                                !
! Emerson Center , Emory University                                           !
!=============================================================================!

module opt
  use pes_shell
  implicit none
  double precision::opt_det=1.0E-4
  double precision,dimension(:),allocatable::opt_x
  character(len=2),dimension(:),allocatable::opt_symb
  double precision,dimension(3,10)::xx
  double precision,parameter::auang=0.5291772083
   
  contains
  
  !================================================================!
  ! print the molecule to a file using a coordinates vector        !
  ! the output file includes the geometry and enrgy in xyz         !
  ! format                                                         !
  !================================================================!
  subroutine prtxmol(x,symbs,enrg,file)
    double precision,dimension(:),intent(in)::x
    character(len=2),dimension(:),intent(in)::symbs
    double precision,intent(in)::enrg
    integer::file
    ! ::::::::::::::::::::
    integer::i,natm
    
    natm=size(x,1)/3
    
    write(file,'(I2)') natm
    write(file,*) enrg
    
    do i=1,natm
       write(file,'(A,3F15.9)') symbs(i),x(3*i-2:3*i)*auang
    end do

    return
  end subroutine prtxmol

  !===============================================!
  ! optimize the geometry using the newton        !
  ! method p=p-LD^{-1}L'g                         !
  !===============================================!
  subroutine optg(p,stat,file)
    double precision,dimension(:),intent(inout)::p
    integer,intent(in)::stat,file
    ! ::::::::::::::::::::
    double precision,dimension(1:size(p))::grad,w,disp,tp,tdisp
    double precision,dimension(1:size(p),1:size(p))::H,invw
    double precision::maxd,enrg,maxg,imaxg,lambda
    integer::n,i
    n=0

!    open(13,status='unknown',file='last.opt')
 
    call gradient(p,grad)       !calculate the gradient vector
    maxg=maxval(abs(grad))
    imaxg=maxg

    do i=1,10
       xx(:,i)=p(3*i-2:3*i)
    end do
    enrg=f(xx)
    call prtxmol(p,opt_symb,enrg,file)

    write(*,'(A,I3,A,F10.7,A,F15.6)') &
         "STEP:",n," MAX_GRAD:",maxg," ENERGY:", enrg
    
    do while(maxg>1.0e-7 .and. n<=30)
       call hessian(p,H)        !calculate the hessian Matrix
       call diag_hessian(H,w)   !diagnolize the hessian Matrix
       call inverse_w(w,stat,invw) !get the inverse matrix of eigenvalues
       
       call displacement(grad,H,invw,disp) !calculae the displacement 
       
       lambda=1.0               !step length in p=p-\lambda*LD^{-1}L'g
       do while(lambda>1.0e-5)
          tdisp=lambda*disp
          tp=p                  
          call updatep(tp,tdisp) !update the point
          call gradient(tp,grad) !calculate the gradient vector
          maxg=maxval(abs(grad)) !maximum gradient
          if(maxg>imaxg) then
             lambda=lambda*0.5
             cycle
          else
             disp=tdisp
             p=tp
             imaxg=maxg
             exit
          end if
       end do
       
       if(lambda<=1.0e-4) then
          write(*,*) "Optimization Terminated."
          return
       end if
       maxd=maxval(abs(disp)) !maximum displacement
       do i=1,10
          xx(:,i)=p(3*i-2:3*i)
       end do
       enrg=f(xx)
       n=n+1
       call prtxmol(p,opt_symb,enrg,file)

       write(*,'(A,I3,A,F10.7,A,F10.7,A,F12.6)') &
            "STEP:",n," MAX_GRAD:",imaxg," MAX_DISP:", maxd," ENERGY:", enrg
       
    end do
    
    write(*,*) "Optimization Completed."
    
    return
  end subroutine optg
  !=====================================================!
  ! calculate the inverse of the eigen values w         !
  ! if the w(i) is too small, set it to 0.              !
  ! in the real case, there should be 6 zeros in        !
  ! the w(i)s. The return should be a diagnonal         !
  ! matrix                                              !
  !=====================================================!
  subroutine inverse_w(w,stat,invw)
    double precision,dimension(:),intent(in)::w
    double precision,dimension(:,:),intent(out)::invw
    integer,intent(in)::stat
    ! ::::::::::::::::::::
    integer::i,dim
   
    dim=size(w,1)

    invw=0
    
    do i=1,stat
       invw(i,i)=1/w(i)
    end do

    do i=stat+7,dim
       invw(i,i)=1/w(i)
    end do
    
    return
  end subroutine inverse_w
  !=====================================================!
  ! update the current point p to the next point        !
  ! p=p+disp                                            !
  !=====================================================!
  subroutine updatep(p,disp)
    double precision,dimension(:),intent(inout)::p
    double precision,dimension(:),intent(in)::disp
    ! ::::::::::::::::::::
    p=p+disp
    return
  end subroutine updatep
  !==========================================================!
  ! calculate the displacement of point p: disp using        !
  ! the diagnolized hessian                                  !
  ! disp=-L*invw*L'g                                         !  
  !==========================================================!
  subroutine displacement(grad,norms,invw,disp)
    double precision,dimension(:),intent(in)::grad
    double precision,dimension(:,:),intent(in)::norms,invw
    double precision,dimension(:),intent(out)::disp
    ! ::::::::::::::::::::
    
    disp=matmul(transpose(norms),grad)
    disp=matmul(invw,disp)
    disp=matmul(norms,disp)

    disp=-disp

    return
  end subroutine displacement
  !==================================================!
  ! diagonalize the hessian matrix and return        !
  ! the eigen value in freq and eigenvectors         !
  ! in norms. The original Hesion matrix             !
  ! will be destroied                                !
  !==================================================!
  subroutine diag_hessian(H,w)
    double precision,dimension(:,:),intent(inout)::H
    double precision,dimension(:),intent(out)::w
    ! ::::::::::::::::::::
    double precision,dimension(:),allocatable::work
    integer::dim,lwork,info,i,j
    
    dim=size(H,1)
    lwork=dim*dim*10;
    allocate(work(1:lwork))
    
    call dsyev('v','u',dim,H,dim,w,work,lwork,info) 
    
    return
  end subroutine diag_hessian
  !===========================================!
  ! calculate the gradient at point p,        !
  ! grad is the output gradient vector        !
  !===========================================!
  subroutine gradient(p,grad)
    double precision,dimension(:),intent(in)::p
    double precision,dimension(:),intent(out)::grad
    ! ::::::::::::::::::::
    double precision,dimension(1:size(p,1))::tp
    double precision::ff,fb,fp
    integer::dim,i,j

    dim=size(p,1)
    tp=p
    do j=1,10
       xx(:,j)=p(3*j-2:3*j)
    end do
    fp=f(xx)

    do i=1,dim
       call pt(tp,i,i,0, 1)
       do j=1,10
          xx(:,j)=tp(3*j-2:3*j)
       end do
       ff=f(xx)
       tp=p !f( 1)

       call pt(tp,i,i,0,-1)
       do j=1,10
          xx(:,j)=tp(3*j-2:3*j)
       end do
       fb=f(xx)
       tp=p !f(-1)
       
       grad(i)=0.5*(ff-fb)/opt_det
    end do
    
    return
  end subroutine gradient
  !==========================================!
  ! calculate the hessian at point p,        !
  ! H is the output Hessian matrix           !
  !==========================================!
  subroutine hessian(p,H)
    double precision,dimension(:),intent(in)::p
    double precision,dimension(:,:),intent(inout)::H
    ! ::::::::::::::::::::
    double precision::f_ff,f_fb,f_bf,f_bb,fp
    double precision,dimension(1:size(p,1))::tp
    integer::dim,i,j,k

    tp=p
    do k=1,10
       xx(:,k)=p(3*k-2:3*k)
    end do
    fp=f(xx)

    dim=size(p,1)
    
    do i=1,dim-1
       do j=i+1,dim
          
          call pt(tp,i,j, 1, 1)
          do k=1,10
             xx(:,k)=tp(3*k-2:3*k)
          end do
          f_ff=f(xx)
          tp=p !f( 1, 1)
          
          call pt(tp,i,j, 1,-1)
          do k=1,10
             xx(:,k)=tp(3*k-2:3*k)
          end do
          f_fb=f(xx)
          tp=p !f( 1,-1)
         
          call pt(tp,i,j,-1, 1)
          do k=1,10
             xx(:,k)=tp(3*k-2:3*k)
          end do
          f_bf=f(xx)
          tp=p !f(-1, 1)

          call pt(tp,i,j,-1,-1)
          do k=1,10
             xx(:,k)=tp(3*k-2:3*k)
          end do
          f_bb=f(xx)
          tp=p !f(-1,-1)
          
          H(i,j)=(f_ff - f_fb - f_bf + f_bb)/(4*opt_det**2)
          
          H(j,i)=H(i,j)

       end do
    end do

    do i=1,dim
       
       call pt(tp,i,i, 0, 1)
       do k=1,10
          xx(:,k)=tp(3*k-2:3*k)
       end do
       f_ff=f(xx)
       tp=p !f( 1)

       call pt(tp,i,i, 0,-1)
       do k=1,10
          xx(:,k)=tp(3*k-2:3*k)
       end do
       f_bb=f(xx)
       tp=p !f(-1)
       
       H(i,i)=(f_ff-2*fp+f_bb)/(opt_det**2)
       
    end do
    return 
  end subroutine hessian

  !==================================================!
  ! move the point p in i,j direction m and n        !
  ! steps(step length is equal to opt_det            !
  !==================================================!
  subroutine pt(p,i,j,m,n)
    double precision,dimension(:),intent(inout)::p        !the original point
    integer::i,j,m,n
    ! ::::::::::::::::::::
    p(i)=p(i)+m*opt_det
    p(j)=p(j)+n*opt_det
    return
  end subroutine pt  
end module opt
