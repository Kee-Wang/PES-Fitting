module dmc_proc
  use pes_shell
  implicit none
  
  ! GLOBAL VARIABLES
  real::stepsize,init_q,alpha
  integer::fstate,fnormal,birth_flag

contains
  !==================================================
  ! convert mass weighted Cartesian displacement in bohr
  ! to Cartesian coor in bohr
  !==================================================
  subroutine mwc2c(pm,x,refx,rmass)
    real,dimension(:),intent(in)::pm,refx,rmass
    real,dimension(:),intent(inout)::x

    integer::i,dim

    dim=size(pm)

    do i=1,dim
       x(i)=pm(i)/rmass(ceiling(dble(i)/3.0))
       x(i)=x(i)+refx(i)
    end do

    return
  end subroutine mwc2c

  !==================================================
  ! convert 3N-6 scaled mw normal mode
  ! to 3N mw Cartesian displacement
  !==================================================
  subroutine q2c(q,B,x)
    real,dimension(:),intent(in)::q
    real,dimension(:,:),intent(in)::B
    real,dimension(:),intent(out)::x

    !q=diff+BQ
    x=matmul(B,q)

    return
  end subroutine q2c

  !==================================================
  ! get potential at one normal mode coor of 3N-6 dim
  !==================================================
  subroutine pot(v,refx,rmass,refv,q)
    real,intent(out)::v
    real,dimension(:),intent(in)::refx,rmass
    real,dimension(:,:),intent(in)::refv
    real,dimension(:),intent(in)::q
    !::::::::::::::::::::
    real,dimension(1:size(refx))::pm,x
    integer::dim,ii,jj
    real,dimension(3,5)::xyz

    call q2c(q,refv,pm)
    call mwc2c(pm,x,refx,rmass)
    do ii=1,5
       do jj=1,3
          xyz(jj,ii)=x(3*ii-3+jj)
       end do
    end do
    v=ff(xyz)

    return
  end subroutine pot


  !==================================================
  ! move all alive psips one step forward
  !==================================================
  subroutine walk(psips,dx)
    real,dimension(:,:,:),intent(inout)::psips
    real,dimension(:),intent(in),optional::dx
    !::::::::::::::::::::
    real,dimension(1:size(psips,1))::x
    real::delta
    integer::dim,i

    dim=size(psips,2)
    
    if (present(dx)) then
       do i=1,dim
          call gasdev(x)
          psips(:,i,2)=psips(:,i,1)+x*dx(ceiling(dble(i)/3.d0))
       end do
    else
       delta=sqrt(stepsize)
       do i=1,dim
          call gasdev(x)
          psips(:,i,2)=psips(:,i,1)+x*delta
       end do
    end if
    
    return
  end subroutine walk
  
  !==================================================
  ! birth-death processes
  !==================================================
  subroutine branch(refx,rmass,refv,symb,vmin,psips,psips_f,v_ref)
    real,dimension(:),intent(in)::refx,rmass
    real,dimension(:,:),intent(in)::refv
    character(len=2),dimension(:),intent(in)::symb
    real,intent(in)::vmin
    real,dimension(:,:,:),intent(inout)::psips
    integer,dimension(0:),intent(inout)::psips_f
    real,intent(inout)::v_ref
    !::::::::::::::::::::
    integer::n_rep,i,j,mark
    real::v_tot

    n_rep=psips_f(0)
    mark=n_rep
    v_tot=0.d0

    if (fstate==0) then
       do i=1,n_rep
          call gbranch(refx,rmass,refv,symb,vmin,&
               psips(:,:,2),psips_f,v_ref,v_tot,i,mark)
       end do
    else
       do i=1,n_rep
          call ebranch(refx,rmass,refv,symb,vmin,&
               psips,psips_f,v_ref,v_tot,i,mark)
       end do
    end if


    !##################################################
    ! remove all dead replicas
    !##################################################
    j=0
    do i=1,mark
       if (psips_f(i)==1) then
          j=j+1
          psips(j,:,2)=psips(i,:,2)
          psips_f(j)=psips_f(i)
       end if
    end do
    psips_f(0)=j

    psips(:,:,1)=psips(:,:,2)

    !##################################################
    ! update V_ref to V_ave_step
    !##################################################
    v_ref=v_tot/psips_f(0)+alpha*(1.d0-3.d0*dble(psips_f(0))/dble(size(psips_f)-1))

  end subroutine branch

  !==================================================
  ! birth-death criteria for excited state
  !==================================================
  subroutine ebranch(refx,rmass,refv,symb,vmin,&
       psips,psips_f,v_ref,v_tot,iw,mark)
    real,dimension(:),intent(in)::refx,rmass
    real,dimension(:,:),intent(in)::refv
    character(len=2),dimension(:),intent(in)::symb
    real,intent(in)::vmin
    real,dimension(:,:,:),intent(inout)::psips
    integer,dimension(0:),intent(inout)::psips_f
    real,intent(in)::v_ref
    real,intent(inout)::v_tot
    integer,intent(in)::iw
    integer,intent(inout)::mark
    !::::::::::::::::::::
    real::sigma,probx,dr(2),cross,prod

    if (fnormal==1) then
       cross=psips(iw,1,2)*init_q
    else
       call rpc(dr,psips(iw,:,:))
       cross=dr(2)*init_q
    end if

    if (cross<0.d0) then
       psips_f(iw)=0
    else
       call random_number(sigma)

       if (fnormal==1) then
          prod=psips(iw,1,1)*psips(iw,1,2)
       else
          prod=dr(1)*dr(2)*rmass(1)**2
       end if

       probx=exp(-2.d0*prod/stepsize)
       if (probx > sigma) then
          psips_f(iw)=0
       else
          call gbranch(refx,rmass,refv,symb,vmin,&
               psips(:,:,2),psips_f,v_ref,v_tot,iw,mark)
       end if
    end if

    return
  end subroutine ebranch

  !==================================================
  ! birth-death criteria for ground state
  !==================================================
  subroutine gbranch(refx,rmass,refv,symb,vmin,&
       psips,psips_f,v_ref,v_tot,iw,mark)
    real,dimension(:),intent(in)::refx,rmass
    real,dimension(:,:),intent(in)::refv
    character(len=2),dimension(:),intent(in)::symb
    real,intent(in)::vmin
    real,dimension(:,:),intent(inout)::psips
    integer,dimension(0:),intent(inout)::psips_f
    real,intent(in)::v_ref
    real,intent(inout)::v_tot
    integer,intent(in)::iw
    integer,intent(inout)::mark
    !::::::::::::::::::::
    integer::n_birth
    real::sigma,prob,v_psip
    integer ii,jj
    real,dimension(3,5)::xyz

    if (fnormal==0) then
       do ii=1,5
          do jj=1,3
             xyz(jj,ii)=psips(iw,3*ii-3+jj)
          end do
       end do
       v_psip=ff(xyz)
    else
       call pot(v_psip,refx,rmass,refv,psips(iw,:))
    end if

    v_psip=v_psip-vmin

    if (v_psip<-1.0E-5) then
       call record_err(refx,rmass,refv,symb,psips(iw,:),v_psip)
       psips_f(iw)=0
       return
!       stop
    end if

    prob=exp((v_ref-v_psip)*stepsize)
    call random_number(sigma)

    if ((1.0d0-prob)>sigma) then
       psips_f(iw)=0
    else
       v_tot=v_tot+v_psip
       psips_f(iw)=1
       
       if (prob>1) then
          prob=prob-1.0d0
          n_birth=int(prob)

          call random_number(sigma)
          if ((prob-n_birth)>sigma) then
             n_birth=n_birth+1
          end if
          if (n_birth>2) then
             birth_flag=birth_flag+1
             !n_birth=2      !as Schulten
          end if
          do while (n_birth>0)
             mark=mark+1
             n_birth=n_birth-1
             psips(mark,:)=psips(iw,:)
             psips_f(mark)=1
             v_tot=v_tot+v_psip
          end do
       end if
    end if

    return
  end subroutine gbranch

  !==================================================
  ! dr=Ro-h(1)-Ro-h(2)
  !==================================================
  subroutine rpc(dr,xx)
    real,dimension(:),intent(out)::dr
    real,dimension(:,:),intent(in)::xx
    !::::::::::::::::::::
    real,dimension(size(xx,1)/3,size(xx,1)/3)::rr
    real,dimension(1:3)::vect
    integer::natm,i,j,k

    natm=size(xx,1)/3

    do k=1,size(xx,2)
       vect(:)=xx(3*3-2:3*3,k)-xx(8*3-2:8*3,k)
       rr(3,8)=vect(1)**2+vect(2)**2+vect(3)**2
       rr(3,8)=sqrt(rr(3,8))

       vect(:)=xx(3*3-2:3*3,k)-xx(9*3-2:9*3,k)
       rr(3,9)=vect(1)**2+vect(2)**2+vect(3)**2
       rr(3,9)=sqrt(rr(3,9))
       dr(k)=rr(3,8)-rr(3,9)
    end do

    return
  end subroutine rpc

  !==================================================
  ! write error geometry
  !==================================================
  subroutine record_err(refx,rmass,refv,symb,errq,v)
    real,dimension(:),intent(in)::refx,rmass
    real,dimension(:,:),intent(in)::refv
    character(len=2),dimension(:),intent(in)::symb
    real,dimension(:),intent(in)::errq
    real,intent(in)::v
    !::::::::::::::::::::
    real,dimension(1:size(refx))::pm,errx
    real,parameter::au_ang=0.5291772083d0
    real,parameter::aucm=219474.6313710d0
    integer::i,natm

    if (fnormal==0) then
       errx=errq*au_ang
    else
       call q2c(errq,refv,pm)
       call mwc2c(pm,errx,refx,rmass)
       errx=errx*au_ang
    end if

    natm=size(refx)/3
    write(30,'(I8)') natm
    write(30,*)

    do i=1,natm
       write(30,'(A,2X,3F15.9)') symb(i),errx(3*i-2),errx(3*i-1),errx(3*i)
    end do

    if (fnormal==1) then
       do i=1,natm-2
          write(30,*) errq(3*i-2),errq(3*i-1),errq(3*i)
       end do
    end if

  end subroutine record_err

  !==================================================
  ! Gaussian random # generator 
  ! using Box-Muller transformation
  !==================================================
  subroutine gasdev(y)
    real,dimension(:),intent(out)::y
    !::::::::::::::::::::
    real::w,x(2)
    integer::i
  
    do i=1,size(y)
       w=2.0d0
       do while (w>=1.0)
          call random_number(x)
          x=2.0*x-1.0
          w=x(1)**2+x(2)**2
       end do
       w=sqrt(-2.0*log(w)/w)
       
       y(i)=x(1)*w
    end do

  end subroutine gasdev
  
end module dmc_proc

