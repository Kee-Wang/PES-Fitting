! Input file (mdint) the coordinates are in angstrom.
! After the coordinates are read in the program converts them to bohr.
! Before coordinates are sent to molpro they are converted back to angstrom.
! In Input file (mdint) the max bond length (fmd_mbl) is in bohr.

module fmd_util
  use pes_shell
  implicit none

  ! some constant
  real,parameter::fmd_aukcal=627.51        ! energy: unit (kacl/mol)/au   
  real,parameter::fmd_auang=0.5291772083   ! length: unit ang/au(bohr)    
  real,parameter::fmd_aucm=219474.6313710  ! energy: unit cm^-1/au        
  real,parameter::fmd_aufs=2.4188843265e-2 ! time  : unit femi_second / au

  ! Define the mass of different atoms
  real,parameter::c_mass=21874.66
  real,parameter::h_mass=1837.150
  real,parameter::o_mass=29156.95
  
  real,parameter::fmd_det=1.0E-3

  real,dimension(:),             allocatable :: fmd_x,fmd_v,fmd_g,fmd_mass
  character(len=2),dimension(:), allocatable :: fmd_symb
  
  ! fmd_te   : kinetic energy
  ! fmd_ve   : potential energy
  ! fmd_ae   : total energy
  ! fmd_mbl  : maximum bond length
  ! fmd_ke   : initial kinetic energy
  ! fmd_be   : initial bottom energy(global minium energy)
  ! fmd_stpz : trajectory step size
  ! fmd_natm : number of atoms
  
  real    :: fmd_te,fmd_ve,fmd_time
  
  integer :: fmd_natm
  real    :: fmd_stpz
  real    :: fmd_ke,fmd_be,fmd_ae,fmd_mbl
  
contains
  !===========================!
  ! update the maximum        !
  ! bond length               !
  !===========================!
  subroutine update_mbl(x,mbl)
    real,dimension(:),intent(in)::x
    real,intent(out)::mbl
    ! ::::::::::::::::::::
    integer::i,j,n,dim
    real::tmp
    real,dimension(1:3)::dis
    
    dim=size(x)
    n=dim/3
    
    dis=0
    mbl=0
    do i=1,n-1
       do j=i+1,n
          dis=x(3*i-2:3*i)-x(3*j-2:3*j)
          tmp=sqrt(dot_product(dis,dis))
          if(tmp > mbl) mbl=tmp
       end do
    end do

    return
  end subroutine update_mbl
  !================================================================!
  ! print the molecule to a file using a coordinates vector        !
  ! the output file includes the geometry and enrgy in xyz         !
  ! format                                                         !
  !================================================================!
  subroutine prtxmol(bx,enrg,file)
    real,dimension(:),intent(in)::bx
    real,intent(in)::enrg
    integer,optional,intent(in)::file
    ! ::::::::::::::::::::
    integer::i,natm
    real,dimension(1:size(bx,1))::x
    
    natm=size(bx,1)/3

    x=bx*fmd_auang
    
    if(present(file)) then
       write(file,'(I2)') natm
       write(file,*) enrg
       
       do i=1,natm
          write(file,'(x,A,2X,3F13.8)') fmd_symb(i),x(3*i-2),x(3*i-1),x(3*i)
       end do
    else
       write(*,'(I2)') natm
       write(*,*) enrg
       
       do i=1,natm
          write(*,'(x,A,2X,3F13.8)') fmd_symb(i),x(3*i-2),x(3*i-1),x(3*i)
       end do
    end if
    
    return
  end subroutine prtxmol
  !===========================================!
  ! record the geometry information to        !
  ! a "file"                                  !
  !===========================================!
  subroutine record_geom(file,label,x,v,g,traj,step)
    integer,intent(in):: file,step,traj
    character(len=*), intent(in)::label
    real,dimension(:),intent(in)::x,v,g
    ! ::::::::::::::::::::
    real,dimension(1:3)::cord,grad,velc
    integer::i,n,dim

    dim=size(x,1)
    n=dim/3
    
    write(file,'(I2,3X,A,A,I10.10,A,I10.10)') n,trim(label),"#",traj,"#",step
    write(file,'(F20.13)') fmd_ve
    
    do i=1,n
       cord=x(3*i-2:3*i)*fmd_auang
       grad=g(3*i-2:3*i)
       velc=v(3*i-2:3*i)
       write(file,'(A,2X,3F13.8,3F13.8,3F13.8)') &
           fmd_symb(i),cord,grad,velc
    end do

    call flush(file)
    
    return
  end subroutine record_geom


  !===========================================!
  ! !BEN:  minor modification of Zhen's       !
  ! record_geom subroutine to increase        !
  ! precision of printed numbers              !
  !===========================================!
  subroutine record_acc(file,x,v,traj,step)
    integer,intent(in):: file,step,traj
    real,dimension(:),intent(in)::x,v
    ! ::::::::::::::::::::
    real,dimension(1:3)::cord,velc
    integer::i,n,dim

    dim=size(x,1)
    n=dim/3
 
    write(file,'(I2,3X,I10,3X,I10)') n,traj,step
    write(file,'(F20.13)') fmd_ve
  
    do i=1,n
       cord=x(3*i-2:3*i)*fmd_auang
       velc=v(3*i-2:3*i)
       write(file,'(A,3F15.8,3F15.8)') &
           fmd_symb(i),cord,velc
    end do

    call flush(file)

    return
  end subroutine record_acc

  !==================================!
  ! record the energy info to        !
  ! a "file"                         !
  !==================================!
  subroutine record_energy(file,traj,step)
    integer,intent(in)::file,step,traj
    ! ::::::::::::::::::::
    real::time_fs,aecm,tecm,vecm,rvcm,r_error
    
    time_fs = fmd_time        * fmd_aufs    ! Time in FS
    vecm    = fmd_ve          * fmd_aucm    ! Total Potential Energy
    tecm    = fmd_te          * fmd_aucm    ! Total Kinetic Energy
    rvcm    = (fmd_ve-fmd_be) * fmd_aucm    ! Relative Potential Energy
    aecm    = (rvcm+tecm)                   ! Relative Total Energy
    r_error = (fmd_ae-fmd_te-fmd_ve)/fmd_ae ! relative error, show the energy reservation
    
    write(file,"(I10,1X,I10,1X,5(E14.6))") &
         traj,step,time_fs,rvcm,aecm,tecm,r_error

    call flush(file)
    
    return
  end subroutine record_energy

  !===========================================!
  ! calculate the gradient info of the        !
  ! molecule configuration                    !
  !===========================================!
  subroutine calc_grad(x,g)
    real,dimension(:),intent(in) :: x
    real,dimension(:),intent(out):: g
    ! ::::::::::::::::::::
    real,dimension(1:size(x,1))::tx
    real::ff,fb
    integer::i,dim

    dim=size(x,1)

    do i=1,dim
       tx=x ; tx(i)=x(i)+fmd_det ; ff=f(tx,.true.)
       tx=x ; tx(i)=x(i)-fmd_det ; fb=f(tx,.true.)
       g(i)=0.5*(ff-fb)/fmd_det
    end do
    
    return
  end subroutine calc_grad
  !===================================================!
  ! move the molecule forward one step                !
  ! use velocity velt algorithm                       !
  ! here: te, ve,  time                               !
  ! and   x ,v , g should be updated                  !
  !===================================================!
  subroutine propagating(mass,x,v,g)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::x
    real,dimension(:),intent(inout)::v
    real,dimension(:),intent(inout)::g
    ! ::::::::::::::::::::
    real,dimension(1:size(g,1))::force,tmas
    integer::i,dim

    dim=size(x,1)
    
    do i=1,dim
       tmas(i)=mass(ceiling(i/3.0))
    end do

    force=-g
    
    x = x + v*fmd_stpz + 0.5*force/tmas*fmd_stpz**2
    
    fmd_time = fmd_time + fmd_stpz
    
    call f_eg(x,fmd_ve,g,.true.)
    
    !fmd_ve=f(x,.true.)
    
    if(fmd_ve < fmd_be) then
       call prtxmol(x,fmd_ve-fmd_be)          
    end if
    
    !call calc_grad(x,g)
    
    force = 0.5*(force-g)
    
    v = v + force/tmas*fmd_stpz
    
    fmd_te = calc_kine(mass,v)

    return
  end subroutine propagating

  !================================================!
  ! calculate the kinetic energy of current        !
  ! configuration                                  !
  !================================================!
  function calc_kine(mass,v) result(te)
    real,dimension(:)::mass
    real,dimension(:)::v
    real::te
    ! ::::::::::::::::::::
    integer::i,n,dim
    real,dimension(1:size(v,1))::sqv
    
    dim=size(v)
    
    sqv=v*v
    
    te=0
    do i=1,dim
       te=te+0.5*mass(ceiling(i/3.0))*sqv(i)
    end do
    
    return
  end function calc_kine
  
  !==================================!
  ! scale the kinetic energy         !
  !==================================!
  subroutine scale_kine(mass,v,ke)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::v
    real,intent(in)::ke
    ! ::::::::::::::::::::
    real::te,factor
    
    te=calc_kine(mass,v)
    
    factor=sqrt(ke/te)

    v=factor*v
    
    return
  end subroutine scale_kine

  !=======================================!
  ! diagnolize the interial tensor        !
  ! and set J=0                           !
  !=======================================!
  subroutine diag_it(mass,x,v)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::x
    real,dimension(:),intent(inout)::v
    ! ::::::::::::::::::::
    real,dimension(1:3,1:3)::it
    integer::i,n
    real::tmas
    real,dimension(1:3)::cord,velc,J,omega
    integer::info,lwork
    real,dimension(1:3)::egvu !Eigen Vaule
    real,dimension(1:12)::work !just temporary space
    lwork=12
    
    n=size(mass,1)

    
    it=0    
    do i=1,n
       tmas = mass(i)
       cord = x(3*i-2:3*i)
       it(1,1) = it(1,1) + tmas*( cord(2)**2 + cord(3)**2 )
       it(2,2) = it(2,2) + tmas*( cord(1)**2 + cord(3)**2 )
       it(3,3) = it(3,3) + tmas*( cord(2)**2 + cord(1)**2 )
       it(1,2) = it(1,2) - tmas*( cord(1)*cord(2) )
       it(2,3) = it(2,3) - tmas*( cord(2)*cord(3) )
       it(1,3) = it(1,3) - tmas*( cord(1)*cord(3) )
    end do

    it(2,1)=it(1,2)
    it(3,2)=it(2,3)
    it(3,1)=it(1,3)


    call dsyev('v','u',3,it,3,egvu,work,lwork,info)
   
    do i=1,n

       cord = x(3*i-2:3*i)
       x(3*i-2:3*i) = matmul(transpose(it),cord)

    end do


    !================!
    ! set J=0        !
    !================!
    J=0; velc=0; omega=0; cord=0
    
    do i=1,n
       tmas = mass(i)
       cord = x(3*i-2:3*i)
       velc = v(3*i-2:3*i)
       
       J(1) = J(1) + tmas*(cord(2)*velc(3)-cord(3)*velc(2))
       J(2) = J(2) + tmas*(cord(3)*velc(1)-cord(1)*velc(3))
       J(3) = J(3) + tmas*(cord(1)*velc(2)-cord(2)*velc(1))
    end do
    
    omega=J/egvu

    do i=1,n
       cord = x(3*i-2:3*i)
       velc = v(3*i-2:3*i)
       v(3*i-2) = velc(1) + cord(2)*omega(3) - cord(3)*omega(2)
       v(3*i-1) = velc(2) + cord(3)*omega(1) - cord(1)*omega(3)
       v(3*i)   = velc(3) + cord(1)*omega(2) - cord(2)*omega(1)
    end do

    return
  end subroutine diag_it
  
  !=============================!
  ! set the velocity to         !
  ! center of mass frame        !
  !=============================!
  subroutine vcom(mass,v)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::v
    ! ::::::::::::::::::::
    integer::i,n,dim
    real,dimension(1:3)::com_velc

    dim=size(v,1)
    n=dim/3
    
    com_velc=0
    
    do i=1,n
       com_velc=com_velc+mass(i)*v(3*i-2:3*i)
    end do
    
    com_velc=com_velc/sum(mass)

    do i=1,n
       v(3*i-2:3*i)=v(3*i-2:3*i)-com_velc
    end do
    
    return
  end subroutine vcom
  
  !=============================!
  ! set the molecule to         !
  ! center of mass frame        !
  !=============================!
  subroutine xcom(mass,x)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::x
    ! ::::::::::::::::::::
    integer::i,n,dim
    real,dimension(1:3)::com_cord !com coordinate

    dim=size(x,1)
    n=dim/3
    
    com_cord=0                  !init the com coordinates:set to 0
    
    do i=1,n
       com_cord=com_cord+mass(i)*x(3*i-2:3*i)
    end do
    
    com_cord=com_cord/sum(mass)

    do i=1,n
       x(3*i-2:3*i)=x(3*i-2:3*i)-com_cord
    end do

    return
  end subroutine xcom

end module fmd_util
