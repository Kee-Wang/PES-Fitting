module shell_water
  use pot_monomer_mod
!  use constants
  implicit none

contains
  !==================================================
  ! Monomer potentials using Partridge-Schwenke's
  !==================================================
  subroutine pot1b(natm,xx,pot)
    integer,intent(in)::natm
    real,dimension(3,natm),intent(in)::xx
    real,intent(inout)::pot
    !::::::::::::::::::::
    real,dimension(1:9)::x1
    real,dimension(3,3)::xr
    real,dimension(natm/3,3)::rij
    real,dimension(natm/3)::e1
 
    pot=0.d0

    x1(1:3)=xx(:,1)
    x1(4:6)=xx(:,2)
    x1(7:9)=xx(:,3)
    call bond(3,x1,xr)
    !h2o pot
    rij(1,1)=xr(1,3)  !O3-H1
    rij(1,2)=xr(2,3)  !O3-H2
    rij(1,3)=(xr(1,3)**2+xr(2,3)**2-xr(1,2)**2)*0.5/xr(1,3)/xr(2,3) !cos
    rij(1,3)=dacos(rij(1,3))  !angle H1-O3-H2

    call vibpot(rij,e1,1)
    ! potential is 0.0 when H2O is at equilibrium geometry
    pot=sum(e1)+0.000001910936

    return
  end subroutine pot1b

  !==================================================!
  ! Calculate the internuclear distances             !
  !==================================================!
  subroutine bond(natm,xx,rr)
    integer,intent(in)::natm
    real,dimension(1:natm*3),intent(in)::xx
    real,dimension(1:natm,1:natm),intent(inout)::rr
    !::::::::::::::::::::
    real,dimension(1:3)::vect
    integer::i,j
  
    do i=1,natm-1
       do j=i+1,natm
          vect(:)=xx(i*3-2:i*3)-xx(j*3-2:j*3)
          !rr(i,j)=vect(1)**2+vect(2)**2+vect(3)**2
          !rr(i,j)=sqrt(rr(i,j))
          rr(i,j)=sqrt(sum(vect*vect))
          rr(j,i)=rr(i,j)
       end do
    end do
    return
  end subroutine bond

  !==================================================!
  ! Initializing HBB water potential                 !
  !==================================================!
  subroutine water_init()

    ! monomer init
    call monomer_init()

    return
  end subroutine water_init

  !==================================================!
  ! water potential                                  !
  !   x: cartesian coordinates in bohr               !
  !  im: optional augument should only be used for   !
  !      efficient numerical gradient calculation    !
  ! * To obtain the full potential of water cluster, !
  !   one should call this function with ONLY ONE    ! 
  !   argument, simply as f(x)                       !
  !==================================================!
  function water(x) result(pot)
    real,dimension(:,:),intent(in)::x
    real::pot
    ! ::::::::::::::::::::
    real::p1
    integer::natm

    natm=size(x)/3

    p1=0.d0

    call pot1b(natm,x,p1)
    pot=p1

    return
  end function water

end module shell_water
