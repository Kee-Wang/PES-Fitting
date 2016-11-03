module pes_shell
use pes

  integer,parameter::wp=selected_real_kind(12,300)

contains
  subroutine pes_init()  ! PES initialization
    character (len=*), parameter :: dname='./coefs/'  ! change if needed

    ! this will read the .dat file, not the .out
    call pes0_init (dname)
    call pes1_init (pes_x2y2z1_sysall)
    call pes1_vinit(pes_x2y2z1_sysall)

    return
  end subroutine pes_init

  function f(x)  ! for potential calculation
    real(kind=wp),dimension(:,:)::x
    real(kind=wp)::f

    f=pes_x2y2z1_pot(x)

    return
  end function f

  function dip(x)  ! for dipole calculation
    real(kind=wp),dimension(:,:)::x
    real(kind=wp)::q(size(x,2))
    real(kind=wp)::dip(3)

    q=pes_x2y2z1_vfun(x)
    dip=matmul(x,q)

    return
  end function dip

end module
