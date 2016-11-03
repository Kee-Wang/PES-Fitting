module basis
  implicit none

contains
  function emsav(x,c) result(v)
    implicit none
    real,dimension(1:0)::x
    real,dimension(0:0)::c
    real::v
    ! ::::::::::::::::::::
    real,dimension(0:0)::p

    call bemsav(x,p)
    v = dot_product(p,c)

    return
  end function emsav

  subroutine bemsav(x,p)
    implicit none
    real,dimension(1:0),intent(in)::x
    real,dimension(0:0),intent(out)::p
    ! ::::::::::::::::::::
    real,dimension(0:0)::m

    call evmono(x,m)
    call evpoly(m,p)

    return
  end subroutine bemsav

  subroutine evmono(x,m)
    implicit none
    real,dimension(1:0),intent(in)::x
    real,dimension(0:0),intent(out)::m
    !::::::::::::::::::::

    m(0) = 1.0D0

    return
  end subroutine evmono

  subroutine evpoly(m,p)
    implicit none
    real,dimension(0:0),intent(in)::m
    real,dimension(0:0),intent(out)::p
    !::::::::::::::::::::

    p(0) = m(0)

    return
  end subroutine evpoly

end module basis
