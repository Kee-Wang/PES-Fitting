module pes_shell
contains
  subroutine pes_init()
    use pes
    character (len=*), parameter :: dname='./coef/'

    call pes0_init (dname)
    call pes1_init (pes_x3y2z1_sysall)
    return
  end subroutine pes_init

  function f(x)
    use pes
    integer,parameter::wp=selected_real_kind(12,300)
    real(kind=wp),dimension(3,6)::x
    real(kind=wp)::f

    f=pes_x3y2z1_pot(x)

    return
  end function f

end module
