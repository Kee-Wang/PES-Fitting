module pes_shell !Package the following contentents. Use 'use' to use.
contains                !contains used to inlcude functions, aka, 'subroutine'
  subroutine pes_init()  !The first subroutien
    use pes
    character (len=*), parameter :: dname='./coef/'

    call pes0_init (dname) !the function from 'pes',initialize pes?
    call pes1_init (pes_x3y2z1_sysall) !the function from 'pes'
    return
  end subroutine pes_init

  function f(x)        ! to calculate energy from pes
    use pes
    integer,parameter::wp=selected_real_kind(12,300) !the precision kind which
!gives 12 decimal precision and range between 10^-300 and 10^300
    real(kind=wp),dimension(3,6)::x !?3d coordinate*6 atoms.x is the cooridnate
    real(kind=wp)::f

    f=pes_x3y2z1_pot(x) !to calculate the energy with given position

    return
  end function f 

end module
