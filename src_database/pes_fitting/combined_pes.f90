module pes_shell
contains
  subroutine pes_init()
    use pes
    character (len=*), parameter :: dname='./coef/'

    call pes0_init (dname)
    call pes1_init (pes_x3y2z1_sysall)
    
    call water_init()  

  return
  end subroutine pes_init

  function f(x)
    use pes
    integer,parameter::wp=selected_real_kind(12,300)
    real(kind=wp),dimension(3,6)::xi
    real(kind=wp),allocatable::xw(:,:),xc(:,:)
    real(kind=wp)::f
   
    allocate(xw(3,3))
    allocate(xc(3,3))

    xc(1,:)=x(6,:)
    xc(2,:)=x(1,:)
    xc(3,:)=x(2,:)
    xw(:,1)=x(4,:)
    xw(:,2)=x(5,:)
    x2(:,3)=x(3,:)

    f=pes_x3y2z1_pot(x)+ co2(x)+h2O(x)

    return
  end function f

end module
