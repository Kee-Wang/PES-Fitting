! ========================================
! main module for pes_shell
! ========================================
module pes_shell
  use pes_util
  implicit none
  
contains
  subroutine pes_init()
    return
  end subroutine pes_init
  !=================================================!
  ! a function to calculate the sp energy of        !
  ! a molecule using a real vector x                !
  !=================================================!
  function f(y,bohr)
    real,dimension(:)::y
    logical,optional::bohr
    real::f
    ! ::::::::::::::::::::
    logical::br
    integer::iun
    real,dimension(1:size(y,1))::x !x stores the Angstrom coordinates
    
    if(present(bohr)) then
       br=bohr
    else
       br=.false.
    end if
    
    ! GET THE ANGSTROM BASED COORDINATES
    if(.not. br) then
       x=y
    else
       x=y*0.5291772083
    end if
    

    call crt_mpjob(x,"mp_sp.tmp")
    call system("/usr/local/bin/molprop_2010_1_Linux_x86_64_i8 -n 4 mp.inp")
    call system("perl ext_mpsp.pl")
    
    call get_file(iun)
    open(iun,status='old',file="mp.mpo")
    read(iun,*)
    read(iun,*) f
    close(iun)

    return
  end function f
  !====================================================!
  ! a subroutine a calculate both sp energy and        !
  ! gradient infomation                                !
  !====================================================!
  subroutine f_eg(y,enrg,grad,bohr)
    real,dimension(:),intent(in)::y
    real,intent(inout)::enrg
    real,dimension(:),intent(inout)::grad
    logical,optional,intent(in)::bohr
    ! ::::::::::::::::::::
    logical::br
    integer::iun,dim,natm,i
    character(len=2)::symb
    real,dimension(1:3)::tx
    real,dimension(1:size(y,1))::x !x stores the Angstrom coordinates

    dim=size(y,1)
    natm=dim/3
    
    if(present(bohr)) then
       br=bohr
    else
       br=.false.
    end if
    
    ! GET THE ANGSTROM BASED COORDINATES
    if(.not. br) then
       x=y
    else
       x=y*0.5291772083
    end if
    

    call crt_mpjob(x,"mp_eg.tmp")
    call system("/usr/local/bin/molprop_2010_1_Linux_x86_64_i8 -n 4 mp.inp")
    call system("perl ext_mpeg.pl")
    
    call get_file(iun)
    open(iun,status='old',file="mp.mpo")
    
    read(iun,*)
    read(iun,*) enrg

    do i=1,natm
       read(iun,*) symb,tx,grad(3*i-2:3*i)
    end do
    
    close(iun)

    return
  end subroutine f_eg
  
end module pes_shell
