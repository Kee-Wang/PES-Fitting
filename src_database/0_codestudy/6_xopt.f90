!=================================!
! here comes the main code        !
! to do the optimization          !
!=================================!

program main
  use opt
  implicit none

  double precision::enrg
  integer::flag,i,ii,jj,n
  character(len=32)::mole_name,flag_s,base_name

  call getarg(1,mole_name)
  call getarg(2,flag_s)

  read(flag_s,*) flag

  open(12,status='old',file=trim(mole_name))
  call opt_init(12)
  close(12)

  i=index(mole_name,'.',.true.)
  base_name=mole_name(1:i-1)

  open(11,status='unknown',file=trim(base_name)//'.log')
  call optg(opt_x,flag,11)
  close(11)

  do i=1,10
     xx(:,i)=opt_x(3*i-2:3*i)
  end do
  enrg=f(xx)

  open(13,status='unknown',file=trim(base_name)//'.opt')
  call prtxmol(opt_x,opt_symb,enrg,13)
  close(13)

end program main

  !=========================================!
  ! read in the geometry information        !
  ! of the molecule                         !
  !=========================================!
  subroutine opt_init(file)
    use opt

    integer,intent(in)::file
    ! ::::::::::::::::::::
    integer::natm,i,dim

    read(file,*) natm
    read(file,*)
    dim=3*natm

    call pes_init()

    allocate(opt_x(1:dim))
    allocate(opt_symb(1:natm))

    do i=1,natm
       read(file,*) opt_symb(i),opt_x(3*i-2),opt_x(3*i-1),opt_x(3*i)
    end do
    opt_x = opt_x / auang

    return
  end subroutine opt_init

