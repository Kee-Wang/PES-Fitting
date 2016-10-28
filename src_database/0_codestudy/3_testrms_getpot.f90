program main
use pes_shell
  implicit none

!  integer,parameter::wp=selected_real_kind(12,300)
  real::auang=0.5291772083  ! 1 Bohr to Angstrom
  real::aucm=219474.63      ! 1 hatree to wavenumber    
  real::xx(3,6),pot0,diff   ! cooridnate is in 3*6 matrix
  character::symb(6)
  character(len=32)::filename
  integer::i,j,natm,ierr,n
  real::pot
  real::rms

  call getarg(1,filename)   ! The same file for fitting

  open(21,status='old',file=trim(filename)) ! open previous file
  i=index(filename,'.',.true.)    !find the number of first appearance of  . in filename.xyz
  filename=filename(1:i-1)      !For example, energy.xyz --> energy
  open(22,status='unknown',file=trim(filename)//".eng") !file.xyz --> file.eng

  call pes_init()    ! subroutine from pes_shell, used to intial

  rms=0.0
  n = 0
  do
     read(21,*,iostat=ierr) natm !put number of atom in it???but why not check?1
    ! This read has to before the iostat check, so that 'ierr' can be updated
    if (ierr < 0) exit         !if ierr<0, means reached the end of file
     n = n + 1                  !counting the number of points
     read(21,*) pot0            !read initial potential
     do i=1,natm                !
        read(21,*) symb(i),xx(:,i) !store position matrix
     end do
     xx=xx/auang !convert into atomic unit

     pot=f(xx) ! Calculate the potential, with atomic unit 
     diff = abs(pot-pot0)*aucm !calculate difference between 2 pots
     write(22,'(2F15.8,F13.2)') pot0,pot,diff

     rms=rms+(pot0-pot)**2 ! here rms is still in the unit of hatree

  end do

  rms=sqrt(rms/real(n))*219474.63 !out put is the hatree in the unit of cm
  write(*,*) rms

end program main
