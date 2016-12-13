program main
use pes_shell
implicit none

  character,dimension(:),allocatable:: symb
  integer :: i, j, natm, numero,n,ierr
  real :: v
  real, dimension (:,:), allocatable :: xyz
  real, dimension (:), allocatable :: grad
  character(len=32)::fname, bname
  real,parameter::aucm = 219474.63
  real,parameter::auang= 0.5291772083d0

  call getarg(1,fname)
  i=index(fname,'.',.true.)
  bname=fname(1:i-1)

  open(12,file=trim(fname),status='old')
  open(13,file=trim(bname)//'.out',status='unknown')
  open(14,file=trim(bname)//'.cm',status='unknown')

  call pes_init()
  read (12,*) natm
  read (12,*)
  
  allocate(xyz(3,natm),grad(3*natm),symb(natm))

  do i=1,natm
     read (12,*) symb(i), xyz(:,i)
  end do
  close(12)

 open(12,file=trim(fname),status='old')

do
     read(12,*,iostat=ierr)
     if (ierr < 0) exit
     n = n + 1
     read(12,*) 
     do i=1,natm
        read(12,*) symb(i),xyz(:,i)
     end do

     v=f(xyz/auang)

     write(13,*) natm
     write(13,'(F20.10)') v
     do i=1,natm
       write (13,*) symb(i), xyz(:,i)
     end do

          write(14,*) natm
     write(14,'(F20.10)') v*aucm
     do i=1,natm
       write (14,*) symb(i), xyz(:,i)
     end do
end do


 !grad=g(xyz)

  !write(13,'(A)') "Gradient (in a.u.) from the fit:"
  !do i=1,size(grad)
  !   write(13,'(F20.10)') grad(i)
  !end do

  close (12)
  close (13)
  close (14)

end program
