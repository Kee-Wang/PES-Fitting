!====================================================
! Normal Mode Analysis
! USAGE: $./RUN file.xyz i 
! file.xyz   geometry file
! i          deltax for hessian is i*5.0E-3 bohr
!====================================================
program nma
  use nma_proc
  implicit none
  
  real::det_scal
  integer::dim,i,j,file,natm
  character(len=32)::mole_name,det_s,base_name

  real,dimension(:),allocatable::xb         !Cartesian coor in bohr
  real,dimension(:),allocatable::freq       !Normal mode harmonic freq
  real,dimension(:,:),allocatable::eigv     !Eigenvectors

!  call prepot()

  call getarg(1,mole_name)
  i=index(mole_name,'.',.true.)
  base_name=mole_name(1:i-1)
  
  call getarg(2,det_s)
  read(det_s,*) det_scal
  nma_det=5.0E-3*det_scal
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  ! initialization
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  call get_file(file)
  open(file,status='old',file=trim(mole_name))

  call nma_init(file)
  close(file)

  dim=size(nma_x,1)
  allocate(xb(1:dim))
  allocate(freq(1:dim))
  allocate(eigv(1:dim,1:dim))
  ! allocate(qm(1:dim))
  xb=nma_x/auang

  !~~~~~~~~~~~~~~~~~~~~~~~~~
  ! normal mode analysis
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  call mw_hessian(xb,eigv)
  call diag_hessian(eigv,freq)

  !~~~~~~~~~~~~~~~~~~~~~~~~~
  ! print freq and vectors
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  call get_file(file)
  open(file,status='unknown',file=trim(base_name)//'.log')

  do i=1,dim,3
     write(file,'(3F15.7)') freq(i),freq(i+1),freq(i+2)
  end do
  write(file,*)
  write(file,*)
  do i=1,dim
     write(file,'(A,I8)') "Mode",i
     do j=1,dim,3
        write(file,'(3F15.7)') eigv(j,i),eigv(j+1,i),eigv(j+2,i)
     end do
     write(file,*)
  end do

  close(file)

  open(file,status='unknown',file=trim(base_name)//'.xyz')

  do i=1,dim
     call prtxvec(nma_x,eigv(:,i),nma_symb,i,file)
  end do

  close(file)

contains
  !=========================================
  ! read in geometry 
  !=========================================
  subroutine nma_init(f)
    integer,intent(in)::f

    integer::natm,i,dim
    
    read(f,*) natm
    read(f,*)
    
    dim=3*natm

    allocate(nma_x(1:dim))
    allocate(nma_symb(1:natm))
    allocate(nma_mass(1:natm))

    do i=1,natm
       read(f,*) nma_symb(i),nma_x(3*i-2:3*i)

       select case(nma_symb(i))
       case("C")
          nma_mass(i)=c_mass
       case("H")
          nma_mass(i)=h_mass
!       case("D")
!          nma_mass(i)=d_mass
       case("O")
          nma_mass(i)=o_mass
       end select
       
       nma_mass(i)=nma_mass(i)*emass

    end do

    call pes_init()

    return
  end subroutine nma_init
  
  !==================================================
  ! get a file handle(by zhen)
  !==================================================
  subroutine get_file(iun)
    integer,intent(out)::iun
    integer::k
    logical::b
    
    k=20
    inquire(unit=k,opened=b)

    do while(b .and. k<=100)
       k=k+1
       inquire(unit=k,opened=b)
    end do
    if(.not. b) then
       iun=k
    else
       stop "get_file: no free unit"
    end if

    return
  end subroutine get_file

end program nma
