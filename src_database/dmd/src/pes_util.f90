module pes_util
  implicit none
contains
  !================================================!
  ! create the molpro input file based on a        !
  ! real vector x and the molpro input file        !
  ! templates.                                     !
  ! NOTE: the first line of molpro template        !
  ! file is the atom list.                         !
  !================================================!
  subroutine crt_mpjob(x,mpfile)
    real,dimension(:),intent(in)::x
    character(len=*),intent(in) ::mpfile
    ! ::::::::::::::::::::
    character(len=2),dimension(:),allocatable::symb
    character(len=80)::line
    integer::dim,natm,ios,tmp_iun,out_iun

    dim=size(x,1)
    natm=dim/3
    
    allocate(symb(1:natm))
    
    call get_file(tmp_iun)
    open(tmp_iun,status='old',file=trim(mpfile))
    
    ! read in the atom symbols from mpfile first
    read(tmp_iun,*) symb

    call get_file(out_iun)
    ! open the molpro input file to write
    open(out_iun,status='unknown',file="mp.inp",iostat=ios)
    
    do while(ios==0)
       read(tmp_iun,'(A)',iostat=ios) line
       
       if(ios /= 0) exit
       
       if(index(line,'<XYZ>')==0) then
          write(out_iun,'(A)') line
       else
          call pes_prtxmol(x,symb,0.0,out_iun)
       end if
    end do
    close(tmp_iun)
    close(out_iun)
    
    return
  end subroutine crt_mpjob
  !================================================================!
  ! print the molecule to a file using a coordinates vector        !
  ! the output file includes the geometry and enrgy in xyz         !
  ! format                                                         !
  !================================================================!
  subroutine pes_prtxmol(x,symbs,enrg,file)
    real,dimension(:),intent(in)::x
    character(len=2),dimension(:),intent(in)::symbs
    real,intent(in)::enrg
    integer,intent(in)::file
    ! ::::::::::::::::::::
    integer::i,natm
    
    natm=size(x,1)/3
    
    write(file,'(I2)') natm
    write(file,*) enrg
    
    do i=1,natm
       write(file,'(A,2X,3F13.8)') symbs(i),x(3*i-2),x(3*i-1),x(3*i)
    end do

    return
  end subroutine pes_prtxmol
  !===============================!
  ! get a free file handle        !
  !===============================!
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
  end subroutine get_file

end module pes_util
