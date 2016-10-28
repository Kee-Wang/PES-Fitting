!=============================!
! main code comes here        !
!=============================!
program main
  use fmd_util
!  use pes, wp=>pes_wp
  implicit none
  real,parameter :: error2=-0.0184
  real           :: error1
!BEN added the use pes, wp=>pes_wp above and seed stuff below
  character(len=255)::mdint,base_name,argmnt,label,dname
  integer::i,j,ntraj,nstep,rstep
  real::mbl

  integer::stmp1,stmp2
  integer, dimension(2) :: seed
 
  call getarg(1,mdint)
  
  i=index(mdint,'.',.true.)
  base_name=mdint(1:i-1)

!  call getarg(2,argmnt); read(argmnt,*) ntraj
!  call getarg(3,argmnt); read(argmnt,*) nstep
!  call getarg(4,argmnt); read(argmnt,*) rstep
!  call getarg(5,label );

!  call getarg(2, dname)
  call getarg(2,argmnt); read(argmnt,*) ntraj
  call getarg(3,argmnt); read(argmnt,*) nstep
  call getarg(4,argmnt); read(argmnt,*) rstep
  call getarg(5,argmnt); read(argmnt,*) stmp1
  call getarg(6,argmnt); read(argmnt,*) stmp2
  call getarg(7,label);

  seed(1)=stmp1
  seed(2)=stmp2

!BEN
!  call init_pes(dname)
   
  open(16,status='unknown',file=trim(base_name)//'.gem')
  open(13,status='unknown',file=trim(base_name)//'.eng')

!BEN
  open(14,status='unknown',file=trim(base_name)//'.xyz')
  open(15,status='unknown',file=trim(base_name)//'.seed')

  write(15,*) 0,seed
  write(13,'(A)') &
       "#     TRAJ       STEP      TIME_FS      R_POTE_CM    R_TOTALE_CM     KINE_CM       R_ERROR"
  call flush(13)
  
  do i=1,ntraj
     open(11,status='old',file=trim(mdint))
     call fmd_init(11,i,seed)
     close(11)
     error1=fmd_ke*1.05
!     if (mod(i,100).eq.0.OR.i.eq.1) print*,"error1=1.05*fmd_ke hartree:", error1

     call xcom(fmd_mass,fmd_x)
     call vcom(fmd_mass,fmd_v)
!set J=0
     call diag_it(fmd_mass,fmd_x,fmd_v)
     call scale_kine(fmd_mass,fmd_v,fmd_ke)
     fmd_te = calc_kine(fmd_mass,fmd_v)

     call calc_deriv(fmd_natm,fmd_x,fmd_g)

     call record_energy(13,i,0)
     call record_geom(16,label,fmd_x,fmd_g,fmd_v,i,0)
!BEN
     call record_acc(14,fmd_x,fmd_v,i,0)
     
     do j=1,nstep
        call propagating(fmd_mass,fmd_x,fmd_v,fmd_g)
        
        call update_mbl(fmd_x,mbl)       
        if(mbl > fmd_mbl) then
           call record_energy(13,i,j)
           call record_geom(16,label,fmd_x,fmd_g,fmd_v,i,j)
          if((fmd_ve-fmd_be).gt.error1.OR.(fmd_ve-fmd_be).lt.error2) then
            call record_acc_bad(14,fmd_x,fmd_v,i,j)
          else
!BEN
            call record_acc(14,fmd_x,fmd_v,i,j)
          endif

           write(*,*) "The Molecule has been blowed up!"
           exit 
        end if
!ychan beg
!        if(mod(j,rstep)==0.OR.&
!&          (fmd_ve-fmd_be).gt.error1.OR.(fmd_ve-fmd_be).lt.error2) then
!ychan end
        if(mod(j,rstep)==0) then
           call record_energy(13,i,j)
           call record_geom(16,label,fmd_x,fmd_g,fmd_v,i,j)
           call record_acc(14,fmd_x,fmd_v,i,j)
        end if
     end do
  end do
  
  close(16)
  close(13)
  close(14)
  close(15)
  
contains
  subroutine fmd_init(file,iflag,seed)
    integer, parameter :: wp=selected_real_kind(12,300)
    integer,intent(in)::file
! BEN added iflag
    integer,intent(in)::iflag
    integer, dimension(2) :: seed
    ! ::::::::::::::::::::
    integer::i,dim
    character(len=2)::symb

! BEN
    integer, parameter :: nk=6
    real (kind=wp), dimension(0:2,0:nk-1) :: xn,xnder,gf,xn1
    real (kind=wp), dimension(nk,3) :: xx
    real (kind=wp) :: dfdis = 1.0d-3,fa,fb
    real :: fmd_ve2

    
    integer,dimension(1:8)::dandt

    if(iflag==1) then    
    ! ========================================
    ! iniate the random seed
    ! ::::::::::::::::::::
    call date_and_time(values=dandt)
    seed(1)=dandt(8)*1000000+seed(1)
    seed(2)=seed(2)*10000+dandt(7)*100+dandt(6)+dandt(8)
    call random_seed(put=seed)
    ! ========================================
    else
    call random_seed(get=seed)
    endif
    
    ! Skip the comment line
    read(file,*)          
    
    ! read in some initial parameters
    read(file,*) fmd_natm,fmd_stpz,fmd_ke,fmd_be,fmd_mbl
    
    dim=3*fmd_natm

! BEN added iflag
    if(iflag==1) then    
       allocate(fmd_x(1:dim))
       allocate(fmd_v(1:dim))      
       allocate(fmd_g(1:dim))
       allocate(fmd_symb(1:fmd_natm))
       allocate(fmd_mass(1:fmd_natm))
    endif

    ! initialize the velocity
    call random_number(fmd_v)
    fmd_v=fmd_v-0.5   !why?  The velocity is within [-0.5,0.5],otherwise
                      !there's no negative velocity
    
    ! convert the initial kinetic energy to hartree
    fmd_ke=fmd_ke/fmd_aucm      

    ! skip the energy line of the xyz file format
    read(file,*)                
    
    
    do i=1,fmd_natm
       read(file,*) symb,fmd_x(3*i-2),fmd_x(3*i-1),fmd_x(3*i)
       
       select case(symb)
       case("C")
          fmd_mass(i)=c_mass
       case("H")
          fmd_mass(i)=h_mass
       case("N")
          fmd_mass(i)=n_mass
       case("O")
          fmd_mass(i)=o_mass
       end select
       
       fmd_symb(i)=symb
    end do
    ! convert the coordinates from angstrom to bohr
    fmd_x=fmd_x/fmd_auang
    
    ! calculate the initial kinetic energy
!BEN
    fmd_te=fmd_ke

    ! calculate the initial potential energy

!  BEN
    do i=1,nk
        xx(i,1) = fmd_x(3*i-2)
        xx(i,2) = fmd_x(3*i-1)
        xx(i,3) = fmd_x(3*i)
        xn(0,i-1) = fmd_x(3*i-2)
        xn(1,i-1) = fmd_x(3*i-1)
        xn(2,i-1) = fmd_x(3*i)
    enddo

    call prepot()
    call calcpot(fmd_ve,xx)
!    fmd_ve=getpot(xn)

    call grad(xnder,dfdis,xn)

    do i=1,nk
       fmd_g(3*i-2) = xnder(0,i-1)
       fmd_g(3*i-1) = xnder(1,i-1)
       fmd_g(3*i)   = xnder(2,i-1)
    enddo
    
    fmd_ae=fmd_ke+fmd_ve
    
    fmd_time=0
    
    return
  end subroutine fmd_init



end program main
