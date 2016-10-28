!=============================!
! main code comes here        !
!=============================!
program main
  use fmd_util
  implicit none

  character(len=32)::mdint,base_name,argmnt,label
  integer::i,j,ntraj,nstep,rstep
  real::mbl
  integer :: stmp1,stmp2
  integer, dimension(2) :: seed

  
  call getarg(1,mdint)
  
  i=index(mdint,'.',.true.)
  base_name=mdint(1:i-1)

  call getarg(2,argmnt); read(argmnt,*) ntraj
  call getarg(3,argmnt); read(argmnt,*) nstep
  call getarg(4,argmnt); read(argmnt,*) rstep
  call getarg(5,argmnt); read(argmnt,*) stmp1
  call getarg(6,argmnt); read(argmnt,*) stmp2
  call getarg(7,label );

  seed(1)=stmp1
  seed(2)=stmp2

  call pes_init()
  
  open(12,status='unknown',file=trim(base_name)//'.gem')
  open(13,status='unknown',file=trim(base_name)//'.eng')
  open(14,status='unknown',file=trim(base_name)//'.xyz')
  write(13,'(A)') &
       "#     TRAJ       STEP      TIME_FS      R_POTE_CM    R_TOTALE_CM     KINE_CM       R_ERROR"
  call flush(13)
  
  do i=1,ntraj
     open(11,status='old',file=trim(mdint))
     call fmd_init(11,i,seed)
     close(11)

     call xcom(fmd_mass,fmd_x)
     call vcom(fmd_mass,fmd_v)

     call diag_it(fmd_mass,fmd_x,fmd_v)
     call scale_kine(fmd_mass,fmd_v,fmd_ke)
     fmd_te = calc_kine(fmd_mass,fmd_v)
     call f_eg(fmd_x,fmd_ve,fmd_g,.true.)

     call record_energy(13,i,0)
     call record_geom(12,label,fmd_x,fmd_g,fmd_v,i,0)
     call record_acc(14,fmd_x,fmd_v,i,0)

     do j=1,nstep
        call propagating(fmd_mass,fmd_x,fmd_v,fmd_g)
        if(mod(j,rstep)==0) then
           call record_energy(13,i,j)
           call record_geom(12,label,fmd_x,fmd_g,fmd_v,i,j)
           call record_acc(14,fmd_x,fmd_v,i,j)
        end if
        
        call update_mbl(fmd_x,mbl)       


        if(mbl > fmd_mbl) then
           call record_energy(13,i,j)
           call record_geom(12,label,fmd_x,fmd_g,fmd_v,i,j)
           call record_acc(14,fmd_x,fmd_v,i,j)

           write(*,*) "The Molecule has been blowed up!"
           write(*,*) fmd_mbl,mbl
           exit 
        end if
     end do
  end do
  
  close(12)
  close(13)
  close(14)
  
contains
  subroutine fmd_init(file,iflag,seed)
    integer,intent(in)::file
    integer,intent(in)::iflag
    integer,dimension(2) :: seed
    ! ::::::::::::::::::::
    integer::i,dim
    character(len=2)::symb
    

    if(iflag==1) then    
    ! ========================================
    ! iniate the random seed
    ! ::::::::::::::::::::
    call random_seed(put=seed)
    ! ========================================
    else
    call random_seed(get=seed)
    endif
    
    ! Skip the comment line
    read(file,*)          
    
    ! read in some initial parameters
    read(file,*) fmd_natm,fmd_stpz,fmd_ke,fmd_be,fmd_mbl
    
    if(iflag==1) then
    dim=3*fmd_natm
    
    allocate(fmd_x(1:dim))
    allocate(fmd_v(1:dim))      
    allocate(fmd_g(1:dim))
    allocate(fmd_symb(1:fmd_natm))
    allocate(fmd_mass(1:fmd_natm))
    end if

    ! initialize the velocity
    call random_number(fmd_v)
    fmd_v=fmd_v-0.5
    
    
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
       case("O")
          fmd_mass(i)=o_mass
       end select
       
       fmd_symb(i)=symb
    end do
    
    ! convert the coordinates from angstrom to bohr
    fmd_x=fmd_x/fmd_auang
    
    ! calculate the initial kinetic energy
    fmd_te=calc_kine(fmd_mass,fmd_v)

    ! calculate the initial potential energy
    !fmd_ve=f(fmd_x,.true.)
    !call calc_grad(fmd_x,fmd_g)
    call f_eg(fmd_x,fmd_ve,fmd_g,.true.)
    
    fmd_ae=fmd_ke+fmd_ve
    

    
    fmd_time=0
    
    return
  end subroutine fmd_init
end program main
