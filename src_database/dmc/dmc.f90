program dmc
  use constants
  use pes_shell
  use dmc_proc
  implicit none 

  ! PARAMETERS
  integer::natm,dim,ntraj,nwalker,nstep,cstep,nmode
  real::vmin
  real,dimension(:,:,:),allocatable::psips,psips_w
  integer,dimension(:),allocatable::psips_f,imode,vec
  real,dimension(:),allocatable::x0,q0,xmin
  real,dimension(:),allocatable::mass,deltax
  real,dimension(:,:),allocatable::vibV
  character(len=2),dimension(:),allocatable::symb

  ! VARIBLES
  character(len=32)::filename,base_name,argv
  integer::argc
  real::v_ref,v_tot,v_ave
  integer::i,j,step,traj

  call getarg(1,filename)
  call getarg(2,argv)
  read(argv,*) argc

  open(unit=40,status='old',file=filename)
  call init(40)
  close(40)

  ! LOG file for recording general info.
  i=index(filename,'.',.true.)                                                                                                 
  base_name=filename(1:i-1)
  open(unit=20,status='unknown',file=trim(base_name)//'.log')
  call log_begin(20)

  ! POT file for recording v_ref & v_ave
  open(unit=21,status='unknown',file=trim(base_name)//'.pot')
  open(unit=30,status='unknown',file="err.xyz")
  open(unit=50,status='unknown',file="last.xyz")

  !##################################################
  ! DMC START
  !##################################################  
  v_tot=0
  do traj=1,ntraj
     call init_dmc()
         
     do step=1,nstep
        if (fnormal==0) then
           call walk(psips(1:psips_f(0),:,:),deltax)
        else
           psips_w=psips(:,imode(1:nmode),:)
           call walk(psips_w(1:psips_f(0),:,:))
           psips(:,imode(1:nmode),:)=psips_w
        end if

        call branch(x0,mass,vibV,symb,vmin,psips,psips_f,v_ref)
        
        write(21,'(2I8,F10.5,F15.5,2I6)') step,psips_f(0),v_ref,v_ref*aucm

        if (step>cstep) then
           v_ave=v_ave+v_ref
        end if

        if (step > nstep-5) then
           do j=1,psips_f(0)
              write(50,'(I4)') natm
              write(50,*)
              do i=1,natm
                 write(50,'(A,3F12.5)') symb(i),psips(j,3*i-2:3*i,1)*auang
              end do
           end do
        end if

     end do

     v_ave=v_ave/(nstep-cstep)
     write(20,'(A23,I4,F15.9,F15.5)') " AVERAGE ENERGY OF TRAJ",traj,&
          v_ave,v_ave*aucm
     v_tot=v_tot+v_ave
     write(20,'(A27,F15.9,F15.5)') " CUMULATED AVERAGE ENERGY  ",&
          v_tot/dble(traj),v_tot*aucm/dble(traj)
  end do

  call log_end(20)
  call clean()

contains
  !==================================================
  ! Read in info. in INPUT file *.inp
  !==================================================
  subroutine init(file)
    integer,intent(in)::file
    !::::::::::::::::::::
    integer,dimension(1:8)::seed
    integer::fbohr,figeom,fvec,i,j,ii,jj
    character(len=2)::cdummy
    real::dr(1)
    real,dimension(3,5)::xyz
  
    !INIT RANDOM NUMBER GENERATOR
    call date_and_time(values=seed)
    seed(1)=1000*seed(6)+100*seed(8)+seed(7)+argc
    call random_seed(put=seed)
    
    !INIT MALON PES
    call pes_init()

    !READ PARAMETERS
    read(file,*)
    read(file,*) ntraj,natm,nwalker,stepsize,nstep,cstep,alpha
    read(file,*)
    read(file,*) fnormal,fstate,fbohr,figeom
    read(file,*)

    ! allocate space for walkers 
    dim=natm*3
    allocate(psips_f(0:nwalker*3))
    allocate(x0(1:dim))
    allocate(xmin(1:dim))
    allocate(symb(1:natm))
    allocate(mass(1:natm))

    if (fnormal==0) then
       allocate(deltax(1:natm))
       allocate(psips(1:nwalker*3,1:dim,2))
       figeom=0
       fvec=0
    else
       fvec=1
       allocate(q0(1:dim-6))
       allocate(imode(1:dim-6))
       allocate(psips(1:nwalker*3,1:dim-6,2))
       allocate(vibV(1:dim,1:dim-6))
       q0=0.d0
       read(file,*) nmode,imode(1:nmode)
       allocate(psips_w(1:nwalker*3,1:nmode,2))
    end if

    birth_flag=0

    ! MINIMUM CARTESIAN COOR 
    read(file,*)
    xmin=0.d0
    do i=1,natm
       read(file,*) symb(i),xmin(3*i-2:3*i)
       select case(symb(i))
       case("D")
          mass(i)= d_mass
       case("C")
          mass(i)= c_mass
       case("H")
          mass(i)= h_mass
       case("O")
          mass(i)= o_mass
       end select
       mass(i)= sqrt(mass(i)*emass)
    end do

    ! REFERENCE CARTESIAN COOR 
    read(file,*)
    x0=xmin
    do i=1,natm
       read(file,*) cdummy,x0(3*i-2:3*i)
    end do
    
    if (fbohr==0) then
       x0=x0/0.5291772083
       xmin=xmin/0.5291772083
    end if

    ! GET VMIN
    do ii=1,5
       do jj=1,3
          xyz(jj,ii)=xmin(3*ii-3+jj)
       end do
    end do 
    vmin=f(xyz)

    ! IF NORM COOR PROVIDED 
    read(file,*)
    if (figeom==1) then
       do i=1,(dim-6)/3
          read(file,*) q0(3*i-2),q0(3*i-1),q0(3*i)
       end do
    init_q=q0(1)
    end if

    ! READ IN VIB VECTORS
    read(file,*)
    if (fvec==1) then
       do i=1,dim-6
          do j=1,natm
             read(file,*) vibV(3*j-2,i),vibV(3*j-1,i),vibV(3*j,i)
          end do
       end do
    end if

    if ((fnormal==0).and.(fstate==1)) then
       call rpc(dr,reshape(xmin,(/dim,1/)))
       init_q=dr(1)
    end if

    return
  end subroutine init

  !==================================================
  ! write log file at the beginning
  !==================================================
  subroutine log_begin(file)
    integer,intent(in)::file
    !::::::::::::::::::::
    character(len=32)::fdate

    write(file,'(15X,A,A)') "DMC for ",base_name
    write(file,*)
    write(file,*) "JOB START AT:",fdate()
    write(file,'(A,I5)') " NUMBER OF RANDOM WALKERS = ",nwalker
    write(file,'(A,I5)') " NUMBER OF TRAJS = ",ntraj
    write(file,'(A,I8)') " NUMBER OF TOTAL STEPS = ",nstep
    write(file,'(A,I8)') " NUMBER OF STEPS BEFORE AVERAGING= ",cstep
    write(file,'(A,F8.2)') " INITIAL STEPSIZE = ",stepsize
    write(file,'(A,F8.3)') " INITIAL ALPHA = ",alpha
    write(file,*)

    if (fnormal==1) then
       write(file,'(A)') " USE NORMAL COORDINATES "
    else
       write(file,'(A)') " USE CARTESIAN COORDINATES"
    end if
    write(file,*) 

    if (fstate==0) then
       write(file,'(A)') " FOR GROUND STATE ENERGY"
    else
       write(file,'(A)') " FOR FIRST EXCITED STATE ENERGY"
    end if
    write(file,*)

    return
  end subroutine log_begin

  !==================================================
  ! write log file at the end
  !==================================================
  subroutine log_end(file)
    integer,intent(in)::file
    !::::::::::::::::::::
    character(len=32)::fdate
    real::etime,ta(2)

    write(file,*)
    write(file,*) "Birth replicas over 2: ",birth_flag
    write(file,*)
    write(file,'(A,F15.5,A2)') " TOTAL CPU TIME:",etime(ta)," s"
    write(file,*) "JOB ENDED AT: ",fdate()
    write(file,*) "NORMAL TERMINATION OF DMC"

  end subroutine log_end

  !==================================================
  ! initialize every DMC traj
  !==================================================
  subroutine init_dmc()
    !::::::::::::::::::::
    integer::i,ii,jj
    real,dimension(3,5)::xyz

    if (fnormal==0) then
       do i=1,natm
          deltax(i)=sqrt(stepsize)/mass(i)
       end do
    end if
    

    psips_f=1
    psips_f(0)=nwalker
    psips_f(nwalker+1:)=0

    if (fnormal==0) then
       do i=1,nwalker
          psips(i,:,1)=x0(:)
       end do
       do ii=1,5
          do jj=1,3
             xyz(jj,ii)=psips(1,3*ii-3+jj,1)
          end do
       end do
       v_ref=ff(xyz)

    else
       do i=1,nwalker
          psips(i,:,1)=q0(:)
       end do
       call pot(v_ref,x0,mass,vibV,psips(1,:,1))
    end if

    v_ave=0
    v_ref=v_ref-vmin

    do i=20,21
       !write(i,*)
       write(i,'(A,I5)') " #STARTING TRAJ:",traj
    end do

    write(21,'(2I8,F10.5,F15.5)') 0,psips_f(0),v_ref,v_ref*aucm
    !write(22,'(2I8)') 0,psips_f(0)

    return
  end subroutine init_dmc

  !==================================================
  ! Calculate Stepsize
  !==================================================
  subroutine calc_prm(i,flag)
    integer,intent(in)::i,flag
    !::::::::::::::::::::
    real::ostepsize

    ostepsize=stepsize

    if (flag==0) return
    
    if ((stepsize>1.d0).and.(mod(i,200)==0)) then
       stepsize=ostepsize*0.25d0
    end if
    
    stepsize=max(1.0d0,stepsize)
    
    if (ostepsize/=stepsize) then
       !print *,stepsize,ostepsize,deltax(1)
       deltax=deltax*sqrt(stepsize/ostepsize)
       !print *,deltax(1)
    end if 

    return
  end subroutine calc_prm

  !====================
  ! release memo space
  !====================
  subroutine clean()
    !::::::::::::::::::::
    
    deallocate(xmin,x0,mass,symb,psips,psips_f)
    if (fnormal==0) then
       deallocate(deltax)
    else 
       deallocate(q0,vibV)
    end if

    close(20)
    close(21)
    close(30)
    close(50)

    return
  end subroutine clean
  
end program dmc
