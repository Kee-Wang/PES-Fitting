PROGRAM fit                     !The program is used to fit pes
use pes, wp=>pes_wp   !the function of pes and wp is used
use px                !the fucntion of px is used
implicit none         !start the main program
character (len=255) :: iargv, dname, fname 
real (kind=wp), allocatable :: coef(:)  
call getarg (1, dname)  ! coeff dir, obtain first input
call getarg (2, fname)  ! data file, obatin second input

pes2_dwt=(/0.02_wp,0.2_wp/)   !=======??????? 

call pes0_init (dir=dname)
call pes1_init (pes_x3y2z1_sysold)
call pes2_init (fname, pes_x3y2z1_nk, pes_x3y2z1_pot)
px_pcv(2:6) = (/ &      !used '&' so that they are in the same line
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &  !???????????
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t( 6, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr) /)
allocate (coef(0:px_x3y2z1_nbase()-1))          !????????????

call px_lsq (px_x3y2z1_base, coef) 
call px_errf (px_x3y2z1_base, coef) 
call cx_getcf (pes_x3y2z1_getcf, pes_x3y2z1_nki, pes_x3y2z1_sysnew, px_pcv, coef)
call pes2_reinit (pes_x3y2z1_pot)
call pes2_errf (pes_x3y2z1_nk, 1)

END PROGRAM fit
