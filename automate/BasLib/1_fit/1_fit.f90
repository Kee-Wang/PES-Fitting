PROGRAM fit
use pes, wp=>pes_wp
use px
implicit none
character (len=255) :: iargv, dname, fname
real (kind=wp), allocatable :: coef(:)
call getarg (1, dname)  ! coeff dir
call getarg (2, fname)  ! data file

pes2_dwt=(/0.02_wp,0.2_wp/)

call pes0_init (dir=dname)
call pes1_init (pes_x3y2z1_sysold)
call pes2_init (fname, pes_x3y2z1_nk, pes_x3y2z1_pot)
px_pcv(2:6) = (/ &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t( 6, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr) /)
allocate (coef(0:px_x3y2z1_nbase()-1))

call px_lsq (px_x3y2z1_base, coef)
call px_errf (px_x3y2z1_base, coef)
call cx_getcf (pes_x3y2z1_getcf, pes_x3y2z1_nki, pes_x3y2z1_sysnew, px_pcv, coef)
call pes2_reinit (pes_x3y2z1_pot)
call pes2_errf (pes_x3y2z1_nk, 1)

END PROGRAM fit
