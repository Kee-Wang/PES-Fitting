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
call pes1_init (pes_x7_sysold)
call pes2_init (fname, pes_x7_nk, pes_x7_pot)
px_pcv(2:2) = (/ &
   cx_t(2, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr) /)
allocate (coef(0:px_x7_nbase()-1))

call px_lsq (px_x7_base, coef)
call px_errf (px_x7_base, coef)
call cx_getcf (pes_x7_getcf, pes_x7_nki, pes_x7_sysnew, px_pcv, coef)
call pes2_reinit (pes_x7_pot)
call pes2_errf (pes_x7_nk, 1)

END PROGRAM fit
