PROGRAM fit_dipole
use pes, wp=>pes_wp
use px
implicit none
character (len=255) :: iargv, dname, fname
real (kind=wp), allocatable :: coef(:), vcoef(:)
call getarg (1, dname)  ! coeff dir
call getarg (2, fname)  ! data file

!pes2_dwt=(/0.04_wp,0.4_wp/)
pes2_havedip = .true. ! for dipole

call pes0_init (dir=dname)
call pes1_init (pes_x2y2z1_sysold)
call pes1_vinit (pes_x2y2z1_sysold)
call pes2_init (fname, pes_x2y2z1_nk, pes_x2y2z1_pot, pes_x2y2z1_vfun)
px_pcv(2:5) = (/ &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t( 6, 3, 0, 1.0e+30*pes_bohr,2.0_wp*pes_bohr) /)
px_vpcv(2:5) = (/ &           ! for dipole
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t(-1, 3, 0, 1.0e+06*pes_bohr,2.0_wp*pes_bohr), &
  cx_t( 6, 3, 0, 1.0e+30*pes_bohr,2.0_wp*pes_bohr) /)

allocate (coef(0:px_x2y2z1_nbase()-1))
allocate (vcoef(0:px_x2y2z1_nvbase()-1))  ! dipole

! fit PES
call px_lsq (px_x2y2z1_base, coef)
call px_errf (px_x2y2z1_base, coef)
call cx_getcf (pes_x2y2z1_getcf, pes_x2y2z1_nki, pes_x2y2z1_sysnew, px_pcv, coef)

! fit dipole
call px_vlsq (px_x2y2z1_vbase, vcoef)
call px_errvf (px_x2y2z1_vbase, vcoef)
call cxv_getcf (pes_x2y2z1_getvcf, pes_x2y2z1_nki, pes_x2y2z1_sysnew, px_vpcv, vcoef)

call pes2_reinit (pes_x2y2z1_pot, pes_x2y2z1_vfun)
call pes2_errf (pes_x2y2z1_nk, 1)
call pes2_errvf (pes_x2y2z1_nk, 1) ! dipole

END PROGRAM fit_dipole
