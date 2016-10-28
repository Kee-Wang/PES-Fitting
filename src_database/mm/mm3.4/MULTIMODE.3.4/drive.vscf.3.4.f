C****************************************************************
C****************************************************************
C**DRIVE
C****************************************************************
C****************************************************************
      SUBROUTINE VSCF(W)
C     PROGRAM VSCF (VERSION 3.4)
C     11-11-98
C************************************************************
C                                                           *
C     COPYRIGHT:  S. CARTER AND J.M. BOWMAN                 *
C                                                           *
C                 DEPARTMENT OF CHEMISTRY                   *
C                 EMORY UNIVERSITY                          *
C                 ATLANTA, GEORGIA 30322                    *
C                 U.S.A.                                    *
C                                                           *
C************************************************************
C
C..VERSION 1
C..1.0 VSCF, in Normal or Internal coordinates (1.10.96)
C..1.1 V-CI, including symmetry (14.10.96)
C..1.2 SCF-CI, including symmetry (1.11.96)
C..1.3 Diagonal Adiabatic Rotations (1.10.97)
C
C..VERSION 2
C..2.0 Davidson/Liu (core) algorithm for VCI (1.12.97)
C..2.1 New VCI integration scheme (1.1.98)
C..2.2 Complete (Whitehead-Handy) Rotations (1.2.98)
C..2.3 Choice of REAL*4 or REAL*8 disk grids (1.3.98)
C..2.4 Selective K-diagonal steps in rotations (1.4.98)
C..2.5 Davidson/Liu (disc) algorithm for VCI (1.6.98)
C
C..VERSION 3
C..3.0 Improved VCI algorithm for MAXBAS,MAXSUM (1.10.98)
C..3.1 Contract primitives to give HEG integration (10.10.98)
C..3.2 HEG integration for SCF virtuals (deleted)
C..3.3 Write SCF and VCI coefficients to disk (1.11.98)
C..3.4 Selective VCI symmetries ...
C..    ... and QL/Givens choice for Davidson/Liu (11.11.98)
C
C************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 TITLE(80)
      CHARACTER*40 VERSN
      LOGICAL LGIV,LINEAR,LANCZ,LANZA,LANZB,TRIAT
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/TITLE/TITLE
      COMMON/HERM/IHERM
      COMMON/PATH/ISCFCI
      COMMON/CYCLE/ICYCLE
      COMMON/DUMP/JDUMP,IDUMP
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
      COMMON/CYCLES/NCYCLE
      COMMON/MAXLAN/LANMAX
      COMMON/GIVEN/LGIV,IGIV
      COMMON/TYPE/LINEAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      COMMON/EVL/EVL,CUT
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/VCIMAX/NMAX
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/AXES/MX(3)
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(4),MAXSUM(4)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/MAXPT/MBFMAX,MBFMX1,MBFMX2,MBFMX3,MBFMX4,MBFMIN
      COMMON/DISC/IDISC
      COMMON/MODES/NMODE,NATOM
      COMMON/SIZES/KTEMP,ISIZE1,ISIZE2,ISIZE3,ISIZE4,ISIZE,JSIZE,
     1ISIZMX
      COMMON/MATRIX/NVAL,NVALR,KSTEP,KSIGN
      DATA VERSN/'MULTIMODE VERSION 3.4 (11 November 1998)'/
C**********************************************************************
100   FORMAT(80A1)
199   FORMAT(///,1X,A40,//////)
200   FORMAT(//,1X,80A1)
205   FORMAT(/,1X,'NUMBER OF MODES = ',I4,/,
     1         1X,'NUMBER OF SCF STATES = ',I4,/,
     2         1X,'NUMBER POTENTIAL TERMS = ',I4,/,
     3         1X,'PRINT LEVEL = ',I4,/,
     4         1X,'TOLERANCE FOR SCF = ',D20.12,/,
     5         1X,'COUPLING OF ',I4,' MODES',/)
210   FORMAT(/,1X,'IDISC = ',I4,/,
     1         1X,'(WRITE POTENTIAL AND CORIOLIS INFO. TO DISC (0)',/,
     2         1X,'POTENTIAL AND CORIOLIS DISCS ALREADY EXIST (1))',/)
215   FORMAT(/,1X,'SCF CALCULATION ONLY',/)
220   FORMAT(//,1X,'START SCF CALCULATION',/)
225   FORMAT(//,1X,'SCF PLUS CI CALCULATION')
230   FORMAT(/,1X,'NUMBER OF CI ENERGIES REQUIRED = ',I4,/)
235   FORMAT(//,1X,'START CI CALCULATION',/)
240   FORMAT(/,1X,'ZERO POINT ENERGY = ',F10.2,/)
245   FORMAT(//,1X,'OVERLAPS OF ORIGINAL SCF STATE FUNCTIONS',/)
250   FORMAT(/,1X,'OVERLAPS WITH STATE ',I4)
255   FORMAT(//,1X,'TEST OF SCHMIDT ORTHOGONALISATION',/)
260   FORMAT(//,1X,'FINAL VIBRATIONAL (K-DIAGONAL) CI ENERGIES',/)
265   FORMAT(//,1X,'NORMAL COORDINATE POTENTIAL',/)
270   FORMAT(//,1X,'INTERNAL COORDINATE POTENTIAL',/)
275   FORMAT(//,1X,'CI MATRIX INVOLVES ',I4,' SCF FUNCTIONS',/)
280   FORMAT(//,1X,'MAXIMUM QUANTA IN VCI BASIS',/,1X,I5,//,
     1          1X,'MAXIMUM SUM QUANTA (0: UNRESTRICTED)',/,1X,I5,/)
285   FORMAT(//,1X,'SIZE OF VIBRATIONAL CI MATRIX IS ',I5,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
290   FORMAT(//,1X,'VIBRATIONAL SYMMETRY ',I2)
295   FORMAT(//,1X,'NUMBER OF VIBRATIONAL SYMMETRIES = ',I3,/,
     1          1X,'NUMBER OF MODE SYMMETRIES = ',I3,/)
300   FORMAT(1X,'MODE SYMMETRY ',I3,'   MODES ',20I3)
305   FORMAT(//,1X,'TOTAL ANGULAR MOMENTUM J = ',I3,/)
310   FORMAT(/,1X,'J = 0 SCF CYCLE',/)
315   FORMAT(/,1X,'J = ',I3,' SCF CYCLE',/)
320   FORMAT('*************************')
325   FORMAT(50(1H*))
350   FORMAT(/,1X,'LINK TO PROGRAM NORMALS',/)
355   FORMAT(/,1X,'ILLEGAL USE OF INORM',/)
360   FORMAT(/,1X,'NUMBER OF LANCZOS CYCLES = ',I3)
365   FORMAT(/,1X,'LANCZOS TOLERANCE = ',F10.6)
370   FORMAT(/,1X,'LANCZOS (HALF) MATRIX MAXIMUM ORDER = ',I6,/)
375   FORMAT(//,1X,'ROTATIONAL SYMMETRY ',I2)
380   FORMAT(//,1X,'SIZE OF RO-VIBRATIONAL CI MATRIX IS ',I5,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
385   FORMAT(//,1X,'FINAL RO-VIBRATIONAL CI ENERGIES',/)
390   FORMAT(1X,'BLOCK ',I2,' VIBRATIONAL SYMMETRY ',I2)
395   FORMAT(//,1X,'SIZES OF CI SYMMETRY BLOCKS: ',/,1X,10I7,/)
400   FORMAT(//,1X,'SHOULD NOT OCCUR',/)
405   FORMAT(//,1X,'MINIMIZATION OF POTENTIAL',/)
410   FORMAT(//,1X,'DEGREE OF POLYNOMIAL: ',I3,/)
415   FORMAT(//,1X,'VCI BASIS DETAILS',/)
420   FORMAT(//,1X,'MAXIMUM QUANTA FOR 1- 2- 3- 4-MODE COUPLED STATES',
     1        /,1X,4I5,/)
425   FORMAT(1X,'MAXIMUM SUM QUANTA FOR 1- 2- 3- 4-MODE STATES',/,
     1       1X,4I5,/)
430   FORMAT(/,1X,'FIND MINIMUM',/)
435   FORMAT(/,1X,'REFERENCE ENERGY (J=0): ',F20.10,' CM-1',/)
440   FORMAT(1X,'MAXIMUM QUANTA FOR INDIVIDUAL MODES',/)
445   FORMAT(1X,'COUPLED INTEGRATION FOR INDIVIDUAL MODES',/)
450   FORMAT(/,1X,'INCREMENT OF K IN ROVIBRATIONAL ANALYSIS: ',I3,/)
455   FORMAT(//,1X,'NUMBER OF REQUIRED SYMMETRIES: ',I2,/,
     1          1X,'SPECIFIC SYMMETRIES: ',4I4,/)
460   FORMAT(//,1X,'CALCULATION TO PERFORM PRELIMINARY CHECKS')
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                  INITIAL SET-UP AND READ INPUT DATA
C**********************************************************************
C**********************************************************************
C**********************************************************************
      WRITE(IOUT,199)VERSN
      CALL TIMIT(2)
      LINEAR=.FALSE.
      LANCZ=.FALSE.
      LANZA=.FALSE.
      LANZB=.FALSE.
      LGIV=.FALSE.
      TRIAT=.FALSE.
      WAVENM=0.21947463D+06
      ATTOJ=0.229371D0
      BOHR=0.52917725D0
      ELMASS=1.82288853D+03
      RAD=180/DACOS(-1.D0)
      DO I=1,4
        NBAS(I)=-1
        MAXSUM(I)=0
      END DO
      ISIZE1=0
      ISIZE2=0
      ISIZE3=0
      ISIZE4=0
      ICID=0
      IRET=0
      IORDER=0
      ETA=1.D-78
      EPS=1.D-15
      TOL=ETA/EPS
      READ(INP,*)
C****************************************************************************
C**TITLE
C**         TITLE: Any suitable title...maximum 80 characters
C****************************************************************************
      READ(INP,100)TITLE
      WRITE(IOUT,200)TITLE
      WRITE(IOUT,200)
      READ(INP,*)
C****************************************************************************
C**NATOM,NSTAT,NPOT,CONV,JCOUPL,ISCFCI,IWHICH,IDISC,NROTTR,KMAX,MCHECK,INORM
C**NATOM:   Number of atoms
C**NSTAT:   Number of SCF states required
C**         NSTAT +ve: Input NSTAT specific SCF states depending on NDEF (below)
C**         NSTAT -ve: Program generates -NSTAT SCF states in increasing energy
C**NPOT:    Number of user-supplied Normal coordinate force field terms
C**         IPOT,JPOT,CPOT (below).  If INORM (below) = 1 set NPOT=0
C**         If IWHICH (below) = 1 set NPOT=NMODE (below)
C**CONV:    Convergence threshold for SCF states in cm-1
C**ICOUPL:  Number of modes coupled (1,2,3,4)
C**         If POSITIVE, use REAL*8 grid
C**         If NEGATIVE, use REAL*4 grid (half the disk storage as *8)
C**ISCFCI:  Number of CI energies required (if positive) for SCFCI or VCI
C**         (see ICI below)
C**         ISCFCI +ve: SCF + CI calculation
C**         ISCFCI=0  : SCF only
C**         ISCFCI -ve: Terminates after initial coordinate checks (see MCHECK,
C**         INORM below)
C**IWHICH:  Type of input potential
C**         IWHICH=0: Conventional Normal coordinate force field...Normal mode
C**         analysis and force field are assumed to be supplied by user
C**         IWHICH=1: Internal coordinate potential
C**IDISC:   IDISC=0: Write potential and coriolis grid values to disc
C**         IDISC=1: Assumes grid values already on disc
C**         Rigid-rotor energies written to (61),(62),(63),(64) for
C**         ICOUPL=1,2,3,4
C**         Potential written to units (71),(72),(73),(74) for ICOUPL=1,2,3,4
C**         Coriolis (vibration) written to units (81),(82),(83),(84) for
C**         ICOUPL=1,2,3,4
C**         Coriolis (rotation) written to units (91),(92),(93),(94) for
C**         ICOUPL=1,2,3,4
C**NROTTR:  Number of rotational and translational degrees of freedom to include
C**         in 'vibrational' modes (number of modes = NMODE=3*NATOM-6+NROTTR)
C**JMAX:    Total angular momentum quantum number
C**         If POSITIVE, perform exact Whitehead-Handy if ICI < 0 (below)
C**         If NEGATIVE, perform adiabatic rotation if ICI < 0 (below)
C**         For VSCF only (ISCFCI = 0) or SCFCI (ICI > 0) perform adiabatic
C**         rotation
C**MCHECK:  MCHECK=0: Use input coordinate system as principal axis system
C**         MCHECK=1: Transform equilibrium geometry to principal axis system,
C**         and calculate rotational constants.  Program continues with geometry
C**         in principal axis system (terminate if ISCFCI < 0)
C**         MCHECK=-1: Finds minimum of potential (terminate if ISCFCI < 0)
C**INORM:   INORM=0: Do not perform Normal mode analysis
C**         INORM=1: Perform Normal mode analysis with internal coordinate
C**         potential (IWHICH=1) (terminate if ISCFCI < 0)
C**         Important note:  Normal modes are ordered in increasing energy.
C****************************************************************************
      READ(INP,*)NATOM,NSTAT,NPOT,CONV,JCOUPL,ISCFCI,IWHICH,IDISC,
     1NROTTR,KMAX,MCHECK,INORM
      JMAX=IABS(KMAX)
      ICOUPL=IABS(JCOUPL)
      TRIAT=(NATOM.EQ.3)
      IF(NSTAT.LT.0)THEN
        IORDER=-NSTAT
        NSTAT=IORDER
      END IF
      NMODE=3*NATOM-6+NROTTR
      FACTOR=1.D0
      DO I=2,NMODE
        FACTOR=FACTOR*I
      END DO
      DO I=1,ICOUPL
        FACTOR=FACTOR/I
      END DO
      IDENOM=NMODE-ICOUPL
      DO I=1,IDENOM
        FACTOR=FACTOR/I
      END DO
      READ(INP,*)
C****************************************************************************
C**ICI,NMAX,CUT,EVLJ0,NVALR,KSTEP,IPRINT,MATSIZ
C**ICI:     CI basis definition
C**         ICI +ve: Number of specific SCF states in SCFCI calculation
C**         ICI -ve: -ICI quanta in each mode in VCI calculation (if NMAX.GE.0):
C**NMAX:    NMAX > 0: Maximum sum of quanta in VCI basis
C**         NMAX = 0: Unrestricted
C**         NMAX < 0: Maximum number coupled modes in basis ... maximum -4
C**CUT:     Cutoff energy (cm-1) for CI printout
C**EVLJ0:   Reference energy (cm-1) for J>0 CI calculations.  In general, this
C**         will be the zero point energy for J=0, which is printed at the end
C**         of the output for J=0.  This value will depend on the value of the
C**         user's potential at the minimum.  The user may wish to scale his
C**         potential such that the value at the minimum is zero, in which case
C**         the value printed for J=0 is the true zero point energy in cm-1
C**NVALR:   Number of rovibrational levels required
C**KSTEP:   Increment between successive K (0 to 2*J+1)
C**         For KSTEP = 0, only one K-diagonal step is carried out, for Ka = 0
C**IPRINT:  Level of output
C**         IPRINT=3: Print everything
C**         IPRINT=2: Omit rotation matrix element output if J>0
C**         IPRINT=1: Omit SCF output, and K-diagonal output if J>0
C**         IPRINT=0: Omit integration grid and normal coordinate output
C**         IPRINT.GE.0 Print ONE (largest) CI coefficient
C**         IPRINT < 0 Print -IPRINT largest CI coefficients
C**MATSIZ:  Generate VCI matrix size
C**         MATSIZ=0: Normal run
C**         MATSIZ=1: Calculate current size of VCI matrix and terminate
C****************************************************************************
      READ(INP,*)ICI,NMAX,CUT,EVLJ0,NVALR,KSTEP,JPRINT,MATSIZ
      IPRINT=IABS(JPRINT)
      IF(KSTEP.LE.0)KSTEP=0
      IF(JMAX.NE.0.AND.ICI.LT.0)WRITE(IOUT,450)KSTEP
      IF(KSTEP.EQ.0)KSTEP=100000
C**********************************
      IF(ICI.GE.0.AND.NMAX.LT.0)NMAX=0
      IF(NMAX.EQ.0)NMAX=-ICI*NMODE
      NNMAX=NMAX
CCCC  IF(NMAX.EQ.0)NNMAX=-ICI*NMODE
C**********************************
      WRITE(IOUT,305)JMAX
      WRITE(IOUT,205)NMODE,NSTAT,NPOT,IPRINT,CONV,ICOUPL
      WRITE(IOUT,210)IDISC
      IF(INORM.NE.0)WRITE(IOUT,350)
      IF(MCHECK.LT.0)WRITE(IOUT,430)
      IF(IWHICH.EQ.0)THEN
        WRITE(IOUT,265)
      ELSE
        WRITE(IOUT,270)
      END IF
      IF(INORM.NE.0.AND.IWHICH.EQ.0)THEN
        WRITE(IOUT,355)
        STOP 'ILLEGAL USE OF INORM'
      END IF
C**FOR VCI, RESET ICI TO NUMBER OF FUNCTIONS
      IF(ICI.LT.0)ICI=ICI-1
C**RESET ICI IF TOO BIG (SCF)
      IF(ICI.GT.0.AND.ICI.GT.IABS(NSTAT))ICI=IABS(NSTAT)
C**RESET ISCFCI IF NOT NORMALS RUN OR MINIMUM CHECK
      IF(ISCFCI.LT.0.AND.INORM.EQ.0.AND.MCHECK.EQ.0)ISCFCI=0
C**RESET ICI IF SCF ONLY
      IF(ISCFCI.EQ.0)THEN
        ICI=0
        WRITE(IOUT,215)
      ELSE
C**RESET ISCFCI IF TOO BIG
        IF(ICI.GT.0.AND.ISCFCI.GT.NSTAT)ISCFCI=NSTAT
        IF(ICI.GT.0.AND.ISCFCI.GT.ICI)ISCFCI=ICI
        IF(ISCFCI.LT.0)THEN
          WRITE(IOUT,460)
        ELSE
          WRITE(IOUT,225)
          WRITE(IOUT,230)ISCFCI
        END IF
      END IF
      WRITE(IOUT,435)EVLJ0
C**SELECTIVE BASIS SETS (INITIALLY VCI, BUT SCFCI NEEDS NSTAT MODS)
C**J IS NUMBER OF COUPLED MODES
      J=-NMAX
C**MAXIMUM ALLOWED IS 4 (AT THE MOMENT)
      IF(J.GT.4)J=4
      IF(J.LT.0)J=0
      IF(J.EQ.0.AND.ICI.LT.0)THEN
        JCI=-ICI-1
        WRITE(IOUT,415)
        WRITE(IOUT,280)JCI,NMAX
      END IF
      MAXJ=J
      READ(INP,*)
      IF(J.GT.0.AND.ISCFCI.GT.0.AND.ICI.LT.0)THEN
        WRITE(IOUT,415)
C****************************************************************************
C**OMIT FOLLOWING INPUT IF SCF ONLY (ISCFCI.LE.0)
C**OMIT FOLLOWING INPUT IF SCFCI (ICI > 0)
C**OMIT FOLLOWING INPUT IF UNRESTRICTED OR RESTRICTED VCI (NMAX > 0)
C**
C**MAXSUM(I),  I = 1,-NMAX
C**MAXSUM:  Maximum sum quanta in J-mode basis ... restricted only
C**         Input single record of -NMAX integers (maximum 4)
C**
C**      MAXSUM(1) is maximum sum quanta in one-mode basis sets
C**      MAXSUM(2) is maximum sum quanta in two-mode basis sets
C**      MAXSUM(3) is maximum sum quanta in three-mode basis sets
C**      MAXSUM(4) is maximum sum quanta in four-mode basis sets
C****************************************************************************
        READ(INP,*)(MAXSUM(I),I=1,J)
      END IF
      READ(INP,*)
      IF(J.GT.0.AND.ISCFCI.GT.0.AND.ICI.LT.0)THEN
C****************************************************************************
C**OMIT FOLLOWING INPUT IF SCF ONLY (ISCFCI.LE.0)
C**OMIT FOLLOWING INPUT IF SCFCI (ICI > 0)
C**OMIT FOLLOWING INPUT IF UNRESTRICTED OR RESTRICTED VCI (NMAX > 0)
C**
C**MAXBAS(MODE,-NMAX),  MODE=1,NMODE
C**MAXBAS:  Defines selective quanta in VCI basis
C**         Input -NMAX records (maximum 4) of NMODE integers
C**
C**      MAXBAS(MODE,1) is maximum quanta for mode MODE in one-mode basis sets
C**      MAXBAS(MODE,2) is maximum quanta for mode MODE in two-mode basis sets
C**      MAXBAS(MODE,3) is maximum quanta for mode MODE in three-mode basis sets
C**      MAXBAS(MODE,4) is maximum quanta for mode MODE in four-mode basis sets
C****************************************************************************
        CALL MEMO(1,LMXBAS,NMODE*J,0,0,0,0,0,0,0,0)
        IF(ICI.LT.0)WRITE(IOUT,440)
        CALL MAXIN(W(LMXBAS),NMODE,J,ICI)
        IF(ICI.LT.0)THEN
C**NBAS:    Maximum no. quanta in J-Mode basis
          WRITE(IOUT,420)(NBAS(I),I=1,J)
          WRITE(IOUT,425)(MAXSUM(I),I=1,J)
        END IF
      END IF
C**KEEP COPY
      IF(NMAX.LT.0.AND.ICI.LT.0)ICI=NMAX
      READ(INP,*)
C****************************************************************************
C**NCYCLE,TOLLAN,LANMAX,IGIV
C**NCYCLE:  Number of Lanczos cycles (Lanczos ignored if NCYCLE=0)
C**         +ve: use LANCZA (less memory usage - more disk)
C**         -ve: use LANCZB (more memory usage - less disk)
C**TOLLAN:  Tolerance for Lanczos eigenvalues
C**LANMAX:  Maximum order of half-matrix used to build Lanczos matrix
C**LGIV:    0: Use QL algorithm in Lanczos
C**         1: Use GIVENS in Lanczos
C****************************************************************************
      READ(INP,*)NCYCLE,TOLLAN,LANMAX,IGIV
      LANZA=(NCYCLE.GT.0)
      LANZB=(NCYCLE.LT.0)
      NCYCLE=IABS(NCYCLE)
      IF(ICI.LT.0)WRITE(IOUT,360)NCYCLE
      LANCZ=(NCYCLE.GT.0)
      IF(LANCZ)THEN
        IF(ICI.LT.0)WRITE(IOUT,365)TOLLAN
        IF(ICI.LT.0)WRITE(IOUT,370)LANMAX
      END IF
      READ(INP,*)
C****************************************************************************
C**NVSYM NWSYM
C**NVSYM:   Number of vibrational symmetry species (1 or 2 or 4)
C**NWSYM:   Number of mode symmetry species        (1 or 2 or 3 or 4)
C**         For no symmetry NVSYM=NWSYM=1 (C1)...can be used for molecules of
C**         any symmetry
C**         For four-atom or larger Cs molecules NVSYM=NWSYM=2 (A',A")
C**         For triatomic C2v molecules NVSYM=NWSYM=2 (A1,B2)
C**         For four-atom C2v molecules NVSYM=4 (A1,B2,B1,A2),
C**         NWSYM=3 (A1,B2,B1)
C**         For larger molecules, where reduction of symmetry to C2v is possible
C**         NVSYM=4 (A1,B2,B1,A2), NWSYM=4 (A1,B2,B1,A2)
C****************************************************************************
      READ(INP,*)NVSYM,NWSYM
      WRITE(IOUT,295)NVSYM,NWSYM
      READ(INP,*)
C****************************************************************************
C**NSYM(NWSYM) ISYM(NWSYM,NSYM(NWSYM))    (order of NWSYM:  A1,B2,B1,A2)
C**         Input NWSYM records: NSYM  ISYM(I), I=1,NSYM
C**
C**         For triatomic C2v molecules, NWSYM=1 corresponds to A1,
C**         NWSYM=2 corresponds to B2
C**         For four-atom C2V molecules, NWSYM=1 corresponds to A1,
C**         NWSYM=2 corresponds to B2, NWSYM=3 corresponds to B1
C**         etc.
C**         The user must know the symmetry of the Normal modes which are
C**         ordered either in increasing energy (see INORM=1),
C**         or according to the user's convention (see INORM=0)
C**
C**NSYM:    Number of modes of current (NWSYM)symmetry
C**ISYM:    NSYM specific modes of current (NWSYM) symmetry
C**
C**         For example, in the case of H2O, there are 2 mode vibrations
C**         of A1 symmetry (modes 1 and 2, say) corresponding to NWSYM=1.
C**         There is a single mode vibration of B2 symmetry (mode 3, say)
C**         corresponding to NWSYM=2.
C**         Hence for NWSYM=1 (totally symmetric), the input required is:
C**         NSYM=2,  ISYM(I)=1 2
C**         For I=2, the input required is:
C**         NSYM=1,  ISYM(I)=3
C****************************************************************************
      DO I=1,NWSYM
        READ(INP,*)NSYM(I),(ISYM(I,J),J=1,NSYM(I))
        DO J=1,NSYM(I)
          IF(ISYM(I,J).GT.NMODE.OR.ISYM(I,J).LE.0)STOP 'ERROR IN ISYM'
        END DO
        WRITE(IOUT,300)I,(ISYM(I,J),J=1,NSYM(I))
      END DO
      READ(INP,*)
C****************************************************************************
C**MVSYM    MWSYM(I),  I=1,MVSYM     (order of MWSYM:  A1,B2,B1,A2)
C**         Input a single record of MVSYM+1 integers:
C**MVSYM:   Single integer corresponding to the number of actual vibrational
C**         (JMAX=0) or rovibrational (JMAX>0) symmetry species required
C**MWSYM:   MVSYM integers corresponding to the Specific vibrational (JMAX=0)
C**         or rovibrational (JMAX>0) symmetry species required
C**
C**         For J = 0 get vibrational energies of specified symmetry.
C**         For J > 0 all vibrational symmetries are included in K-blocks, but
C**         final rovibrational calculations are done for symmetries specified.
C**
C**Note:    For a triatomic molecule, the number of possible rovibrational
C**         species is twice the number of possible vibrational species
C********************************************************************
C**         (Currently not yet implemented for J>0)
C********************************************************************
C****************************************************************************
      READ(INP,*)MVSYM,(MWSYM(I),I=1,MVSYM)
      WRITE(IOUT,455)MVSYM,(MWSYM(I),I=1,MVSYM)
      IF(JMAX.EQ.0.AND.MVSYM.GT.NVSYM)MVSYM=NVSYM
      IF(JMAX.GT.0.AND.MVSYM.GT.NVSYM.AND..NOT.TRIAT)MVSYM=NVSYM
      IF(JMAX.GT.0.AND.MVSYM.GT.2*NVSYM.AND.TRIAT)MVSYM=2*NVSYM
      DO I=1,MVSYM
        IF(JMAX.EQ.0.AND.MWSYM(I).GT.NVSYM.OR.MWSYM(I).LT.0)MWSYM(I)=0
        IF(JMAX.GT.0.AND..NOT.TRIAT.AND.(MWSYM(I).GT.NVSYM.OR.
     1  MWSYM(I).LT.0))MWSYM(I)=0
        IF(JMAX.GT.0.AND.TRIAT.AND.(MWSYM(I).GT.2*NVSYM.OR.
     1  MWSYM(I).LT.0))MWSYM(I)=0
      END DO
      JDUMP=0
      READ(INP,*)
      IF(ICI.LT.0.AND.ISCFCI.GT.0)THEN
C****************************************************************************
C**IDUMP:  Number of VCI functions written to disk (input only if ICI < 0)
C**         IDUMP>0: IDUMP consecutive functions for each symmetry
C**         IDUMP<0: NDUMP(IABS(IDUMP)) specific functions for each symmetry
C********************************************************************
C**         (Currently internal use only)
C********************************************************************
C****************************************************************************
        READ(INP,*)JDUMP
      END IF
      IDUMP=IABS(JDUMP)
C**CAN'T BE GREATER THAN MAXIMUM ASKED FOR
      IF(KMAX.LE.0)THEN
        IF(IDUMP.GT.ISCFCI)IDUMP=ISCFCI
      ELSE
        IF(IDUMP.GT.NVALR)IDUMP=NVALR
      END IF
      IF(JDUMP.LT.0)THEN
        READ(INP,*)
C****************************************************************************
C**NDUMP(I,NVSYM),  I=1,-IDUMP
C**        Input only if IDUMP<0
C**        Input NVSYM records of -IDUMP integers
C**NDUMP:  -IDUMP specific functions for symmetry species NVSYM to dump to disc
C********************************************************************
C**        (Currently internal use only........set IDUMP=0)
C********************************************************************
C****************************************************************************

        CALL MEMO(1,LNDUMP,MVSYM*IDUMP,0,0,0,0,0,0,0,0)
C**ANY OUT-OF-RANGE FUNCTION WILL BE GIVEN VALUE OF ZERO
        CALL INDUMP(W(LNDUMP),MVSYM,IDUMP,ISCFCI)
      END IF
C**READ REMAINING STANDARD INPUT
      KNBF=NMODE
      KMBF=NMODE
      KXM=NATOM
      KX0=NATOM*3
      KXL=NATOM*NMODE*3
      CALL MEMO(5,LNBF,KNBF,LMBF,KMBF,LXM,KXM,LX0,KX0,LXL,KXL)
      CALL MEMO(3,LYL,KXL,LY0,KX0,LNVF,KNBF,0,0,0,0)
C**MAXIMUM OF 6 MODES COUPLED!!
C**RESET NPOT TO NMODE IF INTERNAL FORCE FIELD
      IF(IWHICH.NE.0)NPOT=NMODE
      KIPOT=NPOT*6
      KJPOT=NPOT*6
      KCPOT=NPOT
      KOMEGA=NMODE
      KISTAT=NSTAT*NMODE
      CALL MEMO(5,LIPOT,KIPOT,LJPOT,KJPOT,LCPOT,KCPOT,LOMEGA,KOMEGA,
     1LISTAT,KISTAT)
      CALL MEMO(4,LJSTAT,KISTAT,LESTAT,NSTAT,LKSTAT,NMODE,LWSTAT,NSTAT,
     10,0)
      KXX=NATOM*3
      KRR=NATOM*NATOM
      CALL MEMO(2,LXX,KXX,LRR,KRR,0,0,0,0,0,0)
C**RESET ICI IF TOO BIG (VIRTUALS)
      CALL INPUT(NPOT,W(LIPOT),W(LJPOT),W(LCPOT),NMODE,W(LOMEGA),W(LXX)
     1,W(LNBF),W(LMBF),NSTAT,W(LISTAT),NATOM,W(LXM),W(LX0),W(LXL),
     2W(LRR),ICI,W(LMXBAS),MAXJ,
     3W(LJSTAT),W(LKSTAT),W(LESTAT),W(LWSTAT),W(LNVF),INORM,MCHECK,0)
      IF(ICI.GT.0)WRITE(IOUT,275)ICI
C*******************************************************************BUG
      MBFMIN=100000
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(IK2.LT.MBFMIN)MBFMIN=IK2
      END DO
C*******************************************************************BUG
      CALL FLUSH(IOUT)
      IF(IWHICH.NE.0)THEN
        READ(INP,*)
C****************************************************************************
C**USER-DEFINED POTENTIAL
C**At this point, control passes to the user to input any data he requires in
C**the routine USERIN (see Manual)
C****************************************************************************
        WRITE(IOUT,*)
        WRITE(IOUT,*)'***********************START USER PARAMETERS'
        CALL USERIN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'***********************END USER PARAMETERS'
      END IF
      IF(INORM.NE.0)THEN
        CALL GETNOR(NATOM,NMODE,W(LXM),W(LX0),W(LOMEGA),
     1  W(LXL),W(LXX),W(LRR),W(LIPOT),W(LJPOT),W(LCPOT),NPOT)
        IF(ISCFCI.LT.0.AND.MCHECK.EQ.0)STOP 'NORMALS RUN ONLY'
      END IF
C**GET MEMORY FOR INERTIA TENSOR CALCULATIONS
      KAB=NMODE*3
      KB=NMODE*NMODE
      KAA=NMODE*3*3
      KBB=NMODE
      KQQ=NMODE
      CALL MEMO(5,LAB,KAB,LB,KB,LAA,KAA,LBB,KBB,LQQ,KQQ)
C**GET HERMITE INTEGRATION POINTS AND WEIGHTS AND H.O. FUNCTIONS
      KH=0
      KXQ=0
      KXW=0
      KXV=0
C*****************************************************HEG
C*****************************************************HEG
      KTEMP=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        KH=KH+I1
        KXQ=KXQ+I2
        KXW=KXW+I2
        KXV=KXV+I2
        CALL INTARR(W(LNBF),W(LNBF),K,I1,I2,I3)
        IF(I1.GT.KTEMP)KTEMP=I1
      END DO
      CALL MEMO(5,LH,KH,LXQ,KXQ,LXW,KXW,LXV,KXV,LQM,NMODE)
      CALL MEMO(2,LTEMP,KTEMP,LSUP4,KTEMP,0,0,0,0,0,0)
      ICYCLE=0
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)THEN
        REWIND 60
        WRITE(60)NMODE
      END IF
      K2=0
      K3=0
      KX2=0
      KX3=0
      DO K=1,NMODE
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        KXK=I2
C       CALL MEMO(1,LXK,KXK*KXK,0,0,0,0,0,0,0,0)
        CALL MEMO(3,LXK,KXK*KXK,LEVAL,KXK,LYK,KXK*KXK,0,0,0,0)
        CALL HERMPT(I3,I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),W(LXV+K3),
     1  W(LOMEGA+K1),W(LXK),K,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),
     2  W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT))
        CALL DUMPT1(W(LXQ+K3),I2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     1  W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,
     2  W(LXW+K3),W(LXW+K3),0)
C**CONTRACT CURRENT BASIS
        CALL CONTRA(I3,I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),W(LOMEGA+K1),
     1  W(LYK),K,W(LNVF),NMODE,W(LEVAL))
C**FORM HEG POINTS
        CALL HEG0(I3,I2,K,NVB,MVB,W(LNVF),NMODE)
C**USE NVB CONTRACTED PRIMITIVES AT MVB 'HEG' POINTS
        CALL HEG1(I3,I2,W(LH+K2),W(LXQ+K3),W(LXW+K3),
     1  W(LOMEGA+K1),W(LXK),K,W(LTEMP),W(LH+KX2),NVB,MVB,W(LXV+KX3),
     2  W(LXQ+KX3),W(LMBF),W(LNBF),NMODE,W(LYK),W(LSUP4),NATOM,W(LQQ),
     3  W(LRR),W(LXX),W(LX0),
     4  W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT))
C**NBF, MBF MODIFIED
        K2=K2+I1
        K3=K3+I2
        CALL INTARR(W(LNBF),W(LMBF),K,IX1,IX2,IX3)
        KX2=KX2+IX1
        KX3=KX3+IX2
C       CALL MEMO(-1,LXK,KXK*KXK,0,0,0,0,0,0,0,0)
        CALL MEMO(-3,LXK,KXK*KXK,LEVAL,KXK,LYK,KXK*KXK,0,0,0,0)
      END DO
      CALL MEMO(-2,LTEMP,KTEMP,LSUP4,KTEMP,0,0,0,0,0,0)
C*****************************************************HEG
C*****************************************************HEG
C**FIND MINIMUM OF POTENTIAL
      IF(MCHECK.LT.0)THEN
        WRITE(IOUT,405)
C**START WITH QUADRATIC FIT
        N=3
1111    CONTINUE
        IF(N.LE.MBFMIN)WRITE(IOUT,410)N
        INDEX=1
C**4*N IS BIG ENEOUGH FOR E04FBF FOR N=3
        CALL MEMO(4,LXK,N*N,LEVAL,N,LSUP4,N,LWRK,4*N,0,0)
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
          CALL GETQ(W(LXK),W(LEVAL),W(LSUP4),W(LWRK),N,W(LXQ+K3),
     1    W(LXV+K3),I2,W(LQM),NMODE,K,INDEX)
          K3=K3+I2
        END DO
        CALL MEMO(-4,LXK,N*N,LEVAL,N,LSUP4,N,LWRK,4*N,0,0)
        IF(INDEX.EQ.0)THEN
          CALL GETMIN(W(LQM),NMODE,NATOM,W(LXX),W(LX0),W(LXM),W(LXL))
          K3=0
          DO K=1,NMODE
            CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
            CALL GETV(I2,W(LXQ+K3),W(LXV+K3),K,NMODE,NATOM,W(LQQ),
     1      W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     2      NPOT,W(LIPOT),W(LJPOT),W(LCPOT))
            K3=K3+I2
          END DO
          N=N+2
          GO TO 1111
        END IF
        CALL CHECKM(W(LXM),W(LX0),W(LXX),W(LRR),NATOM)
      END IF
      CALL PAXES(W(LXM),W(LX0),NATOM)
C**FIND MINIMUM OF POTENTIAL
      IF(MCHECK.LT.0)THEN
        WRITE(IOUT,350)
        CALL GETNOR(NATOM,NMODE,W(LXM),W(LX0),W(LOMEGA),
     1  W(LXL),W(LXX),W(LRR),W(LIPOT),W(LJPOT),W(LCPOT),NPOT)
        IF(ISCFCI.LT.0)STOP 'MINIMUM POTENTIAL RUN ONLY'
      END IF
      CALL INPUT(NPOT,W(LIPOT),W(LJPOT),W(LCPOT),NMODE,W(LOMEGA),W(LXX)
     1,W(LNBF),W(LMBF),NSTAT,W(LISTAT),NATOM,W(LXM),W(LX0),W(LXL),
     2W(LRR),ICI,W(LMXBAS),MAXJ,
     3W(LJSTAT),W(LKSTAT),W(LESTAT),W(LWSTAT),W(LNVF),INORM,MCHECK,1)
C**CALCULATE ZETA CONSTANTS AND EQUILIBRIUM BOND LENGTHS
      KXZ=NMODE*NMODE*3
      KR0=NATOM*NATOM
      CALL MEMO(2,LXZ,KXZ,LR0,KR0,0,0,0,0,0,0)
      CALL MEMO(1,LYZ,KXZ,0,0,0,0,0,0,0,0)
      CALL RZETA(NATOM,NMODE,W(LXL),W(LXZ),W(LX0),W(LR0),W(LXM))
      CALL FLUSH(IOUT)
C**LOOP OVER VSCF STATES TO BE CALCULATED
      KESCF=NSTAT
      CALL MEMO(1,LESCF,KESCF,0,0,0,0,0,0,0,0)
      WRITE(IOUT,325)
      WRITE(IOUT,325)
      WRITE(IOUT,220)
      CALL FLUSH(IOUT)
      ITIM=-1
      J21=2*JMAX+1
      CALL MEMO(2,LSSX,J21*9*J21,LSS,J21*9*J21,0,0,0,0,0,0)
      CALL MEMO(3,LW21,J21,LX21,J21*J21,LE21,J21,0,0,0,0)
C**GET ROTATIONAL MATRIX ELEMENTS FOR K-DIAGONAL
      CALL ROTEL(J21-1,W(LSS),W(LSSX),J21,1,1,1)
C**GET EQUILIBRIUM ROTATIONAL CONSTANTS
C**TEMPORARY...LINEAR
      IF(.NOT.LINEAR)CALL GETMI0(NATOM,W(LX0),W(LXM))
C**TEMPORARY...LINEAR
      IF(MATSIZ.NE.0.AND.ICI.LT.0)GO TO 4000
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                            DUMP POTENTIAL AND CORIOLIS DATA TO DISC
C**********************************************************************
C**********************************************************************
C**********************************************************************
C******************************
      IF(ICOUPL.EQ.0)GO TO 3000
      MBFMX1=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(IK2.GT.MBFMX1)MBFMX1=IK2
      END DO
      MBFMAX=MBFMX1
      IF(JCOUPL.GT.0)THEN
        KLV1=MBFMX1
        KLC1=MBFMX1
      ELSE
        KLV1=(1+MBFMX1)/2
        KLC1=(1+MBFMX1)/2
      END IF
      CALL MEMO(2,LV1,KLV1,LC1,KLC1,0,0,0,0,0,0)
      IF(JMAX.GT.0)THEN
        IF(JCOUPL.GT.0)THEN
          KEJK1=MBFMX1*J21
        ELSE
          KEJK1=(1+MBFMX1*J21)/2
        END IF
      ELSE
        KEJK1=1
      END IF
      CALL MEMO(1,LEJK1,KEJK1,0,0,0,0,0,0,0,0)
CCCC  IF(IDISC.EQ.0.OR.JMAX.GT.0)THEN
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
CCCC      IF(IDISC.EQ.0)THEN
C**WRITE POTENTIAL GRIDS TO DISC
            CALL DUMPT1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     1      W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,
     2      W(LV1),W(LV1),1)
C**WRITE CORIOLIS GRIDS TO DISC
            CALL DUMCR1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     1      W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,W(LC1),
     2      W(LC1),W(LVM1),W(LVM1),0)
CCCC      END IF
          IF(JMAX.GT.0.AND.ICOUPL.EQ.1)
     1    CALL GETMI1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),
     1    W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     2    W(LXX),W(LX0),
     2    W(LXL),W(LXM),K,W(LEJK1),W(LEJK1),W(LSS),J21,
     3    W(LX21),W(LW21),W(LE21))
3001      CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.1)GO TO 3000
      MBFMX2=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO L=1,K-1
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          IF(IK2*IL2.GT.MBFMX2)MBFMX2=IK2*IL2
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV2=MBFMX2
        KLC2=MBFMX2*6
      ELSE
        KLV2=(1+MBFMX2)/2
        KLC2=(1+MBFMX2*6)/2
      END IF
      CALL MEMO(2,LV2,KLV2,LC2,KLC2,0,0,0,0,0,0)
      IF(JMAX.GT.0)THEN
        IF(JCOUPL.GT.0)THEN
          KEJK2=MBFMX2*J21
        ELSE
          KEJK2=(1+MBFMX2*J21)/2
        END IF
      ELSE
        KEJK2=1
      END IF
      CALL MEMO(1,LEJK2,KEJK2,0,0,0,0,0,0,0,0)
CCCC  IF(IDISC.EQ.0.OR.JMAX.GT.0)THEN
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
CCCC        IF(IDISC.EQ.0)THEN
C**WRITE POTENTIAL GRIDS TO DISC
              CALL DUMPT2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1        W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),
     2        W(LJPOT),W(LCPOT),K,L,W(LV2),W(LV2))
C**WRITE CORIOLIS GRIDS TO DISC
              CALL DUMCR2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1        W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),
     2        W(LXL),W(LXM),K,L,W(LC2),W(LC2),W(LVM2),W(LVM2),0)
CCCC        END IF
            IF(JMAX.GT.0.AND.ICOUPL.EQ.2)
     1      CALL GETMI2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,W(LQQ),
     1      W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     2      W(LXX),W(LX0),W(LXL),W(LXM),K,L,W(LEJK2),W(LEJK2),W(LSS),
     3      J21,W(LX21),W(LW21),W(LE21))
3022        CONTINUE
            L3=L3+IL2
          END DO
3002      CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.2)GO TO 3000
      MBFMX3=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO L=1,K-1
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          DO N=1,L-1
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
            IF(IK2*IL2*IN2.GT.MBFMX3)MBFMX3=IK2*IL2*IN2
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV3=MBFMX3
        KLC3=MBFMX3*10
      ELSE
        KLV3=(1+MBFMX3)/2
        KLC3=(1+MBFMX3*10)/2
      END IF
      CALL MEMO(2,LV3,KLV3,LC3,KLC3,0,0,0,0,0,0)
      IF(JMAX.GT.0)THEN
        IF(JCOUPL.GT.0)THEN
          KEJK3=MBFMX3*J21
        ELSE
          KEJK3=(1+MBFMX3*J21)/2
        END IF
      ELSE
        KEJK3=1
      END IF
      CALL MEMO(1,LEJK3,KEJK3,0,0,0,0,0,0,0,0)
CCCC  IF(IDISC.EQ.0.OR.JMAX.GT.0)THEN
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
CCCC          IF(IDISC.EQ.0)THEN
C**WRITE POTENTIAL GRIDS TO DISC
                CALL DUMPT3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,IN2,
     1          NMODE,NATOM,W(LQQ),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     2          NPOT,W(LIPOT),W(LJPOT),W(LCPOT),K,L,N,W(LV3),W(LV3))
C**WRITE CORIOLIS GRIDS TO DISC
                CALL DUMCR3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,IN2,
     1          NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     2          W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,W(LC3),W(LC3),
     3          W(LVM3),W(LVM3),0)
CCCC          END IF
              IF(JMAX.GT.0.AND.ICOUPL.EQ.3)
     1        CALL GETMI3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,IN2,
     2        NMODE,NATOM,W(LQQ),
     1        W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     2        W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,
     3        W(LEJK3),W(LEJK3),W(LSS),J21,
     4        W(LX21),W(LW21),W(LE21))
3333          CONTINUE
              N3=N3+IN2
            END DO
3033        CONTINUE
            L3=L3+IL2
          END DO
3003      CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.3)GO TO 3000
      MBFMX4=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        DO L=1,K-1
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
          DO N=1,L-1
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
            DO M=1,N-1
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
              IF(IK2*IL2*IN2*IM2.GT.MBFMX4)MBFMX4=IK2*IL2*IN2*IM2
            END DO
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        KLV4=MBFMX4
        KLC4=MBFMX4*15
      ELSE
        KLV4=(1+MBFMX4)/2
        KLC4=(1+MBFMX4*15)/2
      END IF
      CALL MEMO(2,LV4,KLV4,LC4,KLC4,0,0,0,0,0,0)
      IF(JMAX.GT.0)THEN
        IF(JCOUPL.GT.0)THEN
          KEJK4=MBFMX4*J21
        ELSE
          KEJK4=(1+MBFMX4*J21)/2
        END IF
      ELSE
        KEJK4=1
      END IF
      CALL MEMO(1,LEJK4,KEJK4,0,0,0,0,0,0,0,0)
CCCC  IF(IDISC.EQ.0.OR.JMAX.GT.0)THEN
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M3=0
              DO M=1,N-1
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
CCCC            IF(IDISC.EQ.0)THEN
C**WRITE POTENTIAL GRIDS TO DISC
                  CALL DUMPT4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1            IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LRR),W(LXX),
     2            W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     3            K,L,N,M,W(LV4),W(LV4))
C**WRITE CORIOLIS GRIDS TO DISC
                  CALL DUMCR4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1            IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     2            W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,L,
     3            N,M,W(LC4),W(LC4),W(LVM4),W(LVM4),0)
CCCC            END IF
                IF(JMAX.GT.0.AND.ICOUPL.EQ.4)
     1          CALL GETMI4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     2          IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),
     1          W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     2          W(LXX),W(LX0),
     3          W(LXL),W(LXM),K,L,N,M,W(LEJK4),W(LEJK4),W(LSS),J21,
     4          W(LX21),W(LW21),W(LE21))
34444           CONTINUE
                M3=M3+IM2
              END DO
30444         CONTINUE
              N3=N3+IN2
            END DO
30044       CONTINUE
            L3=L3+IL2
          END DO
30004     CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.4)GO TO 3000
C**5-MODE AND HIGHER
3000  CONTINUE
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                                      VSCF PROCEDURE
C**********************************************************************
C**********************************************************************
C**********************************************************************
      DO 1500 ISTATE=1,NSTAT
C**********************************************************************
C**                                        LOOP OVER SCF STATES IN TURN
C**********************************************************************
C     IF(ITIM.LE.0)THEN
C       WRITE(IOUT,*)'Calculating State ',ISTATE
C       CALL TIMIT(1)
C     END IF
      ITIM=ITIM+1
C**DO SCF CYCLE FOR J=0 FIRST TIME
      J21=1
      JTHIS=0
      WRITE(IOUT,320)
      WRITE(IOUT,310)
      WRITE(IOUT,320)
8000  CONTINUE
      count=0.0
      IF(J21.EQ.1)E0=1.D+20
      KCORIG=JTHIS
      DO 7000 KROT=1,J21
C**********************************************************************
C**                                      LOOP OVER ROTATIONAL FUNCTIONS
C**********************************************************************
C**KA,KC FOR VSCF
      KA=KROT/2
      KC=KCORIG-KA
C     IF(MOD(KA,2).NE.0)THEN
C       KC=KC+MOD(KROT,2)
C     ELSE
C       KC=KC+MOD(KROT+1,2)
C     END IF
      IF(MOD(KROT,2).EQ.0)KC=KC+1
      ICYCLE=1
500   CONTINUE
C**********************************************************************
C**                               RETURN TO LABEL 500 UNTIL CONVERGENCE
C**********************************************************************
      ITIM=ITIM+1
C**FIRST STATE IS ZERO POINT LEVEL...DISPLAY ZERO POINT ENERGY
C**ENERGIES OF REMAINING STATES ARE RELATIVE TO ZERO POINT
C**OF PURE VIBRATIONAL LEVEL
      IF(ISTATE.EQ.1.AND.J21.EQ.1)EVL=0
      MNLK2=0
C     MNLK3=0
      REWIND 58
      REWIND 59
      DO 1000 MODE=1,NMODE
C**********************************************************************
C**                 LOOP OVER ALL MODES....FORM NEW FUNCTIONS THIS MODE
C**********************************************************************
      CALL INTARR(W(LNBF),W(LMBF),MODE,IMODE1,IMODE2,IMODE3)
      KXK=IMODE3
      CALL MEMO(3,LXK,KXK*KXK,LYK,KXK*KXK,LOV,KXK*KXK,0,0,0,0)
      CALL MTZERO(W(LXK),KXK)
      IF(ICOUPL.GT.0)THEN
        IF(J21.GT.1)REWIND 61
        REWIND 71
        REWIND 81
      END IF
      IF(ICOUPL.GT.1)THEN
        IF(J21.GT.1)REWIND 62
        REWIND 72
        REWIND 82
      END IF
      IF(ICOUPL.GT.2)THEN
        IF(J21.GT.1)REWIND 63
        REWIND 73
        REWIND 83
      END IF
      IF(ICOUPL.GT.3)THEN
        IF(J21.GT.1)REWIND 64
        REWIND 74
        REWIND 84
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(MODE.EQ.1.AND.K.EQ.1.AND.ITIM.EQ.0)THEN
          ITIM1A=0
          ITIM1B=0
        END IF
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        CALL MEMO(1,LXK0,IK2*4,0,0,0,0,0,0,0,0)
C**'KINETIC ENERGY'
        IF(K.EQ.MODE)THEN
          CALL THISKE(W(LH+K2),W(LXQ+K3),IK3,IK2,W(LXK),
     1    W(LOMEGA+K1))
        ELSE
          CALL THATKE(W(LISTAT),NSTAT,NMODE,ISTATE,K,W(LH+K2),
     1    W(LXQ+K3),IK3,IK2,W(LXK),KXK,W(LOMEGA+K1))
        END IF
        IF(ICOUPL.EQ.0)GO TO 501
C**CORIOLIS AND POTENTIAL
        IF(K.EQ.MODE)THEN
          CALL THIS1(W(LH+K2),W(LXQ+K3),IK3,IK2,W(LXK),
     1    NMODE,K,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR)
     2    ,W(LXX),W(LX0),W(LXL),W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     3    W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),W(LEJK1),J21,KROT)
        ELSE
          CALL THAT1(W(LISTAT),NSTAT,NMODE,ISTATE,K,W(LH+K2),
     1    W(LXQ+K3),IK3,IK2,W(LXK),KXK,NATOM,W(LQQ),W(LXZ),W(LAB),
     2    W(LB),W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     3    NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV1),W(LV1),W(LC1),W(LC1),
     4    W(LEJK1),W(LEJK1),J21,KROT)
        END IF
        IF(ICOUPL.EQ.1)GO TO 501
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(MODE.EQ.1.AND.L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)THEN
            ITIM2A=0
            ITIM2B=0
          END IF
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**CORIOLIS AND POTENTIAL
          CALL MEMO(1,LXL0,IL2*4,0,0,0,0,0,0,0,0)
          IF(L.EQ.MODE)THEN
            CALL THIS2B(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),W(LH+L2),
     1      W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXK0),W(LXL0),
     2      W(LXK),NMODE,K,L,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     3      W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),ISTATE,
     4      NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV2),W(LV2),W(LC2),
     5      W(LC2),W(LEJK2),W(LEJK2),J21,KROT)
          ELSE
            IF(K.EQ.MODE)
     1      CALL THIS2A(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),W(LH+L2),
     2      W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXK0),W(LXL0),
     3      W(LXK),NMODE,K,L,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     4      W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),ISTATE,
     5      NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV2),W(LV2),W(LC2),
     6      W(LC2),W(LEJK2),W(LEJK2),J21,KROT)
            IF(K.NE.MODE)
     1      CALL THAT2(W(LISTAT),NSTAT,NMODE,ISTATE,K,L,W(LH+K2),
     2      W(LXQ+K3),W(LH+L2),W(LXQ+L3),IK3,IK2,IL3,IL2,
     3      W(LXK0),W(LXL0),W(LXK),KXK,NATOM,W(LQQ),W(LXZ),W(LAB),
     4      W(LB),W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     5      NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV2),W(LV2),W(LC2),
     6      W(LC2),W(LEJK2),W(LEJK2),J21,KROT)
          END IF
          IF(ICOUPL.EQ.2)GO TO 502
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(MODE.EQ.1.AND.N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.
     1      ITIM.EQ.0)THEN
              ITIM3A=0
              ITIM3B=0
            END IF
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**CORIOLIS AND POTENTIAL
            CALL MEMO(1,LXN0,IN2*4,0,0,0,0,0,0,0,0)
            IF(N.EQ.MODE)THEN
              CALL THIS3C(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),W(LH+L2),
     1        W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,
     2        IL3,IL2,IN3,IN2,W(LXK0),W(LXL0),W(LXN0),
     3        W(LXK),NMODE,K,L,N,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     4        W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),ISTATE,
     5        NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV3),W(LV3),W(LC3),
     6        W(LC3),W(LEJK3),W(LEJK3),J21,KROT)
            ELSE
              IF(L.EQ.MODE)THEN
                CALL THIS3B(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),W(LH+L2)
     1          ,W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,
     2          IL3,IL2,IN3,IN2,W(LXK0),W(LXL0),W(LXN0),
     3          W(LXK),NMODE,K,L,N,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     4          W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),ISTATE
     5          ,NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV3),W(LV3),W(LC3),
     6          W(LC3),W(LEJK3),W(LEJK3),J21,KROT)
              ELSE
                IF(K.EQ.MODE)
     1          CALL THIS3A(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),W(LH+L2)
     2          ,W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,
     3          IL3,IL2,IN3,IN2,W(LXK0),W(LXL0),W(LXN0),
     4          W(LXK),NMODE,K,L,N,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     5          W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),ISTATE
     6          ,NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV3),W(LV3),W(LC3),
     7          W(LC3),W(LEJK3),W(LEJK3),J21,KROT)
                IF(K.NE.MODE)
     1          CALL THAT3(W(LISTAT),NSTAT,NMODE,ISTATE,K,L,N,W(LH+K2),
     2          W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),
     3          IK3,IK2,IL3,IL2,IN3,IN2,W(LXK0),W(LXL0),W(LXN0),
     4          W(LXK),KXK,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),
     5          W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     6          NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV3),W(LV3),W(LC3),
     7          W(LC3),W(LEJK3),W(LEJK3),J21,KROT)
              END IF
            END IF
            IF(ICOUPL.EQ.3)GO TO 503
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(MODE.EQ.1.AND.M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.
     1        K.EQ.4.AND.ITIM.EQ.0)THEN
                ITIM4A=0
                ITIM4B=0
              END IF
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**CORIOLIS AND POTENTIAL
              CALL MEMO(1,LXM0,IM2*4,0,0,0,0,0,0,0,0)
              IF(M.EQ.MODE)THEN
                CALL THIS4D(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),W(LH+L2)
     1          ,W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),IK3,
     2          IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXK0),W(LXL0),W(LXN0),
     3          W(LXM0),
     4          W(LXK),NMODE,K,L,N,M,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),
     5          W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),ISTATE
     6          ,NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV4),W(LV4),W(LC4),
     7          W(LC4),W(LEJK4),W(LEJK4),J21,KROT)
              ELSE
                IF(N.EQ.MODE)THEN
                  CALL THIS4C(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1            W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2            W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXK0),
     3            W(LXL0),W(LXN0),W(LXM0),
     4            W(LXK),NMODE,K,L,N,M,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB)
     5            ,W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),W(LXM),
     6            ISTATE,NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV4),W(LV4),
     7            W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,KROT)
                ELSE
                  IF(L.EQ.MODE)THEN
                    CALL THIS4B(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     1              W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2              W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXK0),
     3              W(LXL0),W(LXN0),W(LXM0),
     4              W(LXK),NMODE,K,L,N,M,NATOM,W(LQQ),W(LXZ),W(LAB),
     5              W(LB),W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),
     6              W(LXM),ISTATE,NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     7              W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
     8              KROT)
                  ELSE
                    IF(K.EQ.MODE)
     1              CALL THIS4A(W(LISTAT),NSTAT,W(LH+K2),W(LXQ+K3),
     2              W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     3              W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXK0),
     4              W(LXL0),W(LXN0),W(LXM0),
     5              W(LXK),NMODE,K,L,N,M,NATOM,W(LQQ),W(LXZ),W(LAB),
     6              W(LB),W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),
     7              W(LXM),ISTATE,NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     8              W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
     9              KROT)
                    IF(K.NE.MODE)
     1              CALL THAT4(W(LISTAT),NSTAT,NMODE,ISTATE,K,L,N,M,
     2              W(LH+K2),W(LXQ+K3),W(LH+L2),W(LXQ+L3),W(LH+N2),
     3              W(LXQ+N3),W(LH+M2),W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,
     4              IN2,IM3,IM2,W(LXK0),W(LXL0),W(LXN0),W(LXM0),
     5              W(LXK),KXK,NATOM,W(LQQ),W(LXZ),W(LAB),
     6              W(LB),W(LAA),W(LBB),W(LRR),W(LXX),W(LX0),W(LXL),
     7              W(LXM),NPOT,W(LIPOT),W(LJPOT),W(LCPOT),
     8              W(LV4),W(LV4),W(LC4),W(LC4),W(LEJK4),W(LEJK4),J21,
     9              KROT)
                  END IF
                END IF
              END IF
              IF(ICOUPL.EQ.4)GO TO 504
C**5-MODE COUPLING HERE IF NEEDED
504   CONTINUE
              CALL MEMO(-1,LXM0,IM2*4,0,0,0,0,0,0,0,0)
5004          CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
503   CONTINUE
            CALL MEMO(-1,LXN0,IN2*4,0,0,0,0,0,0,0,0)
5003        CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
502   CONTINUE
          CALL MEMO(-1,LXL0,IL2*4,0,0,0,0,0,0,0,0)
5002      CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
501   CONTINUE
        CALL MEMO(-1,LXK0,IK2*4,0,0,0,0,0,0,0,0)
5001    CONTINUE
        K2=K2+IK1
        K3=K3+IK2
      END DO
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
C**DIAGONALISE FOR THIS MODE..NO NEED FOR ENERGIES UNTIL FINAL MODE
      KWK=KXK
      KEVAL=KXK
      CALL MEMO(2,LWK,KWK,LEVAL,KEVAL,0,0,0,0,0,0)
      IF(ITIM.EQ.1)THEN
        WRITE(IOUT,*)'Calculating DIAG'
        CALL TIMIT(1)
      END IF
      IF(MODE.NE.NMODE)THEN
        CALL DIAG(W(LXK),W(LXK),KXK,KXK,-1,W(LWK),W(LEVAL),W(LWK),KXK,
     1  KXK,W(LISTAT),NSTAT,NMODE,W(LXK),W(LXK),W(LXK),IDUM,IDUM,IDUM)
      ELSE
        CALL DIAG(W(LXK),W(LXK),KXK,KXK,0,W(LWK),W(LEVAL),W(LWK),
     1  KXK,KXK,W(LISTAT),NSTAT,NMODE,W(LXK),W(LXK),W(LXK),IDUM,IDUM,
     2  IDUM)
C**FIND ENERGY FOR THIS STATE....CORRESPONDS TO QUANTUM OF FINAL MODE
        CALL ENERGY(W(LISTAT),NSTAT,NMODE,ISTATE,MODE,W(LEVAL),E,
     1  W(LESCF+ISTATE-1),W(LWK))
      END IF
      IF(ITIM.EQ.1)CALL TIMIT(3)
      CALL FLUSH(IOUT)
      CALL MEMO(-1,LEVAL,KEVAL,0,0,0,0,0,0,0,0)
C**FORM NEW FUNCTIONS THIS MODE
      CALL REFORM(W(LH+MNLK2),W(LXK),IMODE3,IMODE3,IMODE2,W(LWK),
     1W(LYK),W(LOV))
      MNLK2=MNLK2+IMODE1
C     MNLK3=MNLK3+IMODE2
      CALL MEMO(-4,LXK,KXK*KXK,LWK,KWK,LYK,KXK*KXK,LOV,KXK*KXK,0,0)
1000  CONTINUE
C**********************************************************************
C**                                           TEST FOR SELF-CONSISTENCY
C**********************************************************************
      IF(ISTATE.EQ.1.AND.J21.EQ.1)THEN
        CALL PRSCF(W(LISTAT),NSTAT,NMODE,ISTATE,E0,EVL,JTHIS,KA,KC,0)
        E0E=E0-EVL
      ELSE
        CALL PRSCF(W(LISTAT),NSTAT,NMODE,ISTATE,E0,E,JTHIS,KA,KC,0)
        E0E=E0-E
      END IF
c  change convergence criterion to converge abs diff in cm-1
c  to conv, eg. conv = 5e-2 cm-2.  5/12/97
      if(count.gt.15.0)then
        write(IOUT,*)ISTATE,'   did not converge'
        go to 999
      end if
      count=count + 1.0
      IF(dabs(e0e).gt.conv)then
        IF(ISTATE.EQ.1.AND.J21.EQ.1)E0=EVL
        IF(ISTATE.NE.1.OR.J21.NE.1)E0=E
        CALL FLUSH(IOUT)
C**MAKE COPY (59) -> (58)
        IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)THEN
          REWIND 58
          REWIND 59
          DO MODE=1,NMODE
            CALL INTARR(W(LNBF),W(LMBF),MODE,IMODE1,IMODE2,IMODE3)
            CALL MEMO(1,LWK,IMODE3,0,0,0,0,0,0,0,0)
            CALL RECYCL(W(LWK),IMODE3,59,58)
            CALL MEMO(-1,LWK,IMODE3,0,0,0,0,0,0,0,0)
          END DO
        END IF
        ICYCLE=ICYCLE+1
        GO TO 500
      END IF
C**********************************************************************
C**                                                 SCF STATE CONVERGED
C**********************************************************************
999   continue
      IF(ISTATE.EQ.1.AND.J21.EQ.1)THEN
        CALL PRSCF(W(LISTAT),NSTAT,NMODE,ISTATE,E0,EVL,JTHIS,KA,KC,40)
      ELSE
        CALL PRSCF(W(LISTAT),NSTAT,NMODE,ISTATE,E0,E,JTHIS,KA,KC,40)
      END IF
C**VIBRATIONAL GROUND STATE K=0 SCF FUNCTION ONLY IF VIRTUALS
      IF(ICI.LT.0.AND.J21.NE.1)GO TO 1600
      IF(J21.NE.1.AND.KROT.NE.J21)WRITE(IOUT,320)
7000  CONTINUE
C**********************************************************************
C**                                          TIDY UP FOR THIS SCF STATE
C**********************************************************************
      IF(J21.EQ.1)THEN
        IF(JMAX.NE.0)THEN
          J21=2*JMAX+1
          JTHIS=JMAX
          WRITE(IOUT,320)
          WRITE(IOUT,315)JMAX
          WRITE(IOUT,320)
          GO TO 8000
        END IF
      END IF
      IF(ISCFCI.NE.0)THEN
        K2=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
C**USE SCF FUNCTIONS
          IF(ICI.GT.0)
     1    CALL OUT50(W(LH+K2),I3,I2,W(LISTAT),NSTAT,NMODE,
     2    ISTATE,K)
C**USE VIRTUAL FUNCTIONS FOR VIBRATIONAL GROUND STATE IN CORE
          K2=K2+I1
        END DO
      END IF
C**FINISHED THIS STATE
C     CALL TIMIT(3)
      ITIM=0
C**VIBRATIONAL GROUND STATE SCF FUNCTION ONLY IF VIRTUALS
      IF(ICI.LT.0)GO TO 1600
C**ONLY NEED FIRST ICI SCF STATES IF CI
      IF(ISTATE.EQ.ICI)GO TO 1600
1500  CONTINUE
      WRITE(IOUT,215)
      RETURN
C**********************************************************************
C**                                                  RETURN IF SCF ONLY
C**********************************************************************
1600  CONTINUE
C**MAKE COPY (59) -> (60) AND WRITE VCI FUNCTIONS + 1ST, 2ND DERIVS.
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)THEN
        REWIND 59
        K2=0
        DO MODE=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),MODE,IMODE1,IMODE2,IMODE3)
          CALL MEMO(1,LWK,IMODE3,0,0,0,0,0,0,0,0)
          CALL RECYCL(W(LWK),IMODE3,59,60)
          CALL OUTVCI(IMODE3,IMODE2,W(LH+K2),W(LWK))
          CALL MEMO(-1,LWK,IMODE3,0,0,0,0,0,0,0,0)
          K2=K2+IMODE1
        END DO
      END IF
      CALL FLUSH(IOUT)
      IF(ISCFCI.EQ.0)THEN
        WRITE(IOUT,215)
        RETURN
      ELSE
        WRITE(IOUT,325)
        WRITE(IOUT,325)
        WRITE(IOUT,235)
      END IF
C**********************************************************************
C**********************************************************************
C**********************************************************************
C**                                          VSCF-CI AND VCI PROCEDURES
C**********************************************************************
C**********************************************************************
C**********************************************************************
      IF(ICI.LT.0)GO TO 4000
C**********************************************************************
C**********************************************************************
C**                                                             VSCF-CI
C**********************************************************************
C**********************************************************************
C**GET FUNCTIONS SPACE FOR TWO SCF STATES
      KHLR=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        KHLR=KHLR+I2*3
      END DO
      CALL MEMO(2,LHL,KHLR,LHR,KHLR,0,0,0,0,0,0)
4000  CONTINUE
      IF(ICI.GT.0)THEN
C**********************************************************************
C**********************************************************************
C**                                                             VSCF-CI
C**********************************************************************
C**********************************************************************
        ISIZE=ICI
        JSIZE=ISIZE
        KIP=ISIZE*NMODE
        KJP=KIP
        CALL MEMO(2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
C**GET CI BASIS HERE IN IP
        CALL MOVEIP(W(LISTAT),NSTAT,W(LIP),ISIZE,NMODE)
        CALL GETIP(W(LIP),W(LJP),ISIZE,JSIZE,NMODE,ICI,ICI)
        WRITE(IOUT,395)(NTOT(K),K=1,NVSYM)
        CALL FLUSH(IOUT)
      END IF
C*******************************************************************BUG
      NBFMIN=100000
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(IK3.LT.NBFMIN)NBFMIN=IK3
      END DO
      IF(ICI.LT.0.AND.NMAX.GE.0)THEN
        JCI=-ICI
        IF(NMAX.GT.0)THEN
C**ONLY GET NMAX +1 (MAXIMUM) INTEGRALS
          IF(NMAX+1.LT.JCI)JCI=NMAX+1
C**OR NBFMIN (MAXIMUM) INTEGRALS (WHICHEVER IS SMALLER)
          IF(JCI.GT.NBFMIN)JCI=NBFMIN
C**CONVERT NMAX FROM QUANTA TO NUMBER OF FUNCTIONS
          NMAX=NMAX+NMODE
        END IF
      END IF
C*******************************************************************BUG
      IF(ICI.LT.0)THEN
C**********************************************************************
C**********************************************************************
C**                                                                 VCI
C**********************************************************************
C**********************************************************************
        NVSYMX=NVSYM
        NVSYM=1
C**NNMAX USED AS INDICATOR AS TO WHICH ALGORITHM
        NNMAX=NMAX
        IF(NNMAX.GE.0)THEN
C**NMAX SAVED IN NMAXMX (RE-LOADED LATER)
          NMAXMX=NMAX
          DO I=1,4
            NBAS(I)=JCI-1
            MAXSUM(I)=-I
          END DO
        ELSE
          JCI1=NBAS(1)+1
          IF(JCI1.GT.NBFMIN)JCI1=NBFMIN
          NBAS(1)=JCI1-1
          JCI2=NBAS(2)+1
          IF(JCI2.GT.NBFMIN)JCI2=NBFMIN
          NBAS(2)=JCI2-1
          JCI3=NBAS(3)+1
          IF(JCI3.GT.NBFMIN)JCI3=NBFMIN
          NBAS(3)=JCI3-1
          JCI4=NBAS(4)+1
          IF(JCI4.GT.NBFMIN)JCI4=NBFMIN
          NBAS(4)=JCI4-1
          JCI=MAX0(JCI1,JCI2,JCI3,JCI4)
          JCIMAX=JCI
C**TEMPORARY USE OF NMAX HOLDS MAX. POSSIBLE MAXSUM
          NMAXMX=MAX0(MAXSUM(1),MAXSUM(2),MAXSUM(3),MAXSUM(4))
        END IF
        IF(ICOUPL.GT.0.OR.(NNMAX.LT.0.AND.ICI.LT.0))THEN
          NMAX=NMAXMX+1
          ISIZE1=JCI
          JSIZE1=ISIZE1
C**FIRST GET SIZE (ISIZE1 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE1,1,JCI)
          KIP1=ISIZE1
          KJP1=KIP1
          CALL MEMO(2,LIP1,KIP1,LJP1,KJP1,0,0,0,0,0,0)
C**GET BASIC BASIS HERE IN IP1
          CALL GETIP(W(LIP1),W(LJP1),ISIZE1,JSIZE1,1,JCI,ICI)
          CALL PUTJP(W(LJP1),ISIZE1,W(LIP1),ISIZE1,1,1)
          CALL MEMO(-1,LJP1,KJP1,0,0,0,0,0,0,0,0)
        END IF
        IF(ICOUPL.GT.1.OR.(NNMAX.LT.0.AND.ICI.LT.-1))THEN
          NMAX=NMAXMX+2
          ISIZE2=JCI**2
          JSIZE2=ISIZE2
C**FIRST GET SIZE (ISIZE2 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE2,2,JCI)
          KIP2=ISIZE2*2
          KJP2=KIP2
          CALL MEMO(2,LIP2,KIP2,LJP2,KJP2,0,0,0,0,0,0)
C**GET BASIC BASIS HERE IN IP2
          CALL GETIP(W(LIP2),W(LJP2),ISIZE2,JSIZE2,2,JCI,ICI)
          CALL PUTJP(W(LJP2),ISIZE2,W(LIP2),ISIZE2,2,1)
          CALL MEMO(-1,LJP2,KJP2,0,0,0,0,0,0,0,0)
        END IF
        IF(ICOUPL.GT.2.OR.(NNMAX.LT.0.AND.ICI.LT.-2))THEN
          NMAX=NMAXMX+3
          ISIZE3=JCI**3
          JSIZE3=ISIZE3
C**FIRST GET SIZE (ISIZE3 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE3,3,JCI)
          KIP3=ISIZE3*3
          KJP3=KIP3
          CALL MEMO(2,LIP3,KIP3,LJP3,KJP3,0,0,0,0,0,0)
C**GET BASIC BASIS HERE IN IP3
          CALL GETIP(W(LIP3),W(LJP3),ISIZE3,JSIZE3,3,JCI,ICI)
          CALL PUTJP(W(LJP3),ISIZE3,W(LIP3),ISIZE3,3,1)
          CALL MEMO(-1,LJP3,KJP3,0,0,0,0,0,0,0,0)
        END IF
        IF(ICOUPL.GT.3.OR.(NNMAX.LT.0.AND.ICI.LT.-3))THEN
          NMAX=NMAXMX+4
          ISIZE4=JCI**4
          JSIZE4=ISIZE4
C**FIRST GET SIZE (ISIZE4 MODIFIED IN GETSZ FOR NMAX)
          CALL GETSZ(ISIZE4,4,JCI)
          KIP4=ISIZE4*4
          KJP4=KIP4
          CALL MEMO(2,LIP4,KIP4,LJP4,KJP4,0,0,0,0,0,0)
C**GET BASIC BASIS HERE IN IP4
          CALL GETIP(W(LIP4),W(LJP4),ISIZE4,JSIZE4,4,JCI,ICI)
          CALL PUTJP(W(LJP4),ISIZE4,W(LIP4),ISIZE4,4,1)
          CALL MEMO(-1,LJP4,KJP4,0,0,0,0,0,0,0,0)
        END IF
        NVSYM=NVSYMX
        IF(NNMAX.GE.0)THEN
C**GET FULL CI BASIS HERE IN IP
          NMAX=NMAXMX
          ISIZE=(JCI)**NMODE
          JSIZE=ISIZE
C**FIRST GET SIZE (ISIZE MODIFIED IN GETSZ FOR NMAX)
          IF(NMAX.GT.0)CALL GETSZ(ISIZE,NMODE,JCI)
          KIP=ISIZE*NMODE
          KJP=KIP
          CALL MEMO(2,LIP,KIP,LJP,KJP,0,0,0,0,0,0)
          CALL GETIP(W(LIP),W(LJP),ISIZE,JSIZE,NMODE,JCI,ICI)
C**RESET JSIZE
          JSIZE=ISIZE
          WRITE(IOUT,395)(NTOT(K),K=1,NVSYM)
          CALL FLUSH(IOUT)
        ELSE
C**GET SELECTIVE CI BASIS HERE IN IP
          ITOT=0
          ITOT1=1
          ITOT2=0
          ITOT3=0
          ITOT4=0
          IF(ICI.LT.0)THEN
            JCI1=NBAS(1)+1
            NMAX=MAXSUM(1)+1
            KSIZE1=JCI1
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
            CALL GETSZ(KSIZE1,1,JCI1)
C**SUBTRACT ZERO POINT BASIS
            KSIZE1=KSIZE1-1
C**ADD 1 TO ONE-MODE MATRIX FOR ZERO POINT BASIS
            JSIZE1=KSIZE1*NMODE+1
            KJP1=JSIZE1*NMODE
            IF(KJP1.GT.0)
     1      CALL MEMO(1,LJP1,KJP1,0,0,0,0,0,0,0,0)
            IF(ICI.LT.-1)THEN
              JCI2=NBAS(2)+1
              NMAX=MAXSUM(2)+1
              KSIZE1=JCI2
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
              CALL GETSZ(KSIZE1,1,JCI2)
              KSIZE1=KSIZE1-1
              NMAX=MAXSUM(2)+2
              KSIZE2=JCI2*JCI2
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
              CALL GETSZ(KSIZE2,2,JCI2)
              KSIZE2=KSIZE2-2*KSIZE1-1
              JSIZE2=KSIZE2*NMODE*(NMODE-1)/2
              KJP2=JSIZE2*NMODE
              IF(KJP2.GT.0)
     1        CALL MEMO(1,LJP2,KJP2,0,0,0,0,0,0,0,0)
              IF(ICI.LT.-2)THEN
                JCI3=NBAS(3)+1
                NMAX=MAXSUM(3)+1
                KSIZE1=JCI3
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(KSIZE1,1,JCI3)
                KSIZE1=KSIZE1-1
                NMAX=MAXSUM(3)+2
                KSIZE2=JCI3*JCI3
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(KSIZE2,2,JCI3)
                KSIZE2=KSIZE2-2*KSIZE1-1
                NMAX=MAXSUM(3)+3
                KSIZE3=JCI3*JCI3*JCI3
C**FIRST GET SIZE (KSIZE3 MODIFIED IN GETSZ FOR NMAX)
                CALL GETSZ(KSIZE3,3,JCI3)
                KSIZE3=KSIZE3-3*KSIZE2-3*KSIZE1-1
                JSIZE3=KSIZE3*NMODE*(NMODE-1)*(NMODE-2)/(3*2)
                KJP3=JSIZE3*NMODE
                IF(KJP3.GT.0)
     1          CALL MEMO(1,LJP3,KJP3,0,0,0,0,0,0,0,0)
                IF(ICI.LT.-3)THEN
                  JCI4=NBAS(4)+1
                  NMAX=MAXSUM(4)+1
                  KSIZE1=JCI4
C**FIRST GET SIZE (KSIZE1 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE1,1,JCI4)
                  KSIZE1=KSIZE1-1
                  NMAX=MAXSUM(4)+2
                  KSIZE2=JCI4*JCI4
C**FIRST GET SIZE (KSIZE2 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE2,2,JCI4)
                  KSIZE2=KSIZE2-2*KSIZE1-1
                  NMAX=MAXSUM(4)+3
                  KSIZE3=JCI4*JCI4*JCI4
C**FIRST GET SIZE (KSIZE3 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE3,3,JCI4)
                  KSIZE3=KSIZE3-3*KSIZE2-3*KSIZE1-1
                  NMAX=MAXSUM(4)+4
                  KSIZE4=JCI4*JCI4*JCI4*JCI4
C**FIRST GET SIZE (KSIZE4 MODIFIED IN GETSZ FOR NMAX)
                  CALL GETSZ(KSIZE4,4,JCI4)
                  KSIZE4=KSIZE4-4*KSIZE3-6*KSIZE2-4*KSIZE1-1
                  JSIZE4=KSIZE4*NMODE*(NMODE-1)*(NMODE-2)*
     1            (NMODE-3)/(4*3*2)
                  KJP4=JSIZE4*NMODE
                  IF(KJP4.GT.0)
     1            CALL MEMO(1,LJP4,KJP4,0,0,0,0,0,0,0,0)
                END IF
              END IF
            END IF
          END IF
          IF(ICI.GT.-1)GO TO 6600
          DO K=1,NMODE
            JCI1=NBAS(1)+1
            NMAX=MAXSUM(1)+1
            IF(JSIZE1.GT.0)
     1      CALL GETJP1(W(LJP1),JSIZE1,NMODE,K,JCI1,NMAX,ITOT1,
     2      W(LMXBAS))
            IF(ICI.GT.-2)GO TO 6601
            DO L=1,K-1
              JCI2=NBAS(2)+1
              NMAX=MAXSUM(2)+2
              IF(JSIZE2.GT.0)
     1        CALL GETJP2(W(LJP2),JSIZE2,NMODE,K,L,JCI2,NMAX,ITOT2,
     2        W(LMXBAS))
              IF(ICI.GT.-3)GO TO 6602
              DO N=1,L-1
                JCI3=NBAS(3)+1
                NMAX=MAXSUM(3)+3
                IF(JSIZE3.GT.0)
     1          CALL GETJP3(W(LJP3),JSIZE3,NMODE,K,L,N,JCI3,NMAX,ITOT3,
     2          W(LMXBAS))
                IF(ICI.GT.-4)GO TO 6603
                DO M=1,N-1
                  JCI4=NBAS(4)+1
                  NMAX=MAXSUM(4)+4
                  IF(JSIZE4.GT.0)
     1            CALL GETJP4(W(LJP4),JSIZE4,NMODE,K,L,N,M,JCI4,NMAX,
     2            ITOT4,W(LMXBAS))
                END DO
6603  CONTINUE
              END DO
6602  CONTINUE
            END DO
6601  CONTINUE
          END DO
          ISIZE=ITOT1+ITOT2+ITOT3+ITOT4
          JSIZE=ISIZE
          KIP=ISIZE*NMODE
          KJP=KIP
          CALL MEMO(1,LIP,KIP,0,0,0,0,0,0,0,0)
          IF(ICI.LT.0.AND.KJP1.GT.0)THEN
            CALL GETJP(W(LIP),ISIZE,W(LJP1),JSIZE1,NMODE,ITOT,ITOT1)
            CALL MEMO(-1,LJP1,KJP1,0,0,0,0,0,0,0,0)
          END IF
          IF(ICI.LT.-1.AND.KJP2.GT.0)THEN
            CALL GETJP(W(LIP),ISIZE,W(LJP2),JSIZE2,NMODE,ITOT,ITOT2)
            CALL MEMO(-1,LJP2,KJP2,0,0,0,0,0,0,0,0)
          END IF
          IF(ICI.LT.-2.AND.KJP3.GT.0)THEN
            CALL GETJP(W(LIP),ISIZE,W(LJP3),JSIZE3,NMODE,ITOT,ITOT3)
            CALL MEMO(-1,LJP3,KJP3,0,0,0,0,0,0,0,0)
          END IF
          IF(ICI.LT.-3.AND.KJP4.GT.0)THEN
            CALL GETJP(W(LIP),ISIZE,W(LJP4),JSIZE4,NMODE,ITOT,ITOT4)
            CALL MEMO(-1,LJP4,KJP4,0,0,0,0,0,0,0,0)
          END IF
          CALL MEMO(1,LJP,KJP,0,0,0,0,0,0,0,0)
          CALL GETIP(W(LIP),W(LJP),ISIZE,JSIZE,NMODE,0,0)
          WRITE(IOUT,395)(NTOT(K),K=1,NVSYM)
6600  CONTINUE
          IF(ITOT.EQ.0)THEN
            WRITE(IOUT,400)
            STOP 'SHOULD NOT OCCUR'
          END IF
          DO I=1,4
            NBAS(I)=JCIMAX-1
          END DO
        END IF
      END IF
      IF(MATSIZ.NE.0)STOP 'VCI MATRIX SIZE'
      NTOTMX=0
      DO NS=1,NVSYM
        IF(NTOT(NS).GT.NTOTMX)NTOTMX=NTOT(NS)
      END DO
      KEL21=J21
      IF(JMAX.NE.0.AND.ICI.LT.0)THEN
        KEL21=0
        DO KROT=1,J21,KSTEP
          KEL21=KEL21+1
        END DO
      END IF
      ISIZMX=NTOTMX
      CALL MEMO(1,LASSIG,ISIZMX*KEL21*NVSYM*3,0,0,0,0,0,0,0,0)
C**SKIP OVERLAP SECTION AND STATE-BY-STATE ALGORITHM IF VIRTUALS
      IF(ICI.LT.0)GO TO 4500
C**********************************************************************
C**********************************************************************
C**                                                             VSCF-CI
C**********************************************************************
C**********************************************************************
      CALL MEMO(-1,LIP,KIP,0,0,0,0,0,0,0,0)
      KS=NMODE
      CALL MEMO(1,LS,KS,0,0,0,0,0,0,0,0)
C**LOOP OVER SYMMETRIES
      ITIM=-1
      DO 6000 MS=1,MVSYM
C**********************************************************************
C**                                               LOOP ROUND SYMMETRIES
C**********************************************************************
      NS=MWSYM(MS)
      WRITE(IOUT,290)NS
      IF(NS.EQ.0)GO TO 6005
      ISIZE=NTOT(NS)
      WRITE(IOUT,285)ISIZE,CUT
      IF(ISIZE.EQ.0)GO TO 6005
      KIP=ISIZE*NMODE
      KOV=ISIZE*ISIZE
      CALL MEMO(2,LIP,KIP,LOV,KOV,0,0,0,0,0,0)
      CALL PUTJP(W(LJP),JSIZE,W(LIP),ISIZE,NMODE,NS)
      KCORIG=JTHIS
      DO 9500 KROT=1,J21
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
C**KA,KC FOR VSCF-CI
      KA=KROT/2
      KC=KCORIG-KA
C     IF(MOD(KA,2).NE.0)THEN
C       KC=KC+MOD(KROT,2)
C     ELSE
C       KC=KC+MOD(KROT+1,2)
C     END IF
      IF(MOD(KROT,2).EQ.0)KC=KC+1
C**ZEROISE OVERLAP MATRIX
      CALL DIAGZ(ISIZE,ISIZE,W(LOV),ISIZE,W(LOV),ICI)
C**SQUARE MATRIX NEEDED
      KXK=ISIZE*ISIZE
      CALL MEMO(1,LXK,KXK,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
      CALL DIAGZ(ISIZE,ISIZE,W(LXK),ISIZE,W(LXK),ICI)
      DO 2500 ISTATE=1,ISIZE
C********************************************
C**LOOP OVER LHS VSCF STATES TO BE CALCULATED
C********************************************
      REWIND 50
      INDEX=0
      DO IDUMMY=1,ICI
        IF(INDEX.EQ.0)THEN
C**READ LHS FUNCTIONS FOR THIS STATE
          K2=0
          DO K=1,NMODE
            CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
            CALL IN50(W(LHL+K2),I2)
            K2=K2+I2*3
          END DO
          CALL COMPARE(W(LISTAT),NSTAT,W(LIP),ISIZE,NMODE,IDUMMY,
     1    ISTATE,INDEX)
        END IF
      END DO
      DO 2000 JSTATE=ISTATE,ISIZE
C********************************************
C**LOOP OVER RHS VSCF STATES TO BE CALCULATED
C********************************************
      IF(ISTATE.NE.JSTATE)THEN
        REWIND 50
        JNDEX=0
        DO JDUMMY=1,ICI
          IF(JNDEX.EQ.0)THEN
C**READ RHS FUNCTIONS FOR THIS STATE
            K2=0
            DO K=1,NMODE
              CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
              CALL IN50(W(LHR+K2),I2)
              K2=K2+I2*3
            END DO
            CALL COMPARE(W(LISTAT),NSTAT,W(LIP),ISIZE,NMODE,JDUMMY,
     1      JSTATE,JNDEX)
          END IF
        END DO
      ELSE
C**COPY LHS FUNCTION
        DO K=1,KHLR
          W(LHR-1+K)=W(LHL-1+K)
        END DO
      END IF
C**************************************
C**GET OVERLAPS FOR THIS PAIR OF STATES
C**************************************
      K2=0
      DO K=1,NMODE
        CALL INTARR(W(LNBF),W(LMBF),K,I1,I2,I3)
        CALL OVERLP(W(LHL+K2),W(LHR+K2),I2,W(LS-1+K))
        CALL OVERLP(W(LHR+K2),W(LHL+K2),I2,W(LS-1+K))
        K2=K2+I2*3
      END DO
C**SET UP OVERLAP MATRIX
      CALL SCFOV(W(LS),NMODE,W(LOV),ISIZE,ISTATE,JSTATE)
      ITIM=ITIM+1
      IF(ICOUPL.GT.0)THEN
        IF(J21.GT.1)REWIND 61
        REWIND 71
        REWIND 81
      END IF
      IF(ICOUPL.GT.1)THEN
        IF(J21.GT.1)REWIND 62
        REWIND 72
        REWIND 82
      END IF
      IF(ICOUPL.GT.2)THEN
        IF(J21.GT.1)REWIND 63
        REWIND 73
        REWIND 83
      END IF
      IF(ICOUPL.GT.3)THEN
        IF(J21.GT.1)REWIND 64
        REWIND 74
        REWIND 84
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**'KINETIC ENERGY'
        CALL CIKE(NMODE,K,W(LHL+K2),W(LHR+K2),W(LXQ+K3),W(LS),W(LXK),
     1  ISIZE,IK2,W(LOMEGA+K1),ISTATE,JSTATE)
        IF(ICOUPL.EQ.0)GO TO 601
C**CORIOLIS AND POTENTIAL
        CALL CI1(NMODE,K,W(LHL+K2),W(LHR+K2),W(LXQ+K3),W(LS),W(LXK),
     1  ISIZE,IK2,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR)
     2  ,W(LXX),W(LX0),W(LXL),W(LXM),ISTATE,JSTATE,
     3  NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV1),W(LV1),W(LC1),W(LC1),
     4  W(LEJK1),W(LEJK1),J21,KROT)
        IF(ICOUPL.EQ.1)GO TO 601
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**CORIOLIS AND POTENTIAL
          CALL CI2(NMODE,K,L,W(LHL+K2),W(LHR+K2),W(LXQ+K3),W(LHL+L2),
     1    W(LHR+L2),W(LXQ+L3),IK2,IL2,W(LS),W(LXK),ISIZE,NATOM,
     2    W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR),W(LXX),W(LX0)
     3    ,W(LXL),W(LXM),ISTATE,JSTATE,
     4    NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV2),W(LV2),W(LC2),W(LC2),
     5    W(LEJK2),W(LEJK2),J21,KROT)
          IF(ICOUPL.EQ.2)GO TO 602
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**CORIOLIS AND POTENTIAL
            CALL CI3(NMODE,K,L,N,W(LHL+K2),W(LHR+K2),W(LXQ+K3),
     1      W(LHL+L2),W(LHR+L2),W(LXQ+L3),W(LHL+N2),W(LHR+N2),
     2      W(LXQ+N3),IK2,IL2,IN2,W(LS),W(LXK),ISIZE,NATOM,W(LQQ),
     3      W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR),W(LXX),
     4      W(LX0),W(LXL),W(LXM),ISTATE,JSTATE,
     5      NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV3),W(LV3),W(LC3),
     6      W(LC3),W(LEJK3),W(LEJK3),J21,KROT)
            IF(ICOUPL.EQ.3)GO TO 603
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**CORIOLIS AND POTENTIAL
              CALL CI4(NMODE,K,L,N,M,W(LHL+K2),W(LHR+K2),W(LXQ+K3),
     1        W(LHL+L2),W(LHR+L2),W(LXQ+L3),W(LHL+N2),W(LHR+N2),
     2        W(LXQ+N3),W(LHL+M2),W(LHR+M2),W(LXQ+M3),IK2,IL2,IN2,IM2,
     3        W(LS),W(LXK),ISIZE,NATOM,W(LQQ),
     4        W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LRR),W(LXX),
     5        W(LX0),W(LXL),W(LXM),ISTATE,JSTATE,
     6        NPOT,W(LIPOT),W(LJPOT),W(LCPOT),W(LV4),W(LV4),W(LC4),
     7        W(LC4),W(LEJK4),W(LEJK4),J21,KROT)
              IF(ICOUPL.EQ.4)GO TO 604
C**5-MODE COUPLING HERE IF NEEDED
604   CONTINUE
              M2=M2+IM2*3
              M3=M3+IM2
            END DO
603   CONTINUE
            N2=N2+IN2*3
            N3=N3+IN2
          END DO
602   CONTINUE
          L2=L2+IL2*3
          L3=L3+IL2
        END DO
601   CONTINUE
        K2=K2+IK2*3
        K3=K3+IK2
      END DO
2000  CONTINUE
2500  CONTINUE
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
      NVAL=ISCFCI
      NVEC=ISCFCI
      IF(ISIZE.LT.ISCFCI)THEN
        NVAL=ISIZE
        NVEC=ISIZE
      END IF
      WRITE(IOUT,260)
      KXA=ISIZE*ISIZE
      KWK=ISIZE
      KSUP4=ISIZE
      KEVAL=ISIZE
CCCC  IF(NS.EQ.1.AND.J21.EQ.1)EVL=0
CCCC  IF(NS.EQ.1.AND.J21.GT.1)EVL=EVLJ0
      IF(MS.EQ.1)EVL=EVLJ0
      CALL MEMO(4,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,LXA,KXA,0,0)
      WRITE(IOUT,*)'Calculating DIAG'
      CALL TIMIT(1)
      CALL DIAG(W(LXK),W(LXK),ISIZE,ISIZE,1,W(LSUP4),W(LEVAL),W(LWK),
     1NVAL,NVEC,W(LIP),ISIZE,NMODE,W(LOV),W(LXA),W(LASSIG),ISIZMX,
     2J21,NS)
      CALL TIMIT(3)
      CALL MEMO(-5,LXA,KXA,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL)
9500  CONTINUE
      CALL MEMO(-2,LIP,KIP,LOV,KOV,0,0,0,0,0,0)
6005  CONTINUE
6000  CONTINUE
      IF(ICOUPL.GT.0)THEN
        CALL MEMO(-1,LEJK1,KEJK1,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV1,KLV1,LC1,KLC1,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.1)THEN
        CALL MEMO(-1,LEJK2,KEJK2,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV2,KLV2,LC2,KLC2,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.2)THEN
        CALL MEMO(-1,LEJK3,KEJK3,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV3,KLV3,LC3,KLC3,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.3)THEN
        CALL MEMO(-1,LEJK4,KEJK4,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV4,KLV4,LC4,KLC4,0,0,0,0,0,0)
      END IF
      GO TO 5000
C**********************************************************************
C**********************************************************************
C**                                                    VCI - VIBRATIONS
C**********************************************************************
C**********************************************************************
4500  CONTINUE
      ITIM=-1
C**********************************************************************
C**                                  GET BASIC INTEGRALS FOR K-DIAGONAL
C**********************************************************************
C     IF(ICOUPL.GT.1)CALL MEMO(1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
      IF(ICOUPL.GE.0)THEN
        REWIND 21
        KXA1=ISIZE1*(ISIZE1+1)/2
      END IF
      IF(ICOUPL.GT.1)THEN
        REWIND 22
        KXA2=ISIZE2*(ISIZE2+1)/2
      END IF
      IF(ICOUPL.GT.2)THEN
        REWIND 23
        KXA3=ISIZE3*(ISIZE3+1)/2
      END IF
      IF(ICOUPL.GT.3)THEN
        REWIND 24
        KXA4=ISIZE4*(ISIZE4+1)/2
      END IF
      DO 3300 KROT=1,J21,KSTEP
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
      ITIM=ITIM+1
      IF(ICOUPL.GT.0)THEN
        IF(J21.GT.1)REWIND 61
        REWIND 71
        REWIND 81
      END IF
      IF(ICOUPL.GT.1)THEN
        IF(J21.GT.1)REWIND 62
        REWIND 72
        REWIND 82
      END IF
      IF(ICOUPL.GT.2)THEN
        IF(J21.GT.1)REWIND 63
        REWIND 73
        REWIND 83
      END IF
      IF(ICOUPL.GT.3)THEN
        IF(J21.GT.1)REWIND 64
        REWIND 74
        REWIND 84
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**ZEROISE MATRIX
        CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        CALL DIAGZ(ISIZE1,ISIZE1,W(LXA1),ISIZE1,W(LXA1),ICI)
C**'KINETIC ENERGY'
        CALL V0CIKE(NMODE,1,W(LH+K2),W(LXQ+K3),W(LXA1),
     1  IK3,IK2,W(LOMEGA+K1),W(LIP1),ISIZE1)
        IF(ICOUPL.EQ.0)THEN
          CALL MATOUT(W(LXA1),W(LXA1),ISIZE1,21)
          IF(ICOUPL.GE.0)CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          GO TO 3301
        END IF
C**CORIOLIS AND POTENTIAL
        CALL V0CI1(NMODE,1,W(LH+K2),W(LXQ+K3),W(LXA1),IK3,
     1  IK2,W(LIP1),ISIZE1,W(LV1),W(LV1),W(LC1),W(LC1),W(LEJK1),
     4  W(LEJK1),J21,KROT)
        CALL MATOUT(W(LXA1),W(LXA1),ISIZE1,21)
        IF(ICOUPL.GE.0)CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        IF(ICOUPL.EQ.1)GO TO 3301
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**ZEROISE MATRIX
          JCI=NBAS(2)+1
          KTEMP=JCI*JCI*3
          CALL MEMO(2,LXA2,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          CALL DIAGZ(ISIZE2,ISIZE2,W(LXA2),ISIZE2,W(LXA2),ICI)
C**CORIOLIS AND POTENTIAL
          CALL V0CI2(NMODE,1,2,W(LH+K2),W(LXQ+K3),W(LH+L2),
     1    W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA2),W(LIP2),ISIZE2,
     2    W(LTEMP),JCI,W(LV2),W(LV2),W(LC2),W(LC2),W(LEJK2),W(LEJK2),
     3    J21,KROT,W(LIP1),ISIZE1)
          CALL MATOUT(W(LXA2),W(LXA2),ISIZE2,22)
          CALL MEMO(-2,LXA2,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          IF(ICOUPL.EQ.2)GO TO 3302
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**ZEROISE MATRIX
            JCI=NBAS(3)+1
            KTEMP=JCI*JCI*JCI*JCI*3
            CALL MEMO(2,LXA3,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            CALL DIAGZ(ISIZE3,ISIZE3,W(LXA3),ISIZE3,W(LXA3),ICI)
C**CORIOLIS AND POTENTIAL
            CALL V0CI3(NMODE,1,2,3,W(LH+K2),W(LXQ+K3),
     1      W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,
     2      IN2,W(LXA3),W(LIP3),ISIZE3,W(LTEMP),JCI,W(LV3),W(LV3),
     3      W(LC3),W(LC3),W(LEJK3),W(LEJK3),J21,KROT,W(LIP2),ISIZE2)
            CALL MATOUT(W(LXA3),W(LXA3),ISIZE3,23)
            CALL MEMO(-2,LXA3,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            IF(ICOUPL.EQ.3)GO TO 3303
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**ZEROISE MATRIX
              JCI=NBAS(4)+1
              KTEMP=JCI*JCI*JCI*JCI*JCI*JCI*3
              CALL MEMO(2,LXA4,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
              CALL DIAGZ(ISIZE4,ISIZE4,W(LXA4),ISIZE4,W(LXA4),ICI)
C**CORIOLIS AND POTENTIAL
              CALL V0CI4(NMODE,1,2,3,4,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     2        IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA4),
     4        W(LIP4),ISIZE4,W(LTEMP),JCI,W(LV4),W(LV4),W(LC4),W(LC4),
     6        W(LEJK4),W(LEJK4),J21,KROT,W(LIP3),ISIZE3)
              CALL MATOUT(W(LXA4),W(LXA4),ISIZE4,24)
              CALL MEMO(-2,LXA4,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
              IF(ICOUPL.EQ.4)GO TO 3304
C**5-MODE COUPLING HERE IF NEEDED
3304  CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
3303  CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
3302  CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
3301  CONTINUE
        K2=K2+IK1
        K3=K3+IK2
      END DO
3300  CONTINUE
C     IF(ICOUPL.GT.1)CALL MEMO(-1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
      IF(ICOUPL.GT.0)THEN
        CALL MEMO(-1,LEJK1,KEJK1,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV1,KLV1,LC1,KLC1,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.1)THEN
        CALL MEMO(-1,LEJK2,KEJK2,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV2,KLV2,LC2,KLC2,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.2)THEN
        CALL MEMO(-1,LEJK3,KEJK3,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV3,KLV3,LC3,KLC3,0,0,0,0,0,0)
      END IF
      IF(ICOUPL.GT.3)THEN
        CALL MEMO(-1,LEJK4,KEJK4,0,0,0,0,0,0,0,0)
        CALL MEMO(-2,LV4,KLV4,LC4,KLC4,0,0,0,0,0,0)
      END IF
      ITIM=-1
C**********************************************************************
C**                                            SET UP K-DIAGONAL MATRIX
C**********************************************************************
      CALL MEMO(-1,LIP,KIP,0,0,0,0,0,0,0,0)
      ISIZMN=200000
      DO NS=1,NVSYM
        IF(NTOT(NS).LT.ISIZMN.AND.NTOT(NS).NE.0)ISIZMN=NTOT(NS)
      END DO
      NVAL=ISCFCI
      NVEC=ISCFCI
      IF(ISIZMN.LT.ISCFCI)THEN
        NVAL=ISIZMN
        NVEC=ISIZMN
      END IF
      ISIZMX=0
      DO NS=1,NVSYM
        IF(NTOT(NS).GT.ISIZMX)ISIZMX=NTOT(NS)
      END DO
      KIP=ISIZMX*NMODE
      CALL MEMO(3,LIP,KIP,LIPL,KIP,LIPR,KIP,0,0,0,0)
      IF(JTHIS.NE.0)THEN
        CALL MEMO(2,LCFS,ISIZMX*NVAL*KEL21*NVSYM,LEVCI,
     1  NVAL*KEL21*NVSYM,0,0,0,0,0,0)
      END IF
      IF(ISCFCI.GT.0.AND.IDUMP.NE.0)WRITE(60)NVSYM
C**SELECTIVE SYMMETRIES IF J=0
      MMSS=MVSYM
C**MUST DO ALL 'VIBRATIONAL' SYMMETRIES IF J>0
      IF(JTHIS.NE.0)MMSS=NVSYM
      DO 4400 MS=1,MMSS
C**********************************************************************
C**                                                LOOP OVER SYMMETRIES
C**********************************************************************
      IF(JTHIS.EQ.0)THEN
        NS=MWSYM(MS)
      ELSE
        NS=MS
      END IF
      WRITE(IOUT,290)NS
      IF(NS.EQ.0)GO TO 4401
      ITIM=ITIM+1
      ISIZE=NTOT(NS)
      WRITE(IOUT,285)ISIZE,CUT
      CALL FLUSH(IOUT)
      IF(ISIZE.EQ.0)GO TO 4401
      CALL PUTJP(W(LJP),JSIZE,W(LIP),ISIZMX,NMODE,NS)
C**DUMP BASIS TO (60)
      IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
        CALL OUTBAS(ISIZE,W(LIP),ISIZMX,NMODE)
      END IF
      IF(ICOUPL.GE.0)THEN
        REWIND 21
      END IF
      IF(ICOUPL.GT.1)THEN
        REWIND 22
      END IF
      IF(ICOUPL.GT.2)THEN
        REWIND 23
      END IF
      IF(ICOUPL.GT.3)THEN
        REWIND 24
      END IF
      KCORIG=JTHIS
      KEL=0
      DO 9000 KROT=1,J21,KSTEP
      KEL=KEL+1
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
C**KA,KC FOR VCI
      KA=KROT/2
      KC=KCORIG-KA
      IF(KMAX.GT.0)THEN
        IF(MOD(KA,2).NE.0)THEN
          KC=KC+MOD(KROT,2)
        ELSE
C         KC=KC+MOD(KROT+1,2)
          IF(KA.NE.0)KC=KC+MOD(KROT,2)
        END IF
      ELSE
        IF(MOD(KROT,2).EQ.0)KC=KC+1
      END IF
C**ISTART, IEND ARE START AND END COLUMNS
      ISTART=1
      IF(LANCZ)THEN
        CALL MEMO(3,LWK,ISIZE,LYK,ISIZE,LZK,ISIZE,0,0,0,0)
        KLAN=LANMAX*(LANMAX+1)/2
        LSIZE=ISIZE
      ELSE
        KXA=ISIZE*(ISIZE+1)/2
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGZ(ISIZE,ISIZE,W(LXA),ISIZE,W(LXA),ICI)
        IEND=ISIZE
      END IF
      ISKIP=0
6666  CONTINUE
      IF(LANCZ)THEN
        IF(ISKIP.NE.0)THEN
          IF(ICOUPL.GE.0)THEN
            REWIND 21
          END IF
          IF(ICOUPL.GT.1)THEN
            REWIND 22
          END IF
          IF(ICOUPL.GT.2)THEN
            REWIND 23
          END IF
          IF(ICOUPL.GT.3)THEN
            REWIND 24
          END IF
C**RE-POSITION DISCS
          DO KKK=1,KROT-KSTEP,KSTEP
            DO K=1,NMODE
              IF(ICOUPL.GT.0)THEN
                CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
                CALL MATIN(W(LXA1),W(LXA1),ISIZE1,21)
                CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
              END IF
              DO L=1,K-1
                IF(ICOUPL.GT.1)THEN
                  CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
                  CALL MATIN(W(LXA2),W(LXA2),ISIZE2,22)
                  CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
                END IF
                DO N=1,L-1
                  IF(ICOUPL.GT.2)THEN
                    CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
                    CALL MATIN(W(LXA3),W(LXA3),ISIZE3,23)
                    CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
                  END IF
                  DO M=1,N-1
                    IF(ICOUPL.GT.3.AND.ICOUPL.NE.4)THEN
                      CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                      CALL MATIN(W(LXA4),W(LXA4),ISIZE4,24)
                      CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
                    END IF
4004  CONTINUE
                  END DO
4003  CONTINUE
                END DO
4002  CONTINUE
              END DO
4001  CONTINUE
            END DO
          END DO
        END IF
        IEND=ISTART-1
        KXA=0
        DO I=ISTART,ISIZE
          KXA=KXA+LSIZE
          IF(KXA.GT.KLAN)THEN
            KXA=KXA-LSIZE
            GO TO 6665
          END IF
          LSIZE=LSIZE-1
          IEND=IEND+1
        END DO
6665    CONTINUE
        IF(IEND.EQ.0)STOP 'LANMAX TOO SMALL'
        WRITE(IOUT,*)'ISTART,IEND,KXA',ISTART,IEND,KXA
        CALL FLUSH(IOUT)
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGL(KXA,W(LXA))
      END IF
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**'KINETIC ENERGY'
        IF(ICOUPL.EQ.0)THEN
C**READ INTO MATRIX
          CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MATIN(W(LXA1),W(LXA1),ISIZE1,21)
          CALL VCIKE(NMODE,K,W(LH+K2),W(LXQ+K3),W(LXA),W(LXA1),
     1    IK3,IK2,W(LOMEGA+K1),W(LIP),ISIZMX,ISIZE,W(LIP1),ISIZE1,
     2    ISTART,IEND)
          IF(ICOUPL.GE.0)CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          GO TO 701
        END IF
C**CORIOLIS AND POTENTIAL
C**READ INTO MATRIX
        CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        CALL MATIN(W(LXA1),W(LXA1),ISIZE1,21)
        CALL VCI1(NMODE,K,W(LH+K2),W(LXQ+K3),W(LXA),W(LXA1),IK3,
     1  IK2,W(LIP),ISIZMX,ISIZE,W(LIP1),ISIZE1,ISTART,IEND)
        IF(ICOUPL.GE.0)CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        IF(ICOUPL.EQ.1)GO TO 701
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**CORIOLIS AND POTENTIAL
C**READ INTO MATRIX
          JCI=NBAS(2)+1
          CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
          CALL MATIN(W(LXA2),W(LXA2),ISIZE2,22)
          CALL VCI2(NMODE,K,L,W(LH+K2),W(LXQ+K3),W(LH+L2),
     1    W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA),W(LXA2),
     3    W(LIP),ISIZMX,ISIZE,W(LIP2),ISIZE2,JCI,ISTART,IEND)
          IF(ICOUPL.GT.1)CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
          IF(ICOUPL.EQ.2)GO TO 702
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**CORIOLIS AND POTENTIAL
C**READ INTO MATRIX
            JCI=NBAS(3)+1
            CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
            CALL MATIN(W(LXA3),W(LXA3),ISIZE3,23)
            CALL VCI3(NMODE,K,L,N,W(LH+K2),W(LXQ+K3),
     1      W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,IL3,IL2,IN3,
     2      IN2,W(LXA),W(LXA3),
     4      W(LIP),ISIZMX,ISIZE,W(LIP3),ISIZE3,JCI,ISTART,IEND)
            IF(ICOUPL.GT.2)CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
            IF(ICOUPL.EQ.3)GO TO 703
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**CORIOLIS AND POTENTIAL
C**READ INTO MATRIX
              JCI=NBAS(4)+1
              CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
              CALL MATIN(W(LXA4),W(LXA4),ISIZE4,24)
              CALL VCI4(NMODE,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),W(LXQ+M3),
     2        IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA),W(LXA4),
     4        W(LIP),ISIZMX,ISIZE,W(LIP4),ISIZE4,JCI,ISTART,IEND)
              IF(ICOUPL.GT.3)CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
              IF(ICOUPL.EQ.4)GO TO 704
C**5-MODE COUPLING HERE IF NEEDED
704   CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
703   CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
702   CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
701   CONTINUE
        K2=K2+IK1
        K3=K3+IK2
      END DO
      ISKIP=0
      IF(LANCZ)THEN
        CALL LANOUT(W(LXA),W(LWK),W(LWK),W(LYK),W(LZK),W(LZK),ISIZE,
     1  ISTART,IEND)
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        IF(IEND.NE.ISIZE)THEN
          ISTART=IEND+1
          ISKIP=1
        END IF
      END IF
      IF(ISKIP.NE.0)GO TO 6666
      IF(LANCZ)CALL MEMO(-1,LZK,ISIZE,0,0,0,0,0,0,0,0)
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
      WRITE(IOUT,260)
      KXK=ISIZE*NVAL
      KWK=ISIZE
      KSUP4=5*ISIZE
      KEVAL=NVAL
CCCC  IF(NS.EQ.1.AND.JTHIS.EQ.0)EVL=0
CCCC  IF(NS.EQ.1.AND.JTHIS.NE.0)EVL=EVLJ0
      IF(MS.EQ.1)EVL=EVLJ0
      CALL TIMIT(1)
      IF(LANCZ)THEN
        LGIV=(IGIV.NE.0)
        KXKL=NVAL*NCYCLE*NVAL*NCYCLE
        IF(LGIV)KXKL=NVAL*NCYCLE*NVAL
        IF(LANZA)THEN
          CALL MEMO(3,LXK,ISIZE,LVK,ISIZE,LZK,ISIZE,0,0,0,0)
          CALL MEMO(5,LSK,ISIZE,LEN,NVAL,LEL,NVAL,LNV,NVAL,LWKL,
     1    ISIZE)
          CALL MEMO(4,LWRK,NVAL*NCYCLE,LSUP4,5*NVAL*NCYCLE,LEVAL,
     1    NVAL*NCYCLE,LXKL,KXKL,0,0)
          KXA=NVAL*NCYCLE*(NVAL*NCYCLE+1)/2
          IF(LGIV)CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
          CALL MEMO(1,LXKLC,KXA,0,0,0,0,0,0,0,0)
          WRITE(IOUT,*)'Calculating LANCZOS'
          IF(LGIV)THEN
            CALL LANCZA(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVAL,W(LIP),ISIZMX,NMODE,W(LASSIG),ISIZMX,KEL21,NS,
     2      W(LEN),W(LEL),W(LNV),W(LXA),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          ELSE
            CALL LANCZA(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVAL,W(LIP),ISIZMX,NMODE,W(LASSIG),ISIZMX,KEL21,NS,
     2      W(LEN),W(LEL),W(LNV),W(LXKL),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          END IF
          IF(JTHIS.NE.0)THEN
C**SAVE VCI ENERGIES AND COEFFICIENTS FOR Ka>0
            CALL ENERCF(W(LCFS),ISIZMX,W(LEVCI),ISIZE,NVAL,KEL21,KEL,
     1      NS,W(LXK),W(LXK),W(LEL))
          END IF
C**DUMP CI COEFFICIENTS TO (60)
          IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
            IF(JTHIS.EQ.0)THEN
              CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,JDUMP,IDUMP)
            ELSE
CCCC          CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,NVAL,NVAL)
            END IF
          END IF
          CALL MEMO(-5,LXK,ISIZE,LVK,ISIZE,LZK,ISIZE,LWK,ISIZE,LYK,
     1    ISIZE)
          CALL MEMO(-5,LSK,ISIZE,LEN,NVAL,LEL,NVAL,LNV,NVAL,LXKL,
     1    KXKL)
          CALL MEMO(-5,LXKLC,KXA,LWRK,NVAL*NCYCLE,
     1    LSUP4,5*NVAL*NCYCLE,LEVAL,NVAL*NCYCLE,LWKL,ISIZE)
        ELSE
          CALL MEMO(3,LXK,ISIZE*NVAL,LVK,ISIZE*NVAL,LZK,ISIZE*NVAL,0,0,
     1    0,0)
          CALL MEMO(5,LSK,ISIZE*NVAL,LEN,NVAL,LEL,NVAL,LNV,NVAL,LWKL,
     1    ISIZE*NVAL)
          CALL MEMO(4,LWRK,NVAL*NCYCLE,LSUP4,5*NVAL*NCYCLE,LEVAL,
     1    NVAL*NCYCLE,LXKL,KXKL,0,0)
          KXA=NVAL*NCYCLE*(NVAL*NCYCLE+1)/2
          IF(LGIV)CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
          CALL MEMO(1,LXKLC,KXA,0,0,0,0,0,0,0,0)
          WRITE(IOUT,*)'Calculating LANCZOS'
          IF(LGIV)THEN
            CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVAL,W(LIP),ISIZMX,NMODE,W(LASSIG),ISIZMX,KEL21,NS,
     2      W(LEN),W(LEL),W(LNV),W(LXA),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          ELSE
            CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVAL,W(LIP),ISIZMX,NMODE,W(LASSIG),ISIZMX,KEL21,NS,
     2      W(LEN),W(LEL),W(LNV),W(LXKL),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          END IF
          IF(JTHIS.NE.0)THEN
C**SAVE VCI ENERGIES AND COEFFICIENTS FOR Ka>0
            CALL ENERCF(W(LCFS),ISIZMX,W(LEVCI),ISIZE,NVAL,KEL21,KEL,
     1      NS,W(LXK),W(LXK),W(LEL))
          END IF
C**DUMP CI COEFFICIENTS TO (60)
          IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
            IF(JTHIS.EQ.0)THEN
              CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,JDUMP,IDUMP)
            ELSE
CCCC          CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,NVAL,NVAL)
            END IF
          END IF
          CALL MEMO(-5,LXK,ISIZE*NVAL,LVK,ISIZE*NVAL,LZK,ISIZE*NVAL,
     1    LWK,ISIZE,LYK,ISIZE)
          CALL MEMO(-5,LSK,ISIZE*NVAL,LEN,NVAL,LEL,NVAL,LNV,NVAL,LXKL,
     1    KXKL)
          CALL MEMO(-5,LXKLC,KXA,LWRK,NVAL*NCYCLE,
     1    LSUP4,5*NVAL*NCYCLE,LEVAL,NVAL*NCYCLE,LWKL,ISIZE*NVAL)
        END IF
        IF(LGIV)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        LGIV=.FALSE.
      ELSE
        LGIV=.TRUE.
        CALL MEMO(4,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,0,0)
        WRITE(IOUT,*)'Calculating DIAG'
        CALL DIAG(W(LXA),W(LXK),ISIZE,ISIZE,1,W(LSUP4),W(LEVAL),W(LWK),
     1  NVAL,NVEC,W(LIP),ISIZMX,NMODE,W(LXK),W(LXK),W(LASSIG),ISIZMX,
     2  KEL21,NS)
        IF(JTHIS.NE.0)THEN
C**SAVE VCI ENERGIES AND COEFFICIENTS FOR Ka>0
          CALL ENERCF(W(LCFS),ISIZMX,W(LEVCI),ISIZE,NVAL,KEL21,KEL,NS,
     1    W(LXK),W(LXK),W(LWK))
        END IF
C**DUMP CI COEFFICIENTS TO (60)
        IF(ISCFCI.GT.0.AND.IDUMP.NE.0)THEN
          IF(JTHIS.EQ.0)THEN
            CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,JDUMP,IDUMP)
          ELSE
CCCC        CALL OUTCIV(ISIZE,W(LXK),W(LXK),W(LNDUMP),NS,NVAL,NVAL)
          END IF
        END IF
        CALL MEMO(-4,LXK,KXK,LWK,KWK,LSUP4,KSUP4,LEVAL,KEVAL,0,0)
        LGIV=.FALSE.
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
      END IF
      CALL TIMIT(3)
      CALL FLUSH(IOUT)
9000  CONTINUE
4401  CONTINUE
4400  CONTINUE
      IF(JTHIS.EQ.0)GO TO 4444
      IF(KMAX.LT.0)GO TO 4444
C**********************************************************************
C**********************************************************************
C**                                                     VCI - ROTATIONS
C**********************************************************************
C**********************************************************************
      CALL VIBROT(W)
4444  CONTINUE
      CALL FLUSH(IOUT)
      CALL MEMO(-1,LIP,KIP,0,0,0,0,0,0,0,0)
      IF(JTHIS.GT.0)
     1CALL MEMO(-2,LCFS,ISIZMX*NVAL*KEL21*NVSYM,LEVCI,NVAL*KEL21*NVSYM,
     20,0,0,0,0,0)
C*************************************************************
C*************************************************************
5000  CONTINUE
      WRITE(IOUT,240)EVL
      WRITE(IOUT,225)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE MAXIN(MAXBAS,NMODE,J,ICI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,J)
      COMMON/BASIS/NBAS(4),MAXSUM(4)
      COMMON/FILASS/IOUT,INP
105   FORMAT(1X,I1,'-MODE STATES: ',/,1X,50I3)
      DO K=1,J
        READ(INP,*)(MAXBAS(I,K),I=1,NMODE)
        IF(ICI.LT.0)WRITE(IOUT,105)K,(MAXBAS(I,K),I=1,NMODE)
C**CONVERT TO NO. FUNCTIONS
        IF(IABS(ICI).NE.100)THEN
          DO I=1,NMODE
            IF(MAXBAS(I,K).LT.0)MAXBAS(I,K)=0
            IF(MAXBAS(I,K).GT.NBAS(K))NBAS(K)=MAXBAS(I,K)
            MAXBAS(I,K)=MAXBAS(I,K)+1
          END DO
        END IF
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE VIBROT(W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 TITLE(80)
      LOGICAL LGIV,LINEAR,LANCZ,LANZA,LANZB,TRIAT
C**************************************************ASSIGN TOTAL STORAGE
      DIMENSION W(1)
      COMMON/CMEMO/NADD
      COMMON/CADDR/LIPOT,LJPOT,LCPOT,LOMEGA,LISTAT,LH,LXK,LXQ,LXW,LNBF,
     1LMBF,LWK,LEVAL,LXM,LX0,LXL,LXZ,LXX,LR0,LRR,
     2LAB,LB,LAA,LBB,LQQ,LESCF,LXA,LS,LHL,LHR,
     3LSUP4,LOV,LV1,LV2,LV3,LV4,LIP,LJP,LTEMP,LC1,
     4LC2,LC3,LC4,LXK0,LXL0,LXN0,LXM0,LEJK1,LEJK2,LEJK3,
     5LEJK4,LW21,LSS,LSSX,LX21,LE21,LJSTAT,LKSTAT,LESTAT,LWSTAT,
     6LVM1,LVM2,LVM3,LVM4,LCFS,LEVCI,LASSIG,LYL,LYZ,LY0,
     7LYK,LSK,LWRK,LVK,LZK,LXKL,LXKLC,LEN,LEL,LNV,
     8LWKL,LIP1,LJP1,LIP2,LJP2,LIP3,LJP3,LIP4,LJP4,LXA1,
     9LXA2,LXA3,LXA4,LIPL,LIPR,LXV,LQM,LMXBAS,LNDUMP,LNVF
C**************************************************ASSIGN TOTAL STORAGE
      COMMON/TITLE/TITLE
      COMMON/HERM/IHERM
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/LANCZO/LANCZ,LANZA,LANZB
      COMMON/CYCLES/NCYCLE
      COMMON/TRIATO/TRIAT
      COMMON/LANTOL/TOLLAN
      COMMON/MAXLAN/LANMAX
      COMMON/GIVEN/LGIV,IGIV
      COMMON/TYPE/LINEAR
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/TOLS/TOL,EPS
      COMMON/EVL/EVL,CUT
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/VCIMAX/NMAX
      COMMON/ROTS/JMAX,KMAX,J21,KEL21,KEL
      COMMON/ESTATE/IORDER
      COMMON/JKAKC/JTHIS,KA,KC
      COMMON/AXES/MX(3)
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/BASIS/NBAS(4),MAXSUM(4)
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/SYMMS/MVSYM,MWSYM(10)
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/MAXPT/MBFMAX,MBFMX1,MBFMX2,MBFMX3,MBFMX4,MBFMIN
      COMMON/DISC/IDISC
      COMMON/MODES/NMODE,NATOM
      COMMON/SIZES/KTEMP,ISIZE1,ISIZE2,ISIZE3,ISIZE4,ISIZE,JSIZE,
     1ISIZMX
      COMMON/MATRIX/NVAL,NVALR,KSTEP,KSIGN
C**********************************************************************
100   FORMAT(80A1)
200   FORMAT(//,1X,80A1)
205   FORMAT(/,1X,'NUMBER OF MODES = ',I4,/,
     1         1X,'NUMBER OF SCF STATES = ',I4,/,
     2         1X,'NUMBER POTENTIAL TERMS = ',I4,/,
     3         1X,'PRINT LEVEL = ',I4,/,
     4         1X,'TOLERANCE FOR SCF = ',D20.12,/,
     5         1X,'COUPLING OF ',I4,' MODES',/)
210   FORMAT(/,1X,'IDISC = ',I4,/,
     1         1X,'(WRITE POTENTIAL AND CORIOLIS INFO. TO DISC (0)',/,
     2         1X,'POTENTIAL AND CORIOLIS DISCS ALREADY EXIST (1))',/)
215   FORMAT(/,1X,'SCF CALCULATION ONLY',/)
220   FORMAT(//,1X,'START SCF CALCULATION',/)
225   FORMAT(//,1X,'SCF PLUS CI CALCULATION')
230   FORMAT(/,1X,'NUMBER OF CI ENERGIES REQUIRED = ',I4,/)
235   FORMAT(//,1X,'START CI CALCULATION',/)
240   FORMAT(/,1X,'ZERO POINT ENERGY = ',F10.2,/)
245   FORMAT(//,1X,'OVERLAPS OF ORIGINAL SCF STATE FUNCTIONS',/)
250   FORMAT(/,1X,'OVERLAPS WITH STATE ',I4)
255   FORMAT(//,1X,'TEST OF SCHMIDT ORTHOGONALISATION',/)
260   FORMAT(//,1X,'FINAL VIBRATIONAL (K-DIAGONAL) CI ENERGIES',/)
265   FORMAT(//,1X,'NORMAL COORDINATE POTENTIAL',/)
270   FORMAT(//,1X,'INTERNAL COORDINATE POTENTIAL',/)
275   FORMAT(//,1X,'CI MATRIX INVOLVES ',I4,' SCF FUNCTIONS',/)
280   FORMAT(//,1X,'CI MATRIX INVOLVES ',I4,' VIRTUAL GS FUNCTIONS',/,
     1          1X,'WITH SUM OF QUANTA < ',I4,/)
285   FORMAT(//,1X,'SIZE OF VIBRATIONAL CI MATRIX IS ',I5,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
290   FORMAT(//,1X,'VIBRATIONAL SYMMETRY ',I2)
295   FORMAT(//,1X,'NUMBER OF VIBRATIONAL SYMMETRIES = ',I3,/,
     1          1X,'NUMBER OF MODE SYMMETRIES = ',I3,/)
300   FORMAT(1X,'MODE SYMMETRY ',I3,'   MODES ',20I3)
305   FORMAT(//,1X,'TOTAL ANGULAR MOMENTUM J = ',I3,/)
310   FORMAT(/,1X,'J = 0 SCF CYCLE',/)
315   FORMAT(/,1X,'J = ',I3,' SCF CYCLE',/)
320   FORMAT('*************************')
325   FORMAT(50(1H*))
350   FORMAT(/,1X,'LINK TO PROGRAM NORMALS',/)
355   FORMAT(/,1X,'ILLEGAL USE OF INORM',/)
360   FORMAT(/,1X,'NUMBER OF LANCZOS CYCLES = ',I3)
365   FORMAT(/,1X,'LANCZOS TOLERANCE = ',F10.6)
370   FORMAT(/,1X,'LANCZOS (HALF) MATRIX MAXIMUM ORDER = ',I6,/)
375   FORMAT(//,1X,'ROTATIONAL SYMMETRY ',I2)
380   FORMAT(//,1X,'SIZE OF RO-VIBRATIONAL CI MATRIX IS ',I5,/,
     1          1X,'CUT-OFF ENERGY IS ',F10.2,/)
385   FORMAT(//,1X,'FINAL RO-VIBRATIONAL CI ENERGIES',/)
390   FORMAT(1X,'BLOCK ',I2,' VIBRATIONAL SYMMETRY ',I2)
395   FORMAT(//,1X,'SIZES OF CI SYMMETRY BLOCKS: ',/,1X,10I7,/)
400   FORMAT(//,1X,'SHOULD NOT OCCUR',/)
405   FORMAT(//,1X,'MINIMIZATION OF POTENTIAL',/)
410   FORMAT(//,1X,'DEGREE OF POLYNOMIAL: ',I3,/)
C**********************************************************************
C**********************************************************************
C**                                                     VCI - ROTATIONS
C**********************************************************************
C**********************************************************************
      ITIM=-1
C**********************************************************************
      IF(ICOUPL.EQ.0)GO TO 3000
      IF(JCOUPL.GT.0)THEN
        KVM1=MBFMX1*6
      ELSE
        KVM1=(1+MBFMX1*6)/2
      END IF
      CALL MEMO(1,LVM1,KVM1,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**WRITE CORIOLIS GRIDS TO DISC (91)
          CALL DUMCR1(W(LXQ+K3),IK2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     1    W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,W(LC1),
     2    W(LC1),W(LVM1),W(LVM1),1)
3001      CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.1)GO TO 3000
      IF(JCOUPL.GT.0)THEN
        KVM2=MBFMX2*12
      ELSE
        KVM2=(1+MBFMX2*12)/2
      END IF
      CALL MEMO(1,LVM2,KVM2,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**WRITE CORIOLIS GRIDS TO DISC (92)
            CALL DUMCR2(W(LXQ+K3),W(LXQ+L3),IK2,IL2,NMODE,NATOM,
     1      W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),W(LXX),W(LX0),
     2      W(LXL),W(LXM),K,L,W(LC2),W(LC2),W(LVM2),W(LVM2),1)
3022        CONTINUE
            L3=L3+IL2
          END DO
3002      CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.2)GO TO 3000
      IF(JCOUPL.GT.0)THEN
        KVM3=MBFMX3*15
      ELSE
        KVM3=(1+MBFMX3*15)/2
      END IF
      CALL MEMO(1,LVM3,KVM3,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**WRITE CORIOLIS GRIDS TO DISC (93)
              CALL DUMCR3(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),IK2,IL2,IN2,
     1        NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),W(LB),W(LAA),W(LBB),
     2        W(LXX),W(LX0),W(LXL),W(LXM),K,L,N,W(LC3),W(LC3),W(LVM3),
     3        W(LVM3),1)
3333          CONTINUE
              N3=N3+IN2
            END DO
3033        CONTINUE
            L3=L3+IL2
          END DO
3003      CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.3)GO TO 3000
      IF(JCOUPL.GT.0)THEN
        KVM4=MBFMX4*18
      ELSE
        KVM4=(1+MBFMX4*18)/2
      END IF
      CALL MEMO(1,LVM4,KVM4,0,0,0,0,0,0,0,0)
      IF(IDISC.EQ.0)THEN
        K3=0
        DO K=1,NMODE
          CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
          L3=0
          DO L=1,K-1
            CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
            N3=0
            DO N=1,L-1
              CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
              M3=0
              DO M=1,N-1
                CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**WRITE CORIOLIS GRIDS TO DISC (94)
                CALL DUMCR4(W(LXQ+K3),W(LXQ+L3),W(LXQ+N3),W(LXQ+M3),
     1          IK2,IL2,IN2,IM2,NMODE,NATOM,W(LQQ),W(LXZ),W(LAB),
     2          W(LB),W(LAA),W(LBB),W(LXX),W(LX0),W(LXL),W(LXM),K,L,
     3          N,M,W(LC4),W(LC4),W(LVM4),W(LVM4),1)
34444           CONTINUE
                M3=M3+IM2
              END DO
30444         CONTINUE
              N3=N3+IN2
            END DO
30044       CONTINUE
            L3=L3+IL2
          END DO
30004     CONTINUE
          K3=K3+IK2
        END DO
      END IF
C**********************************************************************
      IF(ICOUPL.EQ.4)GO TO 3000
C**5-MODE AND HIGHER
3000  CONTINUE
C     IF(ICOUPL.GT.1)CALL MEMO(1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
      IF(ICOUPL.GE.0)THEN
        REWIND 21
        KXA1=ISIZE1*(ISIZE1+1)/2
      END IF
      IF(ICOUPL.GT.1)THEN
        REWIND 22
        KXA2=ISIZE2*(ISIZE2+1)/2
      END IF
      IF(ICOUPL.GT.2)THEN
        REWIND 23
        KXA3=ISIZE3*(ISIZE3+1)/2
      END IF
      IF(ICOUPL.GT.3)THEN
        REWIND 24
        KXA4=ISIZE4*(ISIZE4+1)/2
      END IF
      DO 4444 IABC=1,9
C**********************************************************************
C**   GET BASIC INTEGRALS FOR FILES  (31)-(39) IF J>0
C**********************************************************************
      ITIM=ITIM+1
      IF(ICOUPL.GT.0)REWIND 91
      IF(ICOUPL.GT.1)REWIND 92
      IF(ICOUPL.GT.2)REWIND 93
      IF(ICOUPL.GT.3)REWIND 94
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        K1=K-1
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
C**CORIOLIS
        IF(IABC.LT.7)THEN
C**ZEROISE MATRIX
          CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL DIAGZ(ISIZE1,ISIZE1,W(LXA1),ISIZE1,W(LXA1),ICI)
          CALL V0MI1(NMODE,1,W(LH+K2),W(LXQ+K3),W(LXA1),IK3,
     1    IK2,W(LIP1),ISIZE1,W(LVM1),W(LVM1),J21,IABC)
          CALL MATOUT(W(LXA1),W(LXA1),ISIZE1,21)
          IF(ICOUPL.GE.0)CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        END IF
        IF(ICOUPL.EQ.1)GO TO 7701
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**ZEROISE MATRIX
          JCI=NBAS(2)+1
          KTEMP=JCI*JCI*3
          CALL MEMO(2,LXA2,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          CALL DIAGZ(ISIZE2,ISIZE2,W(LXA2),ISIZE2,W(LXA2),ICI)
C**CORIOLIS
          CALL V0MI2(NMODE,1,2,W(LH+K2),W(LXQ+K3),W(LH+L2),
     1    W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA2),
     2    W(LIP2),ISIZE2,W(LTEMP),JCI,W(LVM2),W(LVM2),J21,IABC,
     3    W(LIP1),ISIZE1)
          CALL MATOUT(W(LXA2),W(LXA2),ISIZE2,22)
          CALL MEMO(-2,LXA2,KXA2,LTEMP,KTEMP,0,0,0,0,0,0)
          IF(ICOUPL.EQ.2)GO TO 7702
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**ZEROISE MATRIX
            JCI=NBAS(3)+1
            KTEMP=JCI*JCI*JCI*JCI*3
            CALL MEMO(2,LXA3,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            CALL DIAGZ(ISIZE3,ISIZE3,W(LXA3),ISIZE3,W(LXA3),ICI)
C**CORIOLIS
            CALL V0MI3(NMODE,1,2,3,W(LH+K2),W(LXQ+K3),
     1      W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,IL3,IL2,
     2      IN3,IN2,W(LXA3),W(LIP3),ISIZE3,W(LTEMP),JCI,W(LVM3),
     3      W(LVM3),J21,IABC,W(LIP2),ISIZE2)
            CALL MATOUT(W(LXA3),W(LXA3),ISIZE3,23)
            CALL MEMO(-2,LXA3,KXA3,LTEMP,KTEMP,0,0,0,0,0,0)
            IF(ICOUPL.EQ.3)GO TO 7703
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**ZEROISE MATRIX
              JCI=NBAS(4)+1
              KTEMP=JCI*JCI*JCI*JCI*JCI*JCI*3
              CALL MEMO(2,LXA4,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
              CALL DIAGZ(ISIZE4,ISIZE4,W(LXA4),ISIZE4,W(LXA4),ICI)
C**CORIOLIS
              CALL V0MI4(NMODE,1,2,3,4,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2        W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA4),
     3        W(LIP4),ISIZE4,W(LTEMP),JCI,W(LVM4),W(LVM4),J21,IABC,
     4        W(LIP3),ISIZE3)
              CALL MATOUT(W(LXA4),W(LXA4),ISIZE4,24)
              CALL MEMO(-2,LXA4,KXA4,LTEMP,KTEMP,0,0,0,0,0,0)
              IF(ICOUPL.EQ.4)GO TO 7704
C**5-MODE COUPLING HERE IF NEEDED
7704  CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
7703  CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
7702  CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
7701  CONTINUE
        K2=K2+IK1
        K3=K3+IK2
      END DO
4444  CONTINUE
C     IF(ICOUPL.GT.1)CALL MEMO(-1,LTEMP,KTEMP,0,0,0,0,0,0,0,0)
      IF(ICOUPL.GT.0)CALL MEMO(-1,LVM1,KVM1,0,0,0,0,0,0,0,0)
      IF(ICOUPL.GT.1)CALL MEMO(-1,LVM2,KVM2,0,0,0,0,0,0,0,0)
      IF(ICOUPL.GT.2)CALL MEMO(-1,LVM3,KVM3,0,0,0,0,0,0,0,0)
      IF(ICOUPL.GT.3)CALL MEMO(-1,LVM4,KVM4,0,0,0,0,0,0,0,0)
      ITIM=-1
C**********************************************************************
C**   SET UP FILES  (31)-(39) FOR COMMON INTEGRALS
C**********************************************************************
      REWIND 31
      REWIND 32
      REWIND 33
      REWIND 34
      REWIND 35
      REWIND 36
      REWIND 37
      REWIND 38
      REWIND 39
      DO 5500 NSL=1,NVSYM
C**********************************************************************
C**                                          LOOP OVER SYMMETRIES (LHS)
C**********************************************************************
      ISIZEL=NTOT(NSL)
      IF(ISIZEL.EQ.0)GO TO 5501
      CALL PUTJP(W(LJP),JSIZE,W(LIPL),ISIZMX,NMODE,NSL)
      DO 5502 NSR=NSL,NVSYM
C**********************************************************************
C**                                          LOOP OVER SYMMETRIES (RHS)
C**********************************************************************
      ISIZER=NTOT(NSR)
      IF(ISIZER.EQ.0)GO TO 5503
      CALL PUTJP(W(LJP),JSIZE,W(LIPR),ISIZMX,NMODE,NSR)
      IF(ICOUPL.GE.0)THEN
        REWIND 21
      END IF
      IF(ICOUPL.GT.1)THEN
        REWIND 22
      END IF
      IF(ICOUPL.GT.2)THEN
        REWIND 23
      END IF
      IF(ICOUPL.GT.3)THEN
        REWIND 24
      END IF
      DO 9999 IABC=1,9
      ITIM=ITIM+1
      KXA=ISIZEL*ISIZER
      CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
      CALL DIAGZ(ISIZEL,ISIZER,W(LXA),IDUM,W(LXA),0)
C**INTEGRATE OVER SINGLE NORMAL COORDINATE
      K2=0
      K3=0
      DO K=1,NMODE
        IF(K.EQ.1.AND.ITIM.EQ.0)ITIM1A=0
        CALL INTARR(W(LNBF),W(LMBF),K,IK1,IK2,IK3)
        IF(ICOUPL.EQ.0)GO TO 801
C**CORIOLIS
        IF(IABC.LT.7)THEN
          IHERM=1
C**READ INTO MATRIX
          CALL MEMO(1,LXA1,KXA1,0,0,0,0,0,0,0,0)
          CALL MATIN(W(LXA1),W(LXA1),ISIZE1,21)
          CALL VMI1(NMODE,K,W(LH+K2),W(LXQ+K3),W(LXA),W(LXA1),IK3,
     1    IK2,W(LIPL),W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP1),ISIZE1)
          IF(ICOUPL.GE.0)CALL MEMO(-1,LXA1,KXA1,0,0,0,0,0,0,0,0)
        ELSE
          IHERM=-1
        END IF
        IF(ICOUPL.EQ.1)GO TO 801
C**INTEGRATE OVER TWO NORMAL COORDINATES
        L2=0
        L3=0
        DO L=1,K-1
          IF(L.EQ.1.AND.K.EQ.2.AND.ITIM.EQ.0)ITIM2A=0
          CALL INTARR(W(LNBF),W(LMBF),L,IL1,IL2,IL3)
C**CORIOLIS
C**READ INTO MATRIX
          JCI=NBAS(2)+1
          CALL MEMO(1,LXA2,KXA2,0,0,0,0,0,0,0,0)
          CALL MATIN(W(LXA2),W(LXA2),ISIZE2,22)
          CALL VMI2(NMODE,K,L,W(LH+K2),W(LXQ+K3),W(LH+L2),
     1    W(LXQ+L3),IK3,IK2,IL3,IL2,W(LXA),W(LXA2),
     2    W(LIPL),W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP2),ISIZE2,JCI)
          IF(ICOUPL.GT.1)CALL MEMO(-1,LXA2,KXA2,0,0,0,0,0,0,0,0)
          IF(ICOUPL.EQ.2)GO TO 802
C**INTEGRATE OVER THREE NORMAL COORDINATES
          N2=0
          N3=0
          DO N=1,L-1
            IF(N.EQ.1.AND.L.EQ.2.AND.K.EQ.3.AND.ITIM.EQ.0)ITIM3A=0
            CALL INTARR(W(LNBF),W(LMBF),N,IN1,IN2,IN3)
C**CORIOLIS
C**READ INTO MATRIX
            JCI=NBAS(3)+1
            CALL MEMO(1,LXA3,KXA3,0,0,0,0,0,0,0,0)
            CALL MATIN(W(LXA3),W(LXA3),ISIZE3,23)
            CALL VMI3(NMODE,K,L,N,W(LH+K2),W(LXQ+K3),
     1      W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),IK3,IK2,IL3,IL2,
     2      IN3,IN2,W(LXA),W(LXA3),
     3      W(LIPL),W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP3),ISIZE3,JCI)
            IF(ICOUPL.GT.2)CALL MEMO(-1,LXA3,KXA3,0,0,0,0,0,0,0,0)
            IF(ICOUPL.EQ.3)GO TO 803
C**INTEGRATE OVER FOUR NORMAL COORDINATES
            M2=0
            M3=0
            DO M=1,N-1
              IF(M.EQ.1.AND.N.EQ.2.AND.L.EQ.3.AND.K.EQ.4.AND.
     1        ITIM.EQ.0)ITIM4A=0
              CALL INTARR(W(LNBF),W(LMBF),M,IM1,IM2,IM3)
C**CORIOLIS
C**READ INTO MATRIX
              JCI=NBAS(4)+1
              CALL MEMO(1,LXA4,KXA4,0,0,0,0,0,0,0,0)
              CALL MATIN(W(LXA4),W(LXA4),ISIZE4,24)
              CALL VMI4(NMODE,K,L,N,M,W(LH+K2),W(LXQ+K3),
     1        W(LH+L2),W(LXQ+L3),W(LH+N2),W(LXQ+N3),W(LH+M2),
     2        W(LXQ+M3),IK3,IK2,IL3,IL2,IN3,IN2,IM3,IM2,W(LXA),
     3        W(LXA4),W(LIPL),W(LIPR),ISIZMX,ISIZEL,ISIZER,W(LIP4),
     4        ISIZE4,JCI)
              IF(ICOUPL.GT.3)CALL MEMO(-1,LXA4,KXA4,0,0,0,0,0,0,0,0)
              IF(ICOUPL.EQ.4)GO TO 804
C**5-MODE COUPLING HERE IF NEEDED
804   CONTINUE
              M2=M2+IM1
              M3=M3+IM2
            END DO
803   CONTINUE
            N2=N2+IN1
            N3=N3+IN2
          END DO
802   CONTINUE
          L2=L2+IL1
          L3=L3+IL2
        END DO
801   CONTINUE
        K2=K2+IK1
        K3=K3+IK2
      END DO
      CALL DUMVM(W(LXA),ISIZEL,ISIZER,30+IABC)
      CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
9999  CONTINUE
5503  CONTINUE
5502  CONTINUE
5501  CONTINUE
5500  CONTINUE
C**********************************************************************
C**                                  SET UP OFF-DIAGONAL MATRIX FOR K>0
C**********************************************************************
      IF(NVSYM.EQ.1)THEN
C**H2CO C1
        NRSYM=1
C**H2O
        IF(TRIAT)NRSYM=2
      END IF
      IF(NVSYM.EQ.2)THEN
C**H2CO Cs
        NRSYM=2
C**H2O
        IF(TRIAT)NRSYM=4
      END IF
      IF(NVSYM.EQ.4)THEN
C**H2CO C2v
        NRSYM=4
      END IF
C**NRSYM IS TOTAL NUMBER ROTATIONAL SYMMETRIES
C**INDEX IS A COUNT OF THE ACTUAL ROTATIONAL SYMMETRIES REQUIRED
      INDEX=1
      DO 2222 NROT=1,NRSYM
C**********************************************************************
C**                                                LOOP OVER SYMMETRIES
C**********************************************************************
C**SKIP IF FINISHED
      IF(INDEX.GT.MVSYM)GO TO 2222
C**SKIP UNTIL FIND NEXT SYMMETRY THAT IS REQUIRED
      IF(NROT.NE.MWSYM(INDEX))GO TO 2222
      INDEX=INDEX+1
      WRITE(IOUT,375)NROT
      IF(TRIAT)THEN
C**START B2 (IF EXISTS)
        NS=2
        IF(NROT.LE.NVSYM)THEN
          ISIZE=(JTHIS+1)*NVAL
          K21=JTHIS+1
C**ROTATION ELEMENTS 1,3,5,...
          KOFF=1
C**START A1 (A')
          IF(NROT.EQ.1)NS=1
        ELSE
          ISIZE=JTHIS*NVAL
          K21=JTHIS
C**ROTATION ELEMENTS 2,4,6,...
          KOFF=2
C**START A1 (A')
          IF(NROT.EQ.NVSYM+1)NS=1
        END IF
      ELSE
        NS=NROT
        K21=2*JTHIS+1
        ISIZE=K21*NVAL
      END IF
      WRITE(IOUT,380)ISIZE,CUT
C**ISTART, IEND ARE START AND END COLUMNS
      ISTART=1
      KSTART=1
      IF(LANCZ)THEN
        CALL MEMO(3,LWK,ISIZE,LYK,ISIZE,LZK,ISIZE,0,0,0,0)
        KLAN=LANMAX*(LANMAX+1)/2
        LSIZE=ISIZE
        KSIZE=ISIZE
      ELSE
        KXA=ISIZE*(ISIZE+1)/2
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGZ(ISIZE,ISIZE,W(LXA),ISIZE,W(LXA),ICI)
        IEND=ISIZE
        KEND=K21
      END IF
      ISKIP=0
7777  CONTINUE
      IF(LANCZ)THEN
        IEND=ISTART-1
        KEND=KSTART-1
        KXA=0
        KCOL=0
        DO I=ISTART,ISIZE
          KXA=KXA+LSIZE
C**CAN'T GET IT ALL IN
          IF(KXA.GT.KLAN)THEN
C**SUBTRACT TOTAL SIZE 'K' BLOCK SO FAR
            KXA=KXA-KSIZE
C**SET LSIZE TO VALUE FOR START OF 'K' BLOCK
            LSIZE=LSIZE+KCOL
C**STARTING VALUE OF KSIZE FOR NEXT 'K' BLOCK (SIZE FIRST COLUMN)
            KSIZE=LSIZE
C**RESET IEND TO LAST COLUMN PREVIOUS BLOCK
            IEND=IEND-KCOL
            GO TO 7776
          END IF
C**UPDATE FOR NEXT COLUMN
          LSIZE=LSIZE-1
C**KSIZE: TOTAL SIZE 'K' BLOCK INCLUDING NEXT COLUMN
          KSIZE=KSIZE+LSIZE
C**IEND: CURRENT END COLUMN
          IEND=IEND+1
C**KCOL: NUMBER COLUMNS IN 'K' BLOCK SO FAR
          KCOL=KCOL+1
          IF(MOD(KCOL,NVAL).EQ.0)THEN
C**STARTING VALUES FOR NEXT BLOCK
            KCOL=0
            KSIZE=LSIZE
C**KEND: CURRENT KROT(LHS) FOR COMPLETED BLOCK
            KEND=KEND+1
          END IF
        END DO
7776    CONTINUE
        IF(IEND.EQ.0)STOP 'LANMAX TOO SMALL'
        WRITE(IOUT,*)'ISTART,IEND,KSTART,KEND,KXA',
     1  ISTART,IEND,KSTART,KEND,KXA
        CALL FLUSH(IOUT)
        CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
C**ZEROISE MATRIX
        CALL DIAGL(KXA,W(LXA))
      END IF
      CALL MEMO(2,LTEMP,ISIZMX*NVAL,LXK,ISIZMX*ISIZMX,0,0,0,0,0,0)
      DO 8888 KROT=KSTART,KEND
      WRITE(IOUT,390)KROT,NS
C**********************************************************************
C**                                                       LOOP ROUND Ka
C**********************************************************************
      IF(TRIAT)THEN
        JROT=2*(KROT-1)+KOFF
      ELSE
        JROT=KROT
      END IF
C**SET UP DIAGONAL ELEMENTS
      CALL DIAGEL(ISIZE,KROT,W(LXA),W(LEVCI),NVAL,KEL21,JROT,NS,ISTART,
     1KSTART)
C**ISIZEL IS SIZE LHS CI BASIS
C**ISIZER IS SIZE RHS CI BASIS
C**NSL IS SYMMETRY LHS
C**NSR IS SYMMETRY RHS
C**KROTL DENOTES LHS BLOCK OFFSET
C**KROTR DENOTES RHS BLOCK OFFSET
C**JROTL IS LHS ROTATIONAL FUNCTION
C**JROTR IS RHS ROTATIONAL FUNCTION
      NSL=NS
      ISIZEL=NTOT(NSL)
      KROTL=KROT
      IF(TRIAT)THEN
        JROTL=2*(KROTL-1)+KOFF
      ELSE
        JROTL=KROTL
      END IF
      JROTR=JROTL
      KSIGN=1
      IF(KSTEP.NE.1)THEN
        NSR=NSL
        ISIZER=ISIZEL
        KROTR=KROTL
        IELL=1
        IELX=(JROTL-1)/KSTEP
        IELL=IELL+IELX*KSTEP
        IELR=IELL
        KSIGN=-1
C**NEGATE OLD DIAGONAL ELEMENTS FOR Jx**2 + Jy**2 + Jz**2 (1, 2 AND 3)
C**Jx**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,IELL,IELR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,1)
C**Jy**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,IELL,IELR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,2)
C**Jz**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,IELL,IELR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,3)
C**SET UP NEW DIAGONAL ELEMENTS FOR Jx**2 + Jy**2 + Jz**2 (1, 2 AND 3)
        KSIGN=1
C**Jx**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,1)
C**Jy**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,2)
C**Jz**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,3)
      END IF
      IF(TRIAT)THEN
        NSR=NSL+1
        IF(NSR.GT.NVSYM)NSR=1
        KINCR=1
        KROTR=KROTL+KINCR
        JROTR=2*(KROTR-1)+KOFF
      ELSE
        NSR=NSL
        IF(NVSYM.EQ.4)THEN
          NSR=NSL+(-1)**(NROT+1)
          IF(MOD(KROT,4).EQ.0)NSR=NSL-(-1)**(NROT+1)
        END IF
        KINCR=3
        KROTR=KROTL+2
        JROTR=JROTL+2
      END IF
      ISIZER=NTOT(NSR)
      IF(MOD(KROT,2).EQ.0.AND.KROT+KINCR.LE.K21)THEN
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jx (7)
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,7)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JyJz]+ (5)
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,5)
        IF(.NOT.TRIAT)THEN
          NSR=NSL+1
          IF(NSR.GT.NVSYM)NSR=1
          IF(NVSYM.EQ.4)THEN
            NSR=NSL+2
            IF(NSR.GT.NVSYM)THEN
              INCR=NSR-NVSYM
              NSR=INCR
            END IF
          END IF
          KROTR=KROTL+KINCR
          JROTR=JROTL+KINCR
          ISIZER=NTOT(NSR)
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jy (8)
          CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1    (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2    ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,8)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJz]+ (4)
          CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1    (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2    ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,4)
        END IF
      END IF
      IF(.NOT.TRIAT)THEN
        NSR=NSL
        IF(NVSYM.EQ.4)THEN
          NSR=NSL-(-1)**(NROT+1)
          IF(MOD(KROT-1,4).EQ.0)NSR=NSL+(-1)**(NROT+1)
        END IF
        KINCR=2
        KROTR=KROTL+KINCR
        JROTR=JROTL+KINCR
        ISIZER=NTOT(NSR)
      END IF
      IF(MOD(KROT,2).EQ.1.AND.KROT+KINCR.LE.K21)THEN
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jx (7)
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,7)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JyJz]+ (5)
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,5)
        IF(.NOT.TRIAT)THEN
          NSR=NSL+1
          IF(NSR.GT.NVSYM)NSR=1
          IF(NVSYM.EQ.4)THEN
            NSR=NSL+2
            IF(NSR.GT.NVSYM)THEN
              INCR=NSR-NVSYM
              NSR=INCR
            END IF
          END IF
          KROTR=KROTL+1
          JROTR=JROTL+1
          ISIZER=NTOT(NSR)
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jy (8)
          CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1    (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2    ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,8)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJz]+ (4)
          CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1    (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2    ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,4)
        END IF
      END IF
      NSR=NSL
      IF(TRIAT)THEN
        KINCR=2
        KROTR=KROTL+KINCR
        JROTR=2*(KROTR-1)+KOFF
      ELSE
        KINCR=4
        KROTR=KROTL+KINCR
        JROTR=JROTL+KINCR
      END IF
      ISIZER=NTOT(NSR)
      IF(KROT+KINCR.LE.K21)THEN
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jx**2 + Jy**2 (1 AND 2)
C**Jx**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,1)
C**Jy**2
        CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1  (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2  ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,JROTR,
     3  JROTL,JROTR,NSL,NSR,ISTART,KSTART,2)
      END IF
      IF(.NOT.TRIAT)THEN
        NSR=NSL+1
        IF(NSR.GT.NVSYM)NSR=1
        IF(NVSYM.EQ.4)THEN
          NSR=NSL+(-1)**NROT
          IF(MOD(KROT,4).EQ.0)NSR=NSL-(-1)**NROT
          IF(NSR.GT.NVSYM)NSR=1
          IF(NSR.EQ.0)NSR=NVSYM
        END IF
        KINCR=1
        KROTR=KROTL+KINCR
        JROTR=JROTL+KINCR
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).EQ.0)THEN
C**SET UP OFF-DIAGONAL ELEMENTS FOR Jz (9)
          CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1    (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2    ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,9)
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJy]+ (6)
          CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1    (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2    ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,6)
          KINCR=5
          KROTR=KROTL+KINCR
          JROTR=JROTL+KINCR
          IF(KROT+KINCR.LE.K21)THEN
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJy]+ (6)
            CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1      (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2      ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3      JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,6)
          END IF
        END IF
        IF(NVSYM.EQ.4)THEN
          NSR=NSL+(-1)**(NROT+1)
          IF(MOD(KROT-1,4).EQ.0)NSR=NSL-(-1)**(NROT+1)
        END IF
        KINCR=4
        KROTR=KROTL+3
        JROTR=JROTL+3
        ISIZER=NTOT(NSR)
        IF(MOD(KROT,2).NE.0.AND.KROT+KINCR.LE.K21)THEN
C**SET UP OFF-DIAGONAL ELEMENTS FOR [JxJy]+ (6)
          CALL OFFDEL(ISIZE,W(LXA),W(LCFS+ISIZMX*NVAL*KEL21*
     1    (NSL-1)),W(LCFS+ISIZMX*NVAL*KEL21*(NSR-1)),ISIZMX,W(LTEMP),
     2    ISIZEL,ISIZER,NVAL,W(LXK),W(LSS),J21,KROTL,KROTR,JROTL,
     3    JROTR,JROTL,JROTR,NSL,NSR,ISTART,KSTART,6)
        END IF
      END IF
C**UPDATE VIBRATIONAL SYMMETRY FOR NEXT TIME IF REQUIRED
      IF((.NOT.TRIAT).AND.NVSYM.EQ.4)THEN
        IF(MOD(KROT,2).NE.0)THEN
          NS=NS+2
        ELSE
          IF(MOD(KROT,4).NE.0)THEN
            NS=NS+(-1)**NROT
          ELSE
            NS=NS-(-1)**NROT
          END IF
        END IF
        IF(NS.GT.NVSYM)THEN
          INCR=NS-NVSYM
          NS=INCR
        END IF
        IF(NS.EQ.0)NS=NVSYM
      ELSE
        NS=NS+1
        IF(NS.GT.NVSYM)NS=1
      END IF
8888  CONTINUE
      ISKIP=0
      IF(LANCZ)THEN
        CALL LANOUT(W(LXA),W(LWK),W(LWK),W(LYK),W(LZK),W(LZK),ISIZE,
     1  ISTART,IEND)
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        IF(IEND.NE.ISIZE)THEN
          ISTART=IEND+1
          KSTART=KEND+1
          ISKIP=1
        END IF
      END IF
      CALL MEMO(-2,LTEMP,ISIZMX*NVAL,LXK,ISIZMX*ISIZMX,0,0,0,0,0,0)
      IF(ISKIP.NE.0)GO TO 7777
      IF(LANCZ)CALL MEMO(-1,LZK,ISIZE,0,0,0,0,0,0,0,0)
C**********************************************************************
C**                                                     MATRIX COMPLETE
C**********************************************************************
      WRITE(IOUT,385)
      IF(NVALR.GT.ISIZE)NVALR=ISIZE
      NVECR=NVALR
      KXK=ISIZE*NVALR
      KSUP4=5*ISIZE
      KWK=ISIZE
      KEVAL=NVALR
      CALL TIMIT(1)
      ICID=1
      IF(LANCZ)THEN
        LGIV=(IGIV.NE.0)
        KXKL=NVALR*NCYCLE*NVALR*NCYCLE
        IF(LGIV)KXKL=NVALR*NCYCLE*NVALR
        IF(LANZA)THEN
          CALL MEMO(3,LXK,ISIZE,LVK,ISIZE,LZK,ISIZE,0,0,0,0)
          CALL MEMO(5,LSK,ISIZE,LEN,NVALR,LEL,NVALR,LNV,NVALR,LWKL,
     1    ISIZE)
          CALL MEMO(4,LWRK,NVALR*NCYCLE,
     1    LSUP4,5*NVALR*NCYCLE,LEVAL,NVALR*NCYCLE,LXKL,KXKL,0,0)
          KXA=NVALR*NCYCLE*(NVALR*NCYCLE+1)/2
          IF(LGIV)CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
          CALL MEMO(1,LXKLC,KXA,0,0,0,0,0,0,0,0)
          WRITE(IOUT,*)'Calculating LANCZOS'
          IF(LGIV)THEN
            CALL LANCZA(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVALR,W(LJP),JSIZE,NMODE,W(LASSIG),ISIZMX,KEL21,NROT,
     2      W(LEN),W(LEL),W(LNV),W(LXA),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          ELSE
            CALL LANCZA(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVALR,W(LJP),JSIZE,NMODE,W(LASSIG),ISIZMX,KEL21,NROT,
     2      W(LEN),W(LEL),W(LNV),W(LXKL),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          END IF
          CALL MEMO(-5,LXK,ISIZE,LVK,ISIZE,LZK,ISIZE,LWK,ISIZE,LYK,
     1    ISIZE)
          CALL MEMO(-5,LSK,ISIZE,LEN,NVALR,LEL,NVALR,LNV,NVALR,LXKL,
     1    KXKL)
          CALL MEMO(-5,LXKLC,KXA,LWRK,NVALR*NCYCLE,
     1    LSUP4,5*NVALR*NCYCLE,LEVAL,NVALR*NCYCLE,LWKL,ISIZE)
        ELSE
          CALL MEMO(3,LXK,ISIZE*NVALR,LVK,ISIZE*NVALR,LZK,ISIZE*NVALR,
     1    0,0,0,0)
          CALL MEMO(5,LSK,ISIZE*NVALR,LEN,NVALR,LEL,NVALR,LNV,NVALR,
     1    LWKL,ISIZE*NVALR)
          CALL MEMO(4,LWRK,NVALR*NCYCLE,
     1    LSUP4,5*NVALR*NCYCLE,LEVAL,NVALR*NCYCLE,LXKL,KXKL,0,0)
          KXA=NVALR*NCYCLE*(NVALR*NCYCLE+1)/2
          IF(LGIV)CALL MEMO(1,LXA,KXA,0,0,0,0,0,0,0,0)
          CALL MEMO(1,LXKLC,KXA,0,0,0,0,0,0,0,0)
          WRITE(IOUT,*)'Calculating LANCZOS'
          IF(LGIV)THEN
            CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVALR,W(LJP),JSIZE,NMODE,W(LASSIG),ISIZMX,KEL21,NROT,
     2      W(LEN),W(LEL),W(LNV),W(LXA),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          ELSE
            CALL LANCZB(W(LXK),W(LVK),W(LZK),W(LWK),W(LYK),W(LSK),
     1      ISIZE,NVALR,W(LJP),JSIZE,NMODE,W(LASSIG),ISIZMX,KEL21,NROT,
     2      W(LEN),W(LEL),W(LNV),W(LXKL),W(LXKLC),W(LXKL),W(LWRK),
     3      W(LSUP4),W(LEVAL),NCYCLE,W(LWKL))
          END IF
          CALL MEMO(-5,LXK,ISIZE*NVALR,LVK,ISIZE*NVALR,LZK,ISIZE*NVALR,
     1    LWK,ISIZE,LYK,ISIZE)
          CALL MEMO(-5,LSK,ISIZE*NVALR,LEN,NVALR,LEL,NVALR,LNV,NVALR,
     1    LXKL,KXKL)
          CALL MEMO(-5,LXKLC,KXA,LWRK,NVALR*NCYCLE,
     1    LSUP4,5*NVALR*NCYCLE,LEVAL,NVALR*NCYCLE,LWKL,ISIZE*NVALR)
        END IF
        IF(LGIV)CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
        LGIV=.FALSE.
      ELSE
        LGIV=.TRUE.
        CALL MEMO(4,LWK,KWK,LEVAL,KEVAL,LSUP4,KSUP4,LXK,KXK,0,0)
        WRITE(IOUT,*)'Calculating DIAG'
        CALL DIAG(W(LXA),W(LXK),ISIZE,ISIZE,0,W(LSUP4),W(LEVAL),W(LWK),
     1  NVALR,NVECR,W(LJP),JSIZE,NMODE,W(LXK),W(LXK),W(LASSIG),ISIZMX,
     2  KEL21,NROT)
        CALL MEMO(-4,LXK,KXK,LWK,KWK,LEVAL,KEVAL,LSUP4,KSUP4,0,0)
        LGIV=.FALSE.
        CALL MEMO(-1,LXA,KXA,0,0,0,0,0,0,0,0)
      END IF
      ICID=0
      CALL TIMIT(3)
2222  CONTINUE
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETQ(XK,EVAL,SOL,WRK,N,XQ,XV,M,QM,NMODE,MODE,INDEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XK(N,N),EVAL(N,1),SOL(N,1),WRK(4*N),XQ(M),XV(M)
      DIMENSION QM(NMODE)
      COMMON/XACC/XACC
      COMMON/NFIT/NFIT
      COMMON/FILASS/IOUT,INP
      EXTERNAL LSFUN1,LPRINT
100   FORMAT(/,1X,'DGESV IFAIL = ',I3)
101   FORMAT(4D20.10)
102   FORMAT(5X,'SHOULD NOT OCCUR')
103   FORMAT(1X,'SUMSQ = ',D20.10,/)
104   FORMAT(/,1X,'E04FBF IFAIL = ',I3)
105   FORMAT(1X,'COEFFICIENTS OF POLYNOMIAL: 1, Q, Q**2,...:')
106   FORMAT(/,1X,'MODE: ',I3,/)
107   FORMAT(/,1X,'VALUE OF Q AT MINIMUM: ',D20.10,/)
108   FORMAT(/,1X,'MODE: ',I3,/)
      IF(M.LT.3)STOP 'TOO FEW POINTS'
      IF(M.LT.N)RETURN
      WRITE(IOUT,108)MODE
      INDEX=0
      ISTART=(M-N)/2
      DO I=1,N
        EVAL(I,1)=XV(ISTART+I)
        DO J=1,N
          XK(I,J)=XQ(ISTART+I)**(J-1)
        END DO
      END DO
      IFAIL=1
CCCC  CALL F04AAF(XK,N,EVAL,N,N,1,SOL,N,WRK,IFAIL)
C**NEW
      CALL DGESV(N,1,XK,N,SOL,EVAL,N,IFAIL)
C**NEW
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,100)IFAIL
CCCC    STOP 'ERROR IN F04AAF'
C**NEW
        STOP 'ERROR IN DGESV'
C**NEW
      END IF
C**NEW
      DO I=1,N
        SOL(I,1)=EVAL(I,1)
      END DO
C**NEW
      WRITE(IOUT,105)
      WRITE(IOUT,101)(SOL(I,1),I=1,N)
      CALL FLUSH(IOUT)
C**************************************************************
      IFAIL=1
      MFIT=1
      NFIT=N
      XACC=1.D-15
      FTOL=XACC
      XTOL=1.D-10
C**TEMPORARY
C     STEP=10.D0
      STEP=1.D0
C**TEMPORARY
      MAXCAL=1000
      KPRINT=1000
      IWXY=2*MFIT*(MFIT+MFIT)+2*MFIT+5*MFIT
      IF(IWXY.GT.4*N)THEN
        WRITE(IOUT,102)
        STOP
      END IF
C**ANSWER NEAR ZERO, SO INITIAL GUESS SET TO ZERO
      QM(MODE)=0
      CALL E04FBF(MFIT,MFIT,QM(MODE),SOL,SUMSQ,FTOL,XTOL,STEP,WRK,4*N,
     1LSFUN1,LPRINT,KPRINT,MAXCAL,IFAIL)
      IF(IFAIL.NE.0)THEN
        WRITE(IOUT,104)IFAIL
        STOP 'ERROR IN E04FBF'
      END IF
      WRITE(IOUT,107)QM(MODE)
C     WRITE(IOUT,103)SUMSQ
C**************************************************************
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE LSFUN1(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(M)
      COMMON/NFIT/NFIT
      F(1)=F(2)
      DO I=3,NFIT
        K=I-1
        J=K-1
        F(1)=F(1)+F(I)*K*X(1)**J
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE LPRINT(M,N,X,F,S,NCALL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**MONITOR PROGRESS OF E04FBF IF REQUIRED
      DIMENSION X(N),F(M)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETMIN(QM,NMODE,NATOM,XX,X0,XM,XL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION QM(NMODE),XX(NATOM,3),X0(NATOM,3),XM(NATOM)
      DIMENSION XL(NATOM,NMODE,3)
      COMMON/FILASS/IOUT,INP
100   FORMAT(50(1H#))
      WRITE(IOUT,100)
C**CARTESIAN COORDINATES
      DO I=1,NATOM
        DO J=1,3
          XX(I,J)=X0(I,J)
          DO K=1,NMODE
            XX(I,J)=XX(I,J)+XL(I,K,J)*QM(K)/SQRT(XM(I))
          END DO
        END DO
      END DO
      DO I=1,NATOM
        DO J=1,3
          X0(I,J)=XX(I,J)
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETV(M,XQ,XV,MODE,NMODE,NATOM,QQ,
     1RR,XX,X0,XL,XM,NPOT,IPOT,JPOT,CPOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XQ(M),XV(M)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/WHICH/IWHICH
      DO J=1,M
C**ENERGY THIS POINT
        DO K=1,NMODE
          QQ(K)=0
        END DO
        QQ(MODE)=XQ(J)
        IF(IWHICH.NE.0)THEN
          DO I=1,NATOM
            DO K=1,3
              XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1        SQRT(XM(I))
            END DO
          END DO
          CALL GETPOT(VPT,NATOM,XX,RR)
        ELSE
          CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,NMODE,QQ)
        END IF
        XV(J)=VPT
C**ENERGY THIS POINT
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE GETNOR(NATOM,NMODE,XM,X0,OMEGA,XL,XX,RR,
     1IPOT,JPOT,CPOT,NPOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      DIMENSION XM(NATOM),X0(NATOM,3),XL(NATOM,NMODE,3),XX(NATOM,3)
      DIMENSION RR(NATOM,NATOM),OMEGA(NMODE)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/PRINT/IPRINT,JPRINT
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F20.12,/)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
250   FORMAT(/,1X,'MODE(k) = ',I4)
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
      CALL NORMAL(NATOM,NMODE,XM,X0,OMEGA,XL,XX,RR)
C*************************************************************
C*************************************************************
C**CALL USER-SUPPLIED ROUTINE TO GET 'CLEAN' VECTORS IF E OR F
      CALL ROTATE(NATOM,X0,OMEGA,NMODE,XL,WAVENM)
C*************************************************************
C*************************************************************
      IF(IPRINT.GT.0)WRITE(IOUT,245)
      DO J=1,NMODE
        IF(IPRINT.GT.0)WRITE(IOUT,250)J
        DO I=1,NATOM
          IF(IPRINT.GT.0)WRITE(IOUT,255)I,(XL(I,J,K),K=1,3)
        END DO
      END DO
      DO K=1,NMODE
        DO I=1,6
          IPOT(K,I)=0
          JPOT(K,I)=0
        END DO
        IPOT(K,1)=K
        JPOT(K,1)=2
        CPOT(K)=OMEGA(K)*OMEGA(K)/2
        IF(IPRINT.GT.0)WRITE(IOUT,220)IPOT(K,1),OMEGA(K)*WAVENM
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE INPUT(NPOT,IPOT,JPOT,CPOT,NMODE,OMEGA,XX,NBF,MBF,NSTAT
     1,ISTAT,NATOM,XM,X0,XL,RR,ICI,MAXBAS,MAXJ,JSTAT,KSTAT,ESTAT,WSTAT,
     2NVF,INORM,MCHECK,INDEX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION JSTAT(NSTAT,NMODE),KSTAT(NMODE),ESTAT(NSTAT)
      DIMENSION WSTAT(NSTAT),RR(NATOM,NATOM)
      DIMENSION NVF(NMODE),MAXBAS(NMODE,MAXJ)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      DIMENSION NBF(NMODE),MBF(NMODE),OMEGA(NMODE),ISTAT(NSTAT,NMODE)
      DIMENSION XM(NATOM),X0(NATOM,3),XL(NATOM,NMODE,3),XX(NATOM,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/WHICH/IWHICH
      COMMON/PATH/ISCFCI
      COMMON/ESTATE/IORDER
      COMMON/FILASS/IOUT,INP
200   FORMAT(/,1X,'FOR MODE ',I4,/,
     1         1X,'NUMBER PRIMITIVE FUNCTIONS = ',I4,/,
     2         1X,'NUMBER INTEGRATION POINTS = ',I4,/,
C*****************************************************HEG
C*****************************************************HEG
     3         1X,'NUMBER CONTRACTED FUNCTIONS = ',I4,/)
C*****************************************************HEG
C*****************************************************HEG
205   FORMAT(1X,'NO. OF PRIMITIVES CHANGED FROM ',I4,' TO ',I4)
210   FORMAT(/,1X,'POTENTIAL TERMS: IPOT, JPOT, CPOT',/)
215   FORMAT(1X,6I2,1X,6I2,1X,D20.12)
220   FORMAT(/,1X,'MODE',I4,'  OMEGA = ',F7.2,/)
225   FORMAT(//,1X,'NORMAL COORDINATE ANALYSIS',/)
230   FORMAT(1X,'ATOM ',I2,':  MASS =',F12.7)
235   FORMAT(//,1X,'EQUILIBRIUM GEOMETRY',/)
240   FORMAT(1X,'ATOM ',I2,':  X0 =',F12.7,'  Y0 =',F12.7,'  Z0 =',
     1F12.7)
245   FORMAT(//,1X,'DISPLACEMENTS OF NORMAL MODES, L(alpha i,k)')
250   FORMAT(/,1X,'MODE(k) = ',I4)
255   FORMAT(1X,'ATOM(i) ',I2,':  x =',F12.6,'  y =',F12.6,'  z =',
     1F12.6)
260   FORMAT(1X,'NO. OF CONTRACTED FUNCTIONS CHANGED FROM ',I4,' TO ',
     1I4)
      IF(INDEX.NE.0)GO TO 9999
C**NO. PRIMITIVE QUANTA AND ASSOCIATED INTEGRATION POINTS
      READ(INP,*)
C****************************************************************************
C**NBF(NMODE),MBF(NMODE),NVF(NMODE)
C**            Input NMODE records of three integers
C**NBF(K):     Number of quanta in harmonic-oscillator basis for mode K
C**MBF(K):     Number of Gauss Hermite integration points for mode K
C**NVF(K):     Number of quanta of contracted numerical functions for mode K
C****************************************************************************
      DO K=1,NMODE
C*****************************************************HEG
C*****************************************************HEG
        READ(INP,*)NBF(K),MBF(K),NVF(K)
C**NUMBER OF FUNCTIONS ONE GREATER THAN NUMBER OF QUANTA
        NBF(K)=NBF(K)+1
        NVF(K)=NVF(K)+1
        IF(ICI.LT.0.AND.ISCFCI.GT.0)THEN
          IF(MAXJ.EQ.0)THEN
            IF(-ICI.GT.NBF(K))ICI=-NBF(K)
            NVF(K)=-ICI
          ELSE
            CALL HEG2(MAXBAS,NMODE,MAXJ,NVF,K)
          END IF
        END IF
        IF(NVF(K).LE.0)NVF(K)=NBF(K)
        IF(NVF(K).GT.NBF(K))THEN
          WRITE(IOUT,260)NBF(K),NVF(K)
          NBF(K)=NVF(K)
        END IF
C**KINETIC ENERGY INTEGRALS (AT LEAST) EXACT IF MBF.GE.NBF
        IF(MBF(K).LT.NBF(K))THEN
          WRITE(IOUT,205)MBF(K),NBF(K)
          MBF(K)=NBF(K)
        END IF
        IF(NVF(K).LT.NBF(K))THEN
          INCR=MBF(K)-NBF(K)
          IF(NBF(K).LT.NVF(K)+INCR)NBF(K)=NVF(K)+INCR
          IF(MBF(K).LT.NBF(K)+INCR)MBF(K)=NBF(K)+INCR
        END IF
C**MAKE NO, POINTS ODD
        IF(MOD(MBF(K),2).EQ.0)MBF(K)=MBF(K)+1
        WRITE(IOUT,200)K,NBF(K),MBF(K),NVF(K)
C*****************************************************HEG
C*****************************************************HEG
      END DO
      READ(INP,*)
C****************************************************************************
C**IPOT(NPOT),JPOT(NPOT),CPOT(NPOT)
C**         Normal coordinate force field, defined by:
C**         Jelski et al in J. Comput. Chem., 17, 1645 (1996).
C**         A maximum of 6 coupled modes is allowed.
C**IPOT:    (first 6 integers): defines the Normal mode(s) involved.
C**JPOT:    (second 6 integers): defines the number of quanta in the
C**         corresponding mode(s)
C**CPOT:    defines the potential contribution F(ijkl..) in hartrees
C**
C**These parameters are input under the following conditions:
C**IWHICH=0: Input NPOT terms IPOT,JPOT,CPOT where the first NMODE terms are the
C**          NMODE harmonic force constants F(ii)=0.5*omega(i)**2
C**IWHICH=1: Input NMODE terms corresponding to the NMODE F(ii) if INORM=0
C**          Input NOTHING if INORM=1
C**          note F(ii) is omega(i) in cm-1
C****************************************************************************
      IF(IWHICH.EQ.0.OR.INORM.EQ.0)THEN
        WRITE(IOUT,210)
        DO I=1,NPOT
C**MAXIMUM OF 6 MODES COUPLED
          READ(INP,*)(IPOT(I,J),J=1,6),(JPOT(I,J),J=1,6),CPOT(I)
          WRITE(IOUT,215)(IPOT(I,J),J=1,6),(JPOT(I,J),J=1,6),
     1    CPOT(I)
          IF(I.LE.NMODE)THEN
            IF(IWHICH.EQ.0)OMEGA(I)=SQRT(2*CPOT(I))
            IF(IWHICH.NE.0)OMEGA(I)=CPOT(I)/WAVENM
            WRITE(IOUT,220)IPOT(I,1),OMEGA(I)*WAVENM
          END IF
        END DO
      END IF
C**********************************************************
C**********************************************************
      READ(INP,*)
      IF(IORDER.EQ.0)THEN
C****************************************************************************
C**NDEF
C**         SCF state definition
C**
        READ(INP,*)NDEF
C**
C**         NSTAT +ve and NDEF -ve:
C**         Input NSTAT records ISTAT(NSTAT,NMODE) corresponding to NSTAT
C**         specific SCF state definitions (see Jelski et al)
C**
C**         NSTAT +ve and NDEF=0 or NDEF +ve:
C**         Input NSTAT records N1 N2 where N1 are the number of quanta for
C**         mode N2
C**         All remaining modes will take on NDEF quanta (see Jelski et al)
C**
C**         NSTAT -ve:
C**         Input NOTHING. In this case NDEF is a dummy parameter.
C**         It is recommended that NSTAT is set -ve for all VCI calculations
C**         (see ICI -ve above)
C****************************************************************************
        READ(INP,*)
C**ORIGINAL JELSKI
        IF(NDEF.LT.0)THEN
          DO I=1,NSTAT
            READ(INP,*)(ISTAT(I,K),K=1,NMODE)
          END DO
        ELSE
          DO K=1,NMODE
            ISTAT(1,K)=NDEF
          END DO
          DO I=2,NSTAT
            READ(INP,*)N1,N2
            DO K=1,NMODE
              IF(K.NE.N2)THEN
                ISTAT(I,K)=NDEF
              ELSE
                ISTAT(I,K)=N1
              END IF
            END DO
          END DO
        END IF
      WRITE(IOUT,*)
      WRITE(IOUT,*)'COMPUTED STATES'
      DO I=1,NSTAT
        WRITE(IOUT,*)(ISTAT(I,K),K=1,NMODE)
      END DO
      END IF
C**********************************************************
C**********************************************************
C**READ NORMAL COORDINATE DETAILS
      WRITE(IOUT,225)
      READ(INP,*)
C****************************************************************************
C**         XM(I),  I=1,NATOM
C**XM:      Nuclear masses (U)
C****************************************************************************
      READ(INP,*)(XM(I),I=1,NATOM)
      DO I=1,NATOM
        WRITE(IOUT,230)I,XM(I)
        XM(I)=XM(I)*ELMASS
      END DO
C**EQUILIBRIUM CARTESIAN COORDINATES
      READ(INP,*)
C****************************************************************************
C**         X0(ATOM,I),   I=1,3
C**X0:      Input NATOM records of equilibrium x,y,z coordinates
C****************************************************************************
      DO I=1,NATOM
        READ(INP,*)(X0(I,J),J=1,3)
      END DO
      WRITE(IOUT,235)
      DO I=1,NATOM
        WRITE(IOUT,240)I,(X0(I,J),J=1,3)
      END DO
C**CHECK C.OF M. AND PRINCIPAL AXES
      IF(MCHECK.GT.0)THEN
        CALL CHECKM(XM,X0,XX,RR,NATOM)
        IF(ISCFCI.LT.0)STOP 'C.M. CHECK'
      END IF
C**NORMAL COORDINATE DISPLACEMENTS
      READ(INP,*)
C****************************************************************************
C**OMIT FOLLOWING INPUT IF INORM = 1
C**         XL(NATOM,NMODE,K),   K=1,3
C**         Mass-weighted Normal Coordinate displacement
C**XL:      Input NMODE groups of NATOM records of x,y,z mass-weighted Normal
C**         Coordinate displacements (INORM=0)
C****************************************************************************
      IF(INORM.EQ.0)THEN
        WRITE(IOUT,245)
        DO J=1,NMODE
          DO I=1,NATOM
            READ(INP,*)(XL(I,J,K),K=1,3)
          END DO
          WRITE(IOUT,250)J
          DO I=1,NATOM
            WRITE(IOUT,255)I,(XL(I,J,K),K=1,3)
          END DO
        END DO
      ELSE
CC      CALL NORMAL(NATOM,NMODE,XM,X0,OMEGA,XL,XX,RR)
CC      WRITE(IOUT,245)
CC      DO J=1,NMODE
CC        WRITE(IOUT,250)J
CC        DO I=1,NATOM
CC          WRITE(IOUT,255)I,(XL(I,J,K),K=1,3)
CC        END DO
CC      END DO
CC      DO K=1,NMODE
CC        DO I=1,6
CC          IPOT(K,I)=0
CC          JPOT(K,I)=0
CC        END DO
CC        IPOT(K,1)=K
CC        JPOT(K,1)=2
CC        CPOT(K)=OMEGA(K)*OMEGA(K)/2
CC        WRITE(IOUT,220)IPOT(K,1),OMEGA(K)*WAVENM
CC      END DO
      END IF
      RETURN
C**********************************************************
C**********************************************************
9999  CONTINUE
      IF(IORDER.NE.0)THEN
C**SET UP NSTAT STATES BASED ON ENERGY
        DO K=1,NMODE
          DO I=1,NSTAT
            ISTAT(I,K)=0
          END DO
        END DO
C**START STACK OFF WITH PURE MODE 1 (LOWEST ENERGY)
        ESTAT(1)=0
        DO I=2,NSTAT
          J=I-1
          ESTAT(I)=OMEGA(1)*J
          ISTAT(I,1)=J
        END DO
C**ADD REMAINING MODES
        DO K=2,NMODE
C**SAVE ISTAT AND ESTAT IN JSTAT AND WSTAT
          DO I=1,NSTAT
            WSTAT(I)=ESTAT(I)
            DO L=1,NMODE
              JSTAT(I,L)=ISTAT(I,L)
            END DO
          END DO
C**FOR MODE K, ADD N1 QUANTA TO ORIGINAL STATE I1 IN ISTAT
          MORES=0
          I1=0
1000      CONTINUE
          I1=I1+1
          MOREQ=0
          N1=0
2000      CONTINUE
          N1=N1+1
C**KSTAT CONTAINS NEW STATE, ETEMP NEW ENERGY
          DO L=1,NMODE
            KSTAT(L)=ISTAT(I1,L)
          END DO
          KSTAT(K)=ISTAT(I1,K)+N1
          ETEMP=0
          DO L=1,NMODE
            ETEMP=ETEMP+KSTAT(L)*OMEGA(L)
          END DO
C**SEE IF NEW STATE REQUIRED
          DO I=1,NSTAT
            IF(WSTAT(I).GT.ETEMP)THEN
C**NEED IT...INSERT IT AND MOVE THE REST UP
C**FIRST MOVE THEM
              I2=NSTAT-I
              DO J=1,I2
                WSTAT(NSTAT+1-J)=WSTAT(NSTAT-J)
                DO L=1,NMODE
                  JSTAT(NSTAT+1-J,L)=JSTAT(NSTAT-J,L)
                END DO
              END DO
C**THEN INSERT
              WSTAT(I)=ETEMP
              DO L=1,NMODE
                JSTAT(I,L)=KSTAT(L)
              END DO
              GO TO 3000
            ELSE
C**MOREQ SET IF CAN'T INSERT NEW STATE
              IF(I.EQ.NSTAT)THEN
                MOREQ=1
                GO TO 3000
              END IF
            END IF
          END DO
3000      CONTINUE
C**TEST IF POSSIBLY ONE MORE QUANTUM
          IF(MOREQ.EQ.0)GO TO 2000
C**TEST IF POSSIBLY ONE MORE STATE
          IF(I1.EQ.NSTAT)MORES=1
          IF(MORES.EQ.0)GO TO 1000
C**MOVE NEW STATES BACK TO ISTAT AND ESTAT
          DO I=1,NSTAT
            ESTAT(I)=WSTAT(I)
            DO L=1,NMODE
              ISTAT(I,L)=JSTAT(I,L)
            END DO
          END DO
        END DO
      WRITE(IOUT,*)
      WRITE(IOUT,*)'ENERGIES AND COMPUTED STATES'
      DO I=1,NSTAT
        WRITE(IOUT,*)ESTAT(I),(ISTAT(I,K),K=1,NMODE)
      END DO
      END IF
C**********************************************************
C**********************************************************
C**RESET QUANTA TO FUNCTION POINTERS
      DO I=1,NSTAT
        DO K=1,NMODE
          ISTAT(I,K)=ISTAT(I,K)+1
        END DO
      END DO
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE CHECKM(XM,X0,XX,RR,NATOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XM(NATOM),X0(NATOM,3),XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION XCM(3),WR(3),E(3)
      COMMON/FILASS/IOUT,INP
      COMMON/VMIN/VMIN
      COMMON/RETURN/IRET
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/MOMI/XK(3,3),XMU(3,3)
100   FORMAT(//,1X,'CENTRE OF MASS CONDITIONS: X,Y,Z ',//,1X,3D20.12,/)
150   FORMAT(//,1X,'PRINCIPAL AXIS CONDITIONS',/)
200   FORMAT(3D20.12)
250   FORMAT(//,1X,'VECTORS',/)
300   FORMAT(//,1X,'C.OF.M COORDINATES: X,Y,Z ',/)
350   FORMAT(//,1X,'ORIGINAL COORDINATES: X,Y,Z ',/)
400   FORMAT(//,1X,'PRINCIPAL AXIS COORDINATES: X,Y,Z ',/)
450   FORMAT(//,1X,'ROTATIONAL CONSTANTS',/)
500   FORMAT(//,1X,'MOMENT OF INERTIA TENSOR',/)
550   FORMAT(//,1X,'CHECK PRINCIPAL AXES COORDINATES',//)
600   FORMAT(//,1X,'ENERGY (Eh) AT MINIMUM ',F12.6,/)
650   FORMAT(//,1X,'ROTATED VECTORS',/)
      WRITE(IOUT,550)
C**RETURN WITH ROTATIONAL CONSTANTS
C**ORIGINAL COORDINATES
      WRITE(IOUT,350)
      DO I=1,NATOM
        WRITE(IOUT,200)(X0(I,J),J=1,3)
      END DO
      CALL ROTC(X0,XM,NATOM,A,B,C,BBAR,1)
      WRITE(IOUT,450)
      WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
CCCC  CALL GETPOT(VMIN,NATOM,X0,RR)
CCCC  WRITE(IOUT,600)VMIN
C**FIND CENTRE OF MASS
      DO I=1,3
        XCM(I)=0
      END DO
      XMASS=0
      DO I=1,NATOM
        XMASS=XMASS+XM(I)
      END DO
      DO I=1,3
        DO J=1,NATOM
          XCM(I)=XCM(I)+XM(J)*X0(J,I)
        END DO
        XCM(I)=XCM(I)/XMASS
      END DO
      DO I=1,3
        DO J=1,NATOM
          XX(J,I)=X0(J,I)-XCM(I)
C**TEMPORARY
          IF(DABS(XX(J,I)).LT.1.D-10)XX(J,I)=0.D0
C**TEMPORARY
        END DO
      END DO
C**COORDINATES RELATIVE TO C. OF M.
      WRITE(IOUT,300)
      DO I=1,NATOM
        WRITE(IOUT,200)(XX(I,J),J=1,3)
      END DO
      CALL ROTC(XX,XM,NATOM,A,B,C,BBAR,1)
      WRITE(IOUT,450)
      WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
C**SET UP MOMENT OF INERTIA MATRIX
      DO IX=1,3
        DO IY=1,IX
          XK(IY,IX)=0
        END DO
      END DO
      DO IX=1,3
        DO I=1,NATOM
          X=0
          DO J=1,3
            IF(J.NE.IX)X=X+XX(I,J)*XX(I,J)
          END DO
          XK(IX,IX)=XK(IX,IX)+XM(I)*X
        END DO
      END DO
      DO IX=1,3
        DO IY=1,IX-1
          DO I=1,NATOM
            XK(IY,IX)=XK(IY,IX)-XM(I)*XX(I,IY)*XX(I,IX)
          END DO
        END DO
      END DO
      DO IX=1,3
        DO IY=1,IX
          XK(IX,IY)=XK(IY,IX)
        END DO
      END DO
C**MOMENT OF INERTIA TENSOR
      WRITE(IOUT,500)
      DO I=1,3
        WRITE(IOUT,200)(XK(J,I),J=1,3)
      END DO
C**DIAGONALISE MOMENT OF INERTIA MATRIX
      IRET=1
      CALL DIAG(XK,XK,3,3,-1,WR,E,WR,3,3,XK,3,3,XK,XK,XK,IDUM,IDUM,
     1IDUM)
      IRET=0
      A=1/(2*E(1))
      B=1/(2*E(2))
      C=1/(2*E(3))
      WRITE(IOUT,450)
      WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
C**WRITE VECTORS
      WRITE(IOUT,250)
      DO I=1,3
        WRITE(IOUT,200)(XK(J,I),J=1,3)
      END DO
C*******************************************
C*******************************************
C**CALL USER-SUPPLIED ROUTINE TO GET 'CLEAN'
C**GEOMETRY IF SYMMETRIC OR SPHERICAL TOP
      CALL RTGEOM(NATOM,XX,E,WAVENM)
C*******************************************
C*******************************************
C**WRITE ROTATED VECTORS
      WRITE(IOUT,650)
      DO I=1,3
        WRITE(IOUT,200)(XK(J,I),J=1,3)
      END DO
C**MOVE TO PRINCIPAL AXES
      DO I=1,NATOM
        DO J=1,3
          X0(I,J)=0
          DO K=1,3
            X0(I,J)=X0(I,J)+XK(K,J)*XX(I,K)
          END DO
        END DO
      END DO
C**MAKE GEOMETRY NICE !
C**TEMPORARY
C**FIRST ZERO COORDINATES
C     DO I=1,NATOM
C       DO J=1,3
C         IF(DABS(X0(I,J)).LT.1.D-8)X0(I,J)=0
C       END DO
C     END DO
C**x,y,z COORDINATES FOR SAME ATOM
C     DO I=1,NATOM
C       DO J1=1,3
C         IS1=1
C         IF(X0(I,J1).LT.0)IS1=-1
C         DO J2=1,J1-1
C           IF(DABS(DABS(X0(I,J1))-DABS(X0(I,J2))).LT.1.D-8)THEN
C             X0(I,J1)=DABS(X0(I,J2))*IS1
C           END IF
C         END DO
C       END DO
C     END DO
C**x,y,z COORDINATES FOR DIFFERENT ATOMS
C     DO I1=1,NATOM
C       DO J1=1,3
C         IS1=1
C         IF(X0(I1,J1).LT.0)IS1=-1
C         DO I2=1,I1-1
C           DO J2=1,3
C             IS2=1
C             IF(X0(I2,J2).LT.0)IS2=-1
C             IF(DABS(DABS(X0(I1,J1))-DABS(X0(I2,J2))).LT.1.D-8)THEN
C               X0(I1,J1)=DABS(X0(I2,J2))*IS1
C             END IF
C           END DO
C         END DO
C       END DO
C     END DO
C**TEMPORARY
C**ROTATED (PRINCIPAL AXIS) COORDINATES
      WRITE(IOUT,400)
      DO I=1,NATOM
        WRITE(IOUT,200)(X0(I,J),J=1,3)
      END DO
C****************************************************
C**CHECK RESULTS
C****************************************************
      DO I=1,3
        XCM(I)=0
      END DO
      DO I=1,3
        DO J=1,NATOM
          XCM(I)=XCM(I)+XM(J)*X0(J,I)
        END DO
      END DO
      WRITE(IOUT,100)(XCM(I),I=1,3)
C**SET UP MOMENT OF INERTIA MATRIX
      DO IX=1,3
        DO IY=1,IX
          XK(IY,IX)=0
        END DO
      END DO
      DO IX=1,3
        DO I=1,NATOM
          X=0
          DO J=1,3
            IF(J.NE.IX)X=X+X0(I,J)*X0(I,J)
          END DO
          XK(IX,IX)=XK(IX,IX)+XM(I)*X
        END DO
      END DO
      DO IX=1,3
        DO IY=1,IX-1
          DO I=1,NATOM
            XK(IY,IX)=XK(IY,IX)-XM(I)*X0(I,IY)*X0(I,IX)
          END DO
        END DO
      END DO
      DO IX=1,3
        DO IY=1,IX
          XK(IX,IY)=XK(IY,IX)
        END DO
      END DO
C**MOMENT OF INERTIA TENSOR
      WRITE(IOUT,150)
      DO I=1,3
        WRITE(IOUT,200)(XK(J,I),J=1,3)
      END DO
      A=1/(2*XK(1,1))
      B=1/(2*XK(2,2))
      C=1/(2*XK(3,3))
      WRITE(IOUT,450)
      WRITE(IOUT,200)A*WAVENM,B*WAVENM,C*WAVENM
CCCC  CALL GETPOT(VMIN,NATOM,X0,RR)
CCCC  WRITE(IOUT,600)VMIN
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE PAXES(XM,X0,NATOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C**TEMPORARY...LINEAR
      LOGICAL LINEAR
C**TEMPORARY...LINEAR
      DIMENSION XM(NATOM),X0(NATOM,3)
      DIMENSION XK(3,3)
      COMMON/AXES/MX(3)
C**TEMPORARY...LINEAR
      COMMON/TYPE/LINEAR
C**TEMPORARY...LINEAR
      COMMON/FILASS/IOUT,INP
200   FORMAT(3D20.12)
250   FORMAT(//,1X,'PAXES - MOMENT OF INERTIA TENSOR',/)
300   FORMAT(/,1X,'PRINCIPAL X,Y,Z AXES CORRESPOND TO INPUT AXES ',3I3,
     1/)
C**TEMPORARY...LINEAR
C**SEARCH FOR LINEARITY ALONG 'Z'
      LINEAR=.TRUE.
      DO I=1,NATOM
C**LOOK AT 'X' AND 'Y'
        DO J=1,2
          IF(X0(I,J).NE.0)LINEAR=.FALSE.
        END DO
      END DO
C**TEMPORARY...LINEAR
C**SET UP MOMENT OF INERTIA MATRIX
      DO IX=1,3
        DO IY=1,IX
          XK(IY,IX)=0
        END DO
      END DO
      DO IX=1,3
        DO I=1,NATOM
          X=0
          DO J=1,3
            IF(J.NE.IX)X=X+X0(I,J)*X0(I,J)
          END DO
          XK(IX,IX)=XK(IX,IX)+XM(I)*X
        END DO
      END DO
      DO IX=1,3
        DO IY=1,IX
          XK(IX,IY)=XK(IY,IX)
        END DO
      END DO
C**MOMENT OF INERTIA TENSOR
      WRITE(IOUT,250)
      DO I=1,3
        WRITE(IOUT,200)(XK(J,I),J=1,3)
      END DO
C**FIND PRINCIPAL 'X' AXIS
      MX(1)=1
      BIG=XK(1,1)
      IF(XK(2,2).GT.BIG)THEN
        BIG=XK(2,2)
        MX(1)=2
      END IF
      IF(XK(3,3).GT.BIG)THEN
        MX(1)=3
      END IF
C**FIND PRINCIPAL 'Z' AXIS
      MX(3)=3
      SMALL=XK(3,3)
      IF(XK(2,2).LT.SMALL)THEN
        SMALL=XK(2,2)
        MX(3)=2
      END IF
      IF(XK(1,1).LT.SMALL)THEN
        MX(3)=1
      END IF
C**FIND PRINCIPAL 'Y' AXIS
      IF(MX(1).EQ.1)THEN
        IF(MX(3).EQ.2)MX(2)=3
        IF(MX(3).EQ.3)MX(2)=2
      END IF
      IF(MX(1).EQ.2)THEN
        IF(MX(3).EQ.1)MX(2)=3
        IF(MX(3).EQ.3)MX(2)=1
      END IF
      IF(MX(1).EQ.3)THEN
        IF(MX(3).EQ.1)MX(2)=2
        IF(MX(3).EQ.2)MX(2)=1
      END IF
C**TEMPORARY
      MX(1)=1
      MX(2)=2
      MX(3)=3
C**TEMPORARY
      WRITE(IOUT,300)MX(1),MX(2),MX(3)
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HERMPT(NV,M,H,XQ,XJ,XV,OMEGA,XQJ,MODE,NMODE,NATOM,QQ,
     1RR,XX,X0,XL,XM,NPOT,IPOT,JPOT,CPOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(NV,M,3),XQ(M),XJ(M),XV(M)
      DIMENSION XQJ(M,M)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/VMIN/VMIN
      COMMON/WHICH/IWHICH
      COMMON/PRINT/IPRINT,JPRINT
C**************************************************************
201   FORMAT(//,1X,'INTEGRATION POINTS AND WEIGHTS FOR MODE ',I2,/)
202   FORMAT(13X,'Q(A0)',11X,'WEIGHT(EXP(-Q**2).D(Q))',3X,
     1'POTENTIAL(CM-1)',/)
203   FORMAT(2X,F20.12,2X,D25.12,2X,F15.2)
C*****************************************************************
      PIQ=1.33133536380038D0
      YLAM=1/SQRT(OMEGA)
      DO 10 J=1,M
      DO 10 I=1,M
10    XQJ(I,J)=0.D0
      M1=M-1
      DO 20 J=1,M1
      I=J+1
      XQJ(J,I)=SQRT(DFLOAT(J)/2.D0)
      XQJ(I,J)=XQJ(J,I)
20    CONTINUE
30    CALL DIAG(XQJ,XQJ,M,M,-1,XJ,XQ,XJ,M,M,XQJ,M,M,XQJ,XQJ,XQJ,IDUM,
     1IDUM,IDUM)
      IF(IPRINT.GT.0)WRITE(IOUT,201)MODE
      IF (IPRINT.GT.0) WRITE(IOUT,202)
      K=0
      DO 60 J=1,M
      XJ(J)=PIQ*PIQ*XQJ(1,J)*XQJ(1,J)
      Q=XQ(J)
      XQ(J)=Q*YLAM
C**ENERGY THIS POINT
      DO K=1,NMODE
        QQ(K)=0
      END DO
      QQ(MODE)=XQ(J)
      IF(IWHICH.NE.0)THEN
          DO I=1,NATOM
            DO K=1,3
              XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1        SQRT(XM(I))
            END DO
          END DO
          CALL GETPOT(VPT,NATOM,XX,RR)
        ELSE
          CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,NMODE,QQ)
      END IF
      XV(J)=VPT
      IF (IPRINT.GT.0) WRITE(IOUT,203)XQ(J),XJ(J),(VPT-VMIN)*WAVENM
C**ENERGY THIS POINT
60    CONTINUE
      DO 100 J=1,M
      Q=(XQ(J))/YLAM
      H(1,J,1)=1.D0/PIQ
      H(2,J,1)=2.D0*Q/(PIQ*SQRT(2.D0))
      IF(NV.GE.3)THEN
      DO 70 I=3,NV
      I1=I-1
      I2=I-2
70    H(I,J,1)=2*(Q*H(I1,J,1)-I2*H(I2,J,1)/SQRT(2.D0*I2))/SQRT(2.D0*I1)
      END IF
      H(1,J,2)=0.D0
      H(2,J,2)=2.D0/(PIQ*SQRT(2.D0))
      H(1,J,3)=0.D0
      H(2,J,3)=0.D0
      IF(NV.GE.3)THEN
      DO 80 I=3,NV
      I1=I-1
      I2=I-2
      H(I,J,2)=2.D0*I1*H(I1,J,1)/SQRT(2.D0*I1)
      H(I,J,3)=4.D0*I1*I2*H(I2,J,1)/SQRT(4.D0*I1*I2)
80    CONTINUE
      END IF
      DO 90 I=1,NV
      H(I,J,3)=H(I,J,3)-2*Q*H(I,J,2)+(Q*Q-1.D0)
     1*H(I,J,1)
90    H(I,J,2)=H(I,J,2)-Q*H(I,J,1)
100   CONTINUE
      DO 110 J=1,M
      XJS=SQRT(XJ(J))
C**NEXT STATEMENT ONLY NEEDED IF RESTORING FOR HEG
      XJ(J)=XQ(J)
      DO 110 I=1,NV
      H(I,J,1)=H(I,J,1)*XJS
      H(I,J,2)=H(I,J,2)*XJS/YLAM
      H(I,J,3)=H(I,J,3)*XJS/(YLAM*YLAM)
110   CONTINUE
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HEG0(NB,M,MODE,NVB,MVB,NVF,NMODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NVF(NMODE)
      COMMON/FILASS/IOUT
      COMMON/PRINT/IPRINT,JPRINT
C**DETERMINE NEW NUMBER OF POINTS FROM NB,M,NVF
      MVB=NVF(MODE)+M-NB
C**MAKE IT ODD
      IF(MOD(MVB,2).EQ.0)MVB=MVB+1
C**CAN'T BE BIGGER THAT NUMBER OF PRIMITIVES
C**IF SO, RETAIN M GAUSS POINTS
      IF(MVB.GT.NB)MVB=M
C**THERE ARE M ORIGINAL GAUSS POINTS
C**THERE ARE NB ORIGINAL PRIMITIVES
C**THERE ARE MVB CONTRACTED FUNCTIONS (AND HEG POINTS)
C**WE USE NVB CONTRACTED FUNCTIONS
      NVB=NVF(MODE)
      IF(IPRINT.GT.0)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)' NO. GAUSS POINTS',M
        WRITE(IOUT,*)' NO. PRIMITIVES',NB
        WRITE(IOUT,*)' NO. HEG POINTS',MVB
        WRITE(IOUT,*)' NO. CONTRACTED FUNCTIONS USED',NVB
        IF(MVB.EQ.M)THEN
          WRITE(IOUT,*)
          WRITE(IOUT,*)' *********RETAIN GAUSS INTEGRATION*********'
        END IF
      END IF
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HEG1(NB,M,H,XQ,XJ,OMEGA,XQJ,MODE,HX,HY,NVB,MVB,
     1XV,XQX,MBF,NBF,NMODE,YK,SUP4,NATOM,QQ,
     2RR,XX,X0,XL,XM,NPOT,IPOT,JPOT,CPOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(NB,M,3),XQ(M),XJ(M),XV(M),MBF(NMODE),NBF(NMODE)
      DIMENSION HX(NB,NB,3),HY(NVB,MVB,3),XQX(NB)
      DIMENSION XQJ(NB,NB),YK(NB,NB),SUP4(NB)
      DIMENSION RR(NATOM,NATOM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      DIMENSION IPOT(NPOT,6),JPOT(NPOT,6),CPOT(NPOT)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/VMIN/VMIN
      COMMON/WHICH/IWHICH
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/PATH/ISCFCI
      COMMON/CYCLE/ICYCLE
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/DUMP/JDUMP,IDUMP
C**************************************************************
200   FORMAT(1X,'IFAIL IN F04AAF = ',I3)
201   FORMAT(//,1X,'HEG INTEGRATION POINTS FOR MODE ',I2,/)
202   FORMAT(13X,'Q(A0)',11X,'POTENTIAL(CM-1)',/)
203   FORMAT(2X,F20.12,2X,F15.2)
C*****************************************************************
      PIQ=1.33133536380038D0
      YLAM=1/SQRT(OMEGA)
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)WRITE(60)YLAM,MVB
      IF(MVB.NE.M)THEN
C**GET HEG POINTS IF SUFFICIENT PRIMITIVES
        DO J=1,MVB
          DO I=1,MVB
            XQJ(I,J)=0
          END DO
        END DO
C**GET EXPECTATION VALUE OF Q FOR MODE 'MODE' WITH GAUSS POINTS
        DO MM=1,M
          Q=XQ(MM)/YLAM
          DO IX=1,MVB
            X=H(IX,MM,1)*Q
            DO IY=1,MVB
              XQJ(IY,IX)=XQJ(IY,IX)+H(IY,MM,1)*X
            END DO
          END DO
        END DO
        CALL DIAG(XQJ,XQJ,NB,MVB,-1,XJ,XQX,XJ,MVB,MVB,XQJ,MVB,MVB,XQJ,
     1  XQJ,XQJ,IDUM,IDUM,IDUM)
        IF(IPRINT.GT.0)WRITE(IOUT,201)MODE
        IF(IPRINT.GT.0)WRITE(IOUT,202)
        DO J=1,MVB
          Q=XQX(J)
          XQX(J)=Q*YLAM
C**ENERGY THIS POINT
      DO K=1,NMODE
        QQ(K)=0
      END DO
      QQ(MODE)=XQX(J)
      IF(IWHICH.NE.0)THEN
          DO I=1,NATOM
            DO K=1,3
              XX(I,K)=X0(I,K)+XL(I,MODE,K)*QQ(MODE)/
     1        SQRT(XM(I))
            END DO
          END DO
          CALL GETPOT(VPT,NATOM,XX,RR)
        ELSE
          CALL GETPT1(VPT,NPOT,IPOT,JPOT,CPOT,NMODE,QQ)
      END IF
      XV(J)=VPT
      IF (IPRINT.GT.0) WRITE(IOUT,203)XQX(J),(VPT-VMIN)*WAVENM
C**ENERGY THIS POINT
        END DO
C**GET PRIMITIVES (NB PRIMITIVES) AT HEG POINTS (MVB POINTS)
        DO 100 J=1,MVB
        Q=(XQX(J))/YLAM
        HX(1,J,1)=1.D0/PIQ
        HX(2,J,1)=2.D0*Q/(PIQ*SQRT(2.D0))
        IF(NB.GE.3)THEN
        DO 70 I=3,NB
        I1=I-1
        I2=I-2
70      HX(I,J,1)=2*(Q*HX(I1,J,1)-I2*HX(I2,J,1)/SQRT(2.D0*I2))/
     1  SQRT(2.D0*I1)
        END IF
        HX(1,J,2)=0.D0
        HX(2,J,2)=2.D0/(PIQ*SQRT(2.D0))
        HX(1,J,3)=0.D0
        HX(2,J,3)=0.D0
        IF(NB.GE.3)THEN
        DO 80 I=3,NB
        I1=I-1
        I2=I-2
        HX(I,J,2)=2.D0*I1*HX(I1,J,1)/SQRT(2.D0*I1)
        HX(I,J,3)=4.D0*I1*I2*HX(I2,J,1)/SQRT(4.D0*I1*I2)
80      CONTINUE
        END IF
        DO 90 I=1,NB
        HX(I,J,3)=HX(I,J,3)-2*Q*HX(I,J,2)+(Q*Q-1.D0)
     1  *HX(I,J,1)
90      HX(I,J,2)=HX(I,J,2)-Q*HX(I,J,1)
100     CONTINUE
        DO 110 J=1,MVB
        DO 110 I=1,NB
        HX(I,J,2)=HX(I,J,2)/YLAM
        HX(I,J,3)=HX(I,J,3)/(YLAM*YLAM)
110     CONTINUE
C**FORM NEW CONTRACTED FUNCTIONS AT HEG PTS IN OLD STORAGE LOCATIONS
        DO K=1,3
          DO J=1,MVB
            DO I=1,NVB
              HY(I,J,K)=0
              DO L=1,NB
                HY(I,J,K)=HY(I,J,K)+YK(L,I)*HX(L,J,K)
              END DO
            END DO
          END DO
        END DO
C**WEIGHTS FOR FUNCTION I (FUNCTION*SQRT(WT) IN XQJ)
        DO I=1,NVB
          DO J=1,MVB
            XJ(J)=XQJ(I,J)/HY(I,J,1)
          END DO
          DO K=1,3
            DO J=1,MVB
              HY(I,J,K)=HY(I,J,K)*XJ(J)
            END DO
          END DO
        END DO
      ELSE
C**RETAIN GAUSS POINTS AND GAUSS CONTRACTED FUNCTIONS
        DO K=1,3
          DO J=1,MVB
            DO I=1,NVB
              HY(I,J,K)=H(I,J,K)
            END DO
          END DO
        END DO
      END IF
      IF(ISCFCI.GT.0.AND.ICI.LT.0.AND.IDUMP.NE.0)THEN
        WRITE(60)(XQ(J),J=1,MVB)
      END IF
C**TEMPORARY (TEST ORTHOGONALITY)
C     DO J=1,NVB
C       DO I=1,NVB
C         XQJ(I,J)=0
C       END DO
C     END DO
C     DO MM=1,MVB
C       DO IX=1,NVB
C         X=HY(IX,MM,1)
C         DO IY=1,NVB
C           XQJ(IY,IX)=XQJ(IY,IX)+HY(IY,MM,1)*X
C         END DO
C       END DO
C     END DO
C     DO IX=1,NVB
C       WRITE(IOUT,*)(XQJ(IY,IX),IY=1,NVB)
C     END DO
C**TEMPORARY (TEST ORTHOGONALITY)
C**CHANGE NUMBER OF POINTS ARRAY
      MBF(MODE)=MVB
C**CHANGE NUMBER OF FUNCTIONS ARRAY
      NBF(MODE)=NVB
      RETURN
      END
C**************************************************************
C**************************************************************
      SUBROUTINE HEG2(MAXBAS,NMODE,J,NVF,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MAXBAS(NMODE,J),NVF(NMODE)
      NVF(MODE)=0
      DO K=1,J
        IF(MAXBAS(MODE,K).GT.NVF(MODE))NVF(MODE)=MAXBAS(MODE,K)
      END DO
      RETURN
      END
