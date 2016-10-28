C****************************************************************
C****************************************************************
C**LANCZOS
C****************************************************************
C****************************************************************
      SUBROUTINE LANCZA(X,V,Z,W,Y,S,ISIZE,NVAL,ISTAT,NSTAT,NMODE,
     1IASSIG,ISIZMX,J21,NS,ENERGY,ELAST,NVEC,XA,XACOPY,XK,WK,SUP4,EVAL,
     2NCYCLE,WRK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT,LGIV
C**STORE X,V,S,Z,WRK ON DISC
      DIMENSION X(ISIZE),V(ISIZE),Z(ISIZE),S(ISIZE),WRK(ISIZE),
     1W(ISIZE),Y(ISIZE)
      DIMENSION ENERGY(NVAL),ELAST(NVAL),NVEC(NVAL)
      DIMENSION XK(NVAL*NCYCLE,NVAL),XA(1),XACOPY(1)
      DIMENSION WK(NVAL*NCYCLE),SUP4(5*NVAL*NCYCLE),EVAL(NVAL*NCYCLE)
      DIMENSION ISTAT(NSTAT,NMODE),IASSIG(ISIZMX,J21,3,1)
      DIMENSION JJMAX(3),XMAX(4)
      COMMON/FILASS/IOUT,INP
      COMMON/TRIATO/TRIAT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/EVL/EVL,CUT
      COMMON/GIVEN/LGIV
      COMMON/LANTOL/TOLLAN
      COMMON/CIDIAG/ICID
      COMMON/ROTS/JMAX,KMAX,IDUM,KEL21,KEL
      COMMON/MATRIX/NVALV,NVALR,KSTEP
      COMMON/JKAKC/JTHIS,KA,KC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
C****************************************************************
100   FORMAT(/,1X,'INITIAL ENERGY = ',D20.12)
105   FORMAT(/,1X,'LANCZOS NOT CONVERGED FOR FUNCTIONS:')
110   FORMAT(/,1X,'LANCZOS FULLY CONVERGED')
115   FORMAT(20I4)
C****************************************************************
C**STORE V ON 52
      REWIND 52
C**STORE Z ON 53
      REWIND 53
C**STORE X ON 54
      REWIND 54
C**STORE S ON 55
      REWIND 55
C**************************************************************
C**************************************************************
C**INITIAL GUESS !!!
C**************************************************************
C**************************************************************
      DO K=1,NVAL
        DO I=1,ISIZE
          X(I)=0
        END DO
        XLAST=1.D+10
        DO I=1,ISIZE
          IF(Y(I).LT.XLAST)THEN
            IF(K.EQ.1)THEN
              XLAST=Y(I)
              NVEC(K)=I
            ELSE
              IGOT=0
              DO J=1,K-1
                IF(I.EQ.NVEC(J))IGOT=1
              END DO
              IF(IGOT.EQ.0)THEN
                XLAST=Y(I)
                NVEC(K)=I
              END IF
            END IF
          END IF
        END DO
C**INITIAL GUESS !!!
        NV=NVEC(K)
        X(NV)=1
C*******************************DISC(54: X - WRITE)
        WRITE(54)X
      END DO
C**************************************************************
C**************************************************************
C**********INITIAL SET-UP.....M(rs) = <psi(r)/H/psi(s)>
C**ASSUME INITIAL FUNCTIONS X0 ORTHONORMAL (VIRTUAL SCF STATES)
C**SET INITIAL FUNCTION V0 = X0 (NORMALISED)
C**************************************************************
C**************************************************************
      REWIND 54
      DO K=1,NVAL
C*******************************DISC(54: X - READ)
        READ(54)X
        DO I=1,ISIZE
          V(I)=X(I)
        END DO
C*******************************DISC(52: V - WRITE)
        WRITE(52)V
C**FORM Z0 = M.V0 (UNNORMALISED)
        CALL MX(V,Z,W,W,ISIZE)
C*******************************DISC(53: Z - WRITE)
        WRITE(53)Z
C**FORM S0 = Z0
        DO I=1,ISIZE
          S(I)=Z(I)
        END DO
C*******************************DISC(55: S - WRITE)
        WRITE(55)S
C**FORM V0.(M.V0) = V0.Z0
        J0=K*(K-1)/2
        REWIND 53
        DO KZ=1,K
C*******************************DISC(53: Z - READ)
          READ(53)Z
          CALL DOT(V,Z,E,ISIZE)
          XA(J0+KZ)=E
          XACOPY(J0+KZ)=E
        END DO
C**V0.(M.V0) = V0.Z0 = <V0/H/V0> = E-VALUE OF V0
        ENERGY(K)=E*WAVENM
        ELAST(K)=ENERGY(K)
      END DO
C**************************************************************
C**************************************************************
C**FORM NEXT BASIS FUNCTION V1 (UNNORMALISED)
C**************************************************************
C**************************************************************
      KGOT=0
      REWIND 54
      REWIND 55
      DO K=1,NVAL
C*******************************DISC(54: X - READ)
        READ(54)X
C*******************************DISC(55: S - READ)
        READ(55)S
        J0=K*(K-1)/2
        NV=NVEC(K)
        V(NV)=0
        DO I=1,ISIZE
          IF(I.NE.NV)THEN
            V(I)=-(S(I)-XA(J0+K)*X(I))/(Y(I)-XA(J0+K))
          END IF
        END DO
C**ORTHONORMALISE V1 TO V0
        CALL SCHMDA(V,K,WRK,NVAL,ISIZE,1)
C**SCHMIDT POSITIONS DISC 52 FOR NEW V(K)
C*******************************DISC(52: V - WRITE)
        WRITE(52)V
        NVEC(K)=100000
        IGOT=0
        DO I=1,ISIZE
          IF(V(I).NE.0)IGOT=1
        END DO
C******************************
        IF(IGOT.EQ.0)THEN
          ELAST(K)=0.D0
          NVEC(K)=2
        END IF
CSC     IF(IGOT.EQ.0)THEN
CSC       ELAST(K)=0.D0
CSC       NVEC(K)=2
CSC     END IF
C******************************
        IF(ELAST(K).NE.0.D0)KGOT=1
      END DO
      IF(KGOT.EQ.0)GO TO 2000
      KSIZE=NVAL
C**************************************************************
C**************************************************************
C**START LANCZOS ITERATIONS
C**************************************************************
C**************************************************************
      DO ILAN=2,NCYCLE
        ILOFF=0
        DO I=1,ILAN-1
C***********************
C**OFF-DIAGONAL ELEMENTS
C***********************
C**POSITION DISC 52 FOR V
          REWIND 52
          DO J=1,ILAN-1
            DO K=1,NVAL
C*******************************DISC(52: V - READ)
              READ(52)V
            END DO
          END DO
          KVV=0
          DO K=1,NVAL
C*******************************DISC(52: V - READ)
            READ(52)V
            IF(NVEC(K).GT.ILAN)THEN
              KVV=KVV+1
              IROFF=KSIZE+KVV
              J0=ILOFF+IROFF*(IROFF-1)/2
              KZZ=0
C**POSITION DISC 53 FOR Z
              REWIND 53
              DO J=1,I-1
                DO KZ=1,NVAL
C*******************************DISC(53: Z - READ)
                  READ(53)Z
                END DO
              END DO
              DO KZ=1,NVAL
C*******************************DISC(53: Z - READ)
                READ(53)Z
                IF(NVEC(KZ).GT.I)THEN
                  KZZ=KZZ+1
                  CALL DOT(V,Z,P,ISIZE)
                  XA(J0+KZZ)=P
                  XACOPY(J0+KZZ)=P
                END IF
              END DO
            END IF
          END DO
          ILOFF=ILOFF+NVAL
          DO K=1,NVAL
            IF(NVEC(K).LE.I)ILOFF=ILOFF-1
          END DO
        END DO
C*******************
C**DIAGONAL ELEMENTS
C*******************
        KVV=0
C**POSITION DISC 52 FOR V
        REWIND 52
        DO J=1,ILAN-1
          DO K=1,NVAL
C*******************************DISC(52: V - READ)
            READ(52)V
          END DO
        END DO
        DO K=1,NVAL
C*******************************DISC(52: V - READ)
          READ(52)V
C**FORM Z = M.V
          IF(NVEC(K).GT.ILAN)CALL MX(V,Z,W,W,ISIZE)
C**DISC 53 ALREADY POSITIONED FOR Z
C*******************************DISC(53: Z - WRITE)
          WRITE(53)Z
          IF(NVEC(K).GT.ILAN)THEN
            KVV=KVV+1
            IROFF=KSIZE+KVV
            J0=KSIZE+IROFF*(IROFF-1)/2
            KZZ=0
C**POSITION DISC 53 FOR Z
            REWIND 53
            DO J=1,ILAN-1
              DO KZ=1,NVAL
C*******************************DISC(53: Z - READ)
                READ(53)Z
              END DO
            END DO
            DO KZ=1,K
C*******************************DISC(53: Z - READ)
              READ(53)Z
              IF(NVEC(KZ).GT.ILAN)THEN
                KZZ=KZZ+1
                CALL DOT(V,Z,P,ISIZE)
                XA(J0+KZZ)=P
                XACOPY(J0+KZZ)=P
              END IF
            END DO
          END IF
        END DO
C**GET CURRENT MATRIX SIZE FROM PREVIOUS SIZE
        KSIZE=KSIZE+NVAL
        DO K=1,NVAL
          IF(NVEC(K).LE.ILAN)KSIZE=KSIZE-1
        END DO
C**SET UP SQUARE MATRIX FOR QL
        IF(.NOT.LGIV)THEN
          CALL QLCOPY(XACOPY,KSIZE,XK,NVAL*NCYCLE)
        END IF
C**DIAGONALISE MATRIX ORDER KSIZE (NO OVERLAP)
        CALL DIAG(XA,XK,NVAL*NCYCLE,KSIZE,-1,SUP4,EVAL,WK,
     1  NVAL,NVAL,IP,ISIZMX,NMODE,XK,XK,IASSIG,ISIZMX,J21,NS)
C**************************************************************
C**************************************************************
C**FORM NEW FUNCTION X.....LOWEST-ENERGY EIGENFUNCTIONS K (NORMALISED)
C**************************************************************
C**************************************************************
        REWIND 54
        DO K=1,NVAL
C*******************************DISC(54: X - READ)
          READ(54)X
CNCH      IF(ELAST(K).NE.0.D0)THEN
            DO I=1,ISIZE
              X(I)=0
            END DO
C*******************************DISC(54: X - WRITE)
            BACKSPACE 54
            WRITE(54)X
CNCH      END IF
        END DO
        JOFF=0
        DO J=1,ILAN
          REWIND 54
          DO K=1,NVAL
C*******************************DISC(54: X - READ)
            READ(54)X
CNCH        IF(ELAST(K).NE.0.D0)THEN
              KOFF=JOFF
C**POSITION DISC 52 FOR V
              REWIND 52
              DO I=1,J-1
                DO KV=1,NVAL
C*******************************DISC(52: V - READ)
                  READ(52)V
                END DO
              END DO
              DO N=1,NVAL
C*******************************DISC(52: V - READ)
                READ(52)V
                IF(NVEC(N).GT.J)THEN
                  KOFF=KOFF+1
                  DO I=1,ISIZE
                    X(I)=X(I)+XK(KOFF,K)*V(I)
                  END DO
                END IF
              END DO
C*******************************DISC(54: X - WRITE)
              BACKSPACE 54
              WRITE(54)X
CNCH        END IF
          END DO
          JOFF=KOFF
        END DO
C**************************************************************
C**************************************************************
C**FORM NEW FUNCTION S
C**************************************************************
C**************************************************************
        REWIND 55
        DO K=1,NVAL
CNCH      IF(ELAST(K).NE.0.D0)ENERGY(K)=EVAL(K)*WAVENM
          ENERGY(K)=EVAL(K)*WAVENM
          DO I=1,ISIZE
            S(I)=0
          END DO
C*******************************DISC(55: S - WRITE)
          WRITE(55)S
        END DO
        JOFF=0
        DO J=1,ILAN
          REWIND 55
          DO K=1,NVAL
C*******************************DISC(55: S - READ)
            READ(55)S
            KOFF=JOFF
C**POSITION DISC 53 FOR Z
            REWIND 53
            DO I=1,J-1
              DO KV=1,NVAL
C*******************************DISC(53: Z - READ)
                READ(53)Z
              END DO
            END DO
            DO N=1,NVAL
C*******************************DISC(53: Z - READ)
              READ(53)Z
              IF(NVEC(N).GT.J)THEN
                KOFF=KOFF+1
                DO I=1,ISIZE
                  S(I)=S(I)+XK(KOFF,K)*Z(I)
                END DO
              END IF
            END DO
C*******************************DISC(55: S - WRITE)
            BACKSPACE 55
            WRITE(55)S
          END DO
          JOFF=KOFF
        END DO
        IF(ILAN.EQ.NCYCLE)GO TO 1000
        KGOT=0
        DO K=1,NVAL
          IF(ELAST(K).NE.0.D0)THEN
            IF(DABS(ENERGY(K)-ELAST(K)).LT.TOLLAN)THEN
              ELAST(K)=0.D0
C**NVEC(K) IS FIRST CYCLE IN WHICH FUNCTION 'K' IS NOT INCLUDED
              NVEC(K)=ILAN+1
            ELSE
              ELAST(K)=ENERGY(K)
            END IF
          END IF
          IF(ELAST(K).NE.0.D0)KGOT=1
        END DO
        IF(KGOT.EQ.0)GO TO 2000
C**************************************************************
C**************************************************************
C**FORM NEXT BASIS FUNCTION V(ILAN+1) (UNNORMALISED)
C**************************************************************
C**************************************************************
        KGOT=0
        REWIND 54
        REWIND 55
        DO K=1,NVAL
C*******************************DISC(54: X - READ)
          READ(54)X
C*******************************DISC(55: S - READ)
          READ(55)S
          IF(ELAST(K).NE.0.D0)THEN
            DO I=1,ISIZE
              V(I)=-(S(I)-EVAL(K)*X(I))/
     1        (Y(I)-EVAL(K))
            END DO
C**ORTHONORMALISE V(ILAN+1) TO V(ILAN)
            CALL SCHMDA(V,K,WRK,NVAL,ISIZE,ILAN)
            IGOT=0
            DO I=1,ISIZE
              IF(V(I).NE.0)IGOT=1
            END DO
C******************************************
            IF(IGOT.EQ.0)THEN
              ELAST(K)=0.D0
              NVEC(K)=ILAN+1
            END IF
CSC         IF(IGOT.EQ.0)THEN
CSC           ELAST(K)=0.D0
CSC           NVEC(K)=ILAN+1
CSC         END IF
C******************************************
          ELSE
            DO I=1,ISIZE
              V(I)=0
            END DO
          END IF
C**SCHMIDT POSITIONS DISC 52 FOR NEW V(K)
C**OTHERWISE DISC 52 ALREADY POSITIONED FORMING NEW X
C*******************************DISC(52: V - WRITE)
          WRITE(52)V
          IF(ELAST(K).NE.0.D0)KGOT=1
        END DO
        IF(KGOT.EQ.0)GO TO 2000
C**RE-LOAD MATRIX ORDER ILAN IN PREPARATION FOR MATRIX ORDER ILAN+1
        DO I=1,KSIZE*(KSIZE+1)/2
          XA(I)=XACOPY(I)
        END DO
      END DO
C**************************************************************
C**************************************************************
C**END LANCZOS ITERATIONS
C**************************************************************
C**************************************************************
1000  CONTINUE
      WRITE(IOUT,105)
      N=0
      DO K=1,NVAL
        IF(ELAST(K).NE.0)THEN
          N=N+1
          NVEC(N)=K
        END IF
      END DO
      WRITE(IOUT,115)(NVEC(K),K=1,N)
      GO TO 3000
2000  CONTINUE
      WRITE(IOUT,110)
3000  CONTINUE
C**************************************************************
C**************************************************************
C**ASSIGN LEVELS
C**************************************************************
C**************************************************************
      MT=1
      IF(EVL.NE.0.0D0)GO TO 4
      MT=2
4     CONTINUE
      IF(EVL.EQ.0)EVL=ENERGY(1)
      DO 7 I=1,NVAL
      ELAST(I)=ENERGY(I)/WAVENM
      ENERGY(I)=ENERGY(I)-EVL
7     CONTINUE
      XMAX(1)=2
      IF(ICID.EQ.0)THEN
        NSIZE=ISIZE
        KKC=KC
        KKA=KA
      ELSE
CC      J21=2*JTHIS+1
        NSIZE=NVALV
      END IF
      REWIND 54
      ICYCL=1
      IF(JPRINT.LT.-1.AND.ICID.EQ.0)ICYCL=2
      IF(JPRINT.LT.-2.AND.ICID.EQ.0)ICYCL=3
      DO 10 I=1,NVAL
C*******************************DISC(54: X - READ)
      READ(54)X
      IF(ENERGY(I).GT.CUT)GO TO 1001
      DO 101 ICYC=1,ICYCL
      NSX=NS
      IOFF=NSX-NVSYM
      IF(IOFF.GT.0)NSX=IOFF
      IF(ICID.NE.0)THEN
        KKC=JTHIS
        KKA=0
        IF(TRIAT.AND.NS.GT.NVSYM)THEN
          KKC=JTHIS-1
          KKA=1
        END IF
      END IF
      JX=I+1-MT
      if (icyc.lt.3) then
c  find two max ci coeffs
        XMAX(ICYC+1)=0.D0
        DO 9 J=1,ISIZE
        IF(X(J).EQ.XMAX(ICYC))GO TO 9
        IF(ABS(X(J)).LE.ABS(XMAX(ICYC+1)))GO TO 9
        XMAX(ICYC+1)=X(J)
        IMAX=J
9       CONTINUE
      else
c  get third biggest coeff
        xmax(4)=0.0
        do 90 j=1,isize
        if((x(j).eq.xmax(2)).or.(x(j).eq.xmax(3))) go to 90
        if (abs(x(j)).lt.abs(xmax(4))) go to 90
        xmax(4)=x(j)
        imax=j
90      continue
      endif
      IIA=0
12    CONTINUE
      IF(NSIZE.GE.IMAX)GO TO 13
      IMAX=IMAX-NSIZE
      IIA=IIA+1
      IF(TRIAT)THEN
        KKA=KKA+1
        IF(MOD(IIA,2).EQ.0)KKC=KKC-2
      ELSE
        IA=MOD(IIA,2)
        KKA=KKA+IA
        KKC=JTHIS-KKA
        KKC=KKC+MOD(IIA+1,2)
      END IF
C**UPDATE VIBRATIONAL SYMMETRY FOR NEXT TIME IF REQUIRED
      IF((.NOT.TRIAT).AND.NVSYM.EQ.4)THEN
        IF(MOD(IIA,2).NE.0)THEN
          NSX=NSX+2
        ELSE
          IF(MOD(IIA,4).NE.0)THEN
            NSX=NSX+(-1)**NS
          ELSE
            NSX=NSX-(-1)**NS
          END IF
        END IF
        IF(NSX.GT.NVSYM)THEN
          INCR=NSX-NVSYM
          NSX=INCR
        END IF
        IF(NSX.EQ.0)NSX=NVSYM
      ELSE
        NSX=NSX+1
        IF(NSX.GT.NVSYM)NSX=1
      END IF
      GO TO 12
13    CONTINUE
      IF(ICID.EQ.0)THEN
C**STORE CI ASSIGNMENTS
        KROT=KEL
        IASSIG(I,KROT,ICYC,NSX)=IMAX
        JJMAX(ICYC)=IMAX
      ELSE
C**RECALL CI ASSIGNMENTS
        IOFF=0
        IF(NSX.GT.1)THEN
          DO K=1,NSX-1
            IOFF=IOFF+NTOT(K)
          END DO
        END IF
        IF(TRIAT)THEN
          KROT=2*KKA
          KOFF=1
          IF(NS.GT.NVSYM)KOFF=0
          JROT=KROT+KOFF
        ELSE
          KROT=2*KKA
          IF(KKA.EQ.0.OR.KKA+KKC.NE.JTHIS)KROT=KROT+1
          JROT=KROT
        END IF
        KKEL=(JROT-1)/KSTEP+1
        JJMAX(1)=IOFF+IASSIG(IMAX,KKEL,1,NSX)
        IF(JPRINT.LT.-1)JJMAX(2)=IOFF+IASSIG(IMAX,KKEL,2,NSX)
        IF(JPRINT.LT.-2)JJMAX(3)=IOFF+IASSIG(IMAX,KKEL,3,NSX)
      END IF
101   CONTINUE
      CALL PRCI(ISTAT,NSTAT,NMODE,JJMAX(1),JJMAX(2),JJMAX(3),
     1XMAX(2),XMAX(3),XMAX(4),JX,ENERGY(I),EVL,JTHIS,KKA,KKC,45)
1001  CONTINUE
10    CONTINUE
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE SCHMDA(V,KV,WRK,NVAL,ISIZE,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(ISIZE),WRK(ISIZE)
      COMMON/FILASS/IOUT,INP
C**SCHMIDT ORTHOGONALISE
      N=NV+1
C**NEW FUNCTION IS V(ISIZE,KV,N)
C**NORMALISE NEW FUNCTION
      CALL DOT(V,V,S,ISIZE)
      SQ=SQRT(S)
C**BEWARE!!! IS 10**-6 CORRECT?
C**WE WILL MAKE THE PROGRAM PRINT...
C**EACH TIME THIS IS OBEYED
      IF(SQ.GT.1.D-6)THEN
        S=1/SQ
      ELSE
        WRITE(IOUT,*)'ZERO VECTOR NO. ',KV
        S=0
C**POSITION DISC 52 FOR V
        REWIND 52
        DO I=1,N
          NXVAL=NVAL
          IF(I.EQ.N)NXVAL=KV-1
          DO K=1,NXVAL
            READ(52)WRK
          END DO
        END DO
      END IF
      DO I=1,ISIZE
        V(I)=V(I)*S
      END DO
      IF(S.EQ.0)RETURN
C**IF KV=1, OLD FUNCTIONS ARE V(ISIZE,NVAL,NV)
C**IF KV>1, OLD FUNCTIONS ALSO INCLUDE V(ISIZE,KV-1,N)
      REWIND 52
      DO I=1,N
        NXVAL=NVAL
        IF(I.EQ.N)NXVAL=KV
        DO K=1,NXVAL
C*******************************DISC(52: V - READ (USE WRK))
          IF(I.NE.N.OR.K.NE.KV)READ(52)WRK
          IF(I.EQ.N.AND.K.EQ.KV)THEN
            DO J=1,ISIZE
              WRK(J)=V(J)
            END DO
          END IF
          CALL DOT(V,WRK,S,ISIZE)
          IF(I.EQ.N.AND.K.EQ.KV)THEN
            SQ=SQRT(S)
C**BEWARE!!! IS 10**-6 CORRECT?
C**WE WILL MAKE THE PROGRAM PRINT...
C**EACH TIME THIS IS OBEYED
            IF(SQ.GT.1.D-6)THEN
              S=1.D0-1.D0/SQ
            ELSE
              WRITE(IOUT,*)'REJECTED VECTOR NO. ',K
              CALL FLUSH(IOUT)
              DO J=1,ISIZE
                V(J)=0
                WRK(J)=0
              END DO
            END IF
          END IF
          DO J=1,ISIZE
            V(J)=V(J)-WRK(J)*S
          END DO
        END DO
      END DO
C**THE ABOVE IS SCHMIDT OF W.M.
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE LANCZB(X,V,Z,W,Y,S,ISIZE,NVAL,ISTAT,NSTAT,NMODE,
     1IASSIG,ISIZMX,J21,NS,ENERGY,ELAST,NVEC,XA,XACOPY,XK,WK,SUP4,EVAL,
     2NCYCLE,WRK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL TRIAT,LGIV
C**STORE X,V,Z ON DISC
      DIMENSION X(ISIZE,NVAL),V(ISIZE,NVAL),Z(ISIZE,NVAL),
     1WRK(ISIZE,NVAL),
     2W(ISIZE),Y(ISIZE),S(ISIZE,NVAL)
      DIMENSION ENERGY(NVAL),ELAST(NVAL),NVEC(NVAL)
      DIMENSION XK(NVAL*NCYCLE,NVAL),XA(1),XACOPY(1)
      DIMENSION WK(NVAL*NCYCLE),SUP4(5*NVAL*NCYCLE),EVAL(NVAL*NCYCLE)
      DIMENSION ISTAT(NSTAT,NMODE),IASSIG(ISIZMX,J21,3,1)
      DIMENSION JJMAX(3),XMAX(4)
      COMMON/FILASS/IOUT,INP
      COMMON/TRIATO/TRIAT
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/PRINT/IPRINT,JPRINT
      COMMON/EVL/EVL,CUT
      COMMON/GIVEN/LGIV
      COMMON/LANTOL/TOLLAN
      COMMON/CIDIAG/ICID
      COMMON/ROTS/JMAX,KMAX,IDUM,KEL21,KEL
      COMMON/MATRIX/NVALV,NVALR,KSTEP
      COMMON/JKAKC/JTHIS,KA,KC
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
C****************************************************************
100   FORMAT(/,1X,'INITIAL ENERGY = ',D20.12)
105   FORMAT(/,1X,'LANCZOS NOT CONVERGED FOR FUNCTIONS:')
110   FORMAT(/,1X,'LANCZOS FULLY CONVERGED')
115   FORMAT(20I4)
C****************************************************************
C**STORE V ON 52
      REWIND 52
C**STORE Z ON 53
      REWIND 53
C**************************************************************
C**************************************************************
C**INITIAL GUESS !!!
C**************************************************************
C**************************************************************
      DO K=1,NVAL
        DO I=1,ISIZE
          X(I,K)=0
        END DO
        XLAST=1.D+10
        DO I=1,ISIZE
          IF(Y(I).LT.XLAST)THEN
            IF(K.EQ.1)THEN
              XLAST=Y(I)
              NVEC(K)=I
            ELSE
              IGOT=0
              DO J=1,K-1
                IF(I.EQ.NVEC(J))IGOT=1
              END DO
              IF(IGOT.EQ.0)THEN
                XLAST=Y(I)
                NVEC(K)=I
              END IF
            END IF
          END IF
        END DO
C**INITIAL GUESS !!!
        NV=NVEC(K)
        X(NV,K)=1
      END DO
C**************************************************************
C**************************************************************
C**********INITIAL SET-UP.....M(rs) = <psi(r)/H/psi(s)>
C**ASSUME INITIAL FUNCTIONS X0 ORTHONORMAL (VIRTUAL SCF STATES)
C**SET INITIAL FUNCTION V0 = X0 (NORMALISED)
C**************************************************************
C**************************************************************
      DO K=1,NVAL
        DO I=1,ISIZE
          V(I,K)=X(I,K)
        END DO
      END DO
C*******************************DISC
      WRITE(52)V
C*******************************DISC
C**FORM Z0 = M.V0 (UNNORMALISED)
      DO K=1,NVAL
        CALL MX(V(1,K),Z(1,K),W,W,ISIZE)
C**FORM S0 = Z0
        DO I=1,ISIZE
          S(I,K)=Z(I,K)
        END DO
      END DO
C*******************************DISC
      WRITE(53)Z
C*******************************DISC
C**FORM V0.(M.V0) = V0.Z0
      DO KV=1,NVAL
        J0=KV*(KV-1)/2
        DO KZ=1,KV
          CALL DOT(V(1,KV),Z(1,KZ),E,ISIZE)
          XA(J0+KZ)=E
          XACOPY(J0+KZ)=E
        END DO
C**V0.(M.V0) = V0.Z0 = <V0/H/V0> = E-VALUE OF V0
        ENERGY(KV)=E*WAVENM
        ELAST(KV)=ENERGY(KV)
      END DO
C**************************************************************
C**************************************************************
C**FORM NEXT BASIS FUNCTION V1 (UNNORMALISED)
C**************************************************************
C**************************************************************
      DO KV=1,NVAL
        J0=KV*(KV-1)/2
        NV=NVEC(KV)
        V(NV,KV)=0
        DO I=1,ISIZE
          IF(I.NE.NV)THEN
            V(I,KV)=-(S(I,KV)-XA(J0+KV)*X(I,KV))/(Y(I)-XA(J0+KV))
          END IF
        END DO
C**ORTHONORMALISE V1 TO V0
        CALL SCHMDB(V,KV,WRK,NVAL,ISIZE,1)
      END DO
C**SCHMIDT POSITIONS DISC 52 FOR NEW V(K)
C*******************************DISC
      WRITE(52)V
C*******************************DISC
C******************************************************
      KGOT=0
      DO K=1,NVAL
        NVEC(K)=100000
        IGOT=0
        DO I=1,ISIZE
          IF(V(I,K).NE.0)IGOT=1
        END DO
C******************************
        IF(IGOT.EQ.0)THEN
          ELAST(K)=0.D0
          NVEC(K)=2
        END IF
CSC     IF(IGOT.EQ.0)THEN
CSC       ELAST(K)=0.D0
CSC       NVEC(K)=2
CSC     END IF
C******************************
        IF(ELAST(K).NE.0.D0)KGOT=1
      END DO
      IF(KGOT.EQ.0)GO TO 2000
      KSIZE=NVAL
C**************************************************************
C**************************************************************
C**START LANCZOS ITERATIONS
C**************************************************************
C**************************************************************
      DO ILAN=2,NCYCLE
        ILOFF=0
C*******************************DISC
        REWIND 53
C*******************************DISC
        DO I=1,ILAN-1
C*******************************DISC
          READ(53)Z
C*******************************DISC
          KVV=0
          DO KV=1,NVAL
            IF(NVEC(KV).GT.ILAN)THEN
              KVV=KVV+1
              IROFF=KSIZE+KVV
              J0=ILOFF+IROFF*(IROFF-1)/2
              KZZ=0
              DO KZ=1,NVAL
                IF(NVEC(KZ).GT.I)THEN
                  KZZ=KZZ+1
                  CALL DOT(V(1,KV),Z(1,KZ),P,ISIZE)
                  XA(J0+KZZ)=P
                  XACOPY(J0+KZZ)=P
                END IF
              END DO
            END IF
          END DO
          ILOFF=ILOFF+NVAL
          DO K=1,NVAL
            IF(NVEC(K).LE.I)ILOFF=ILOFF-1
          END DO
        END DO
C**FORM Z = M.V
        DO K=1,NVAL
          IF(NVEC(K).GT.ILAN)CALL MX(V(1,K),Z(1,K),W,W,ISIZE)
        END DO
C*******************************DISC
        WRITE(53)Z
C*******************************DISC
        KVV=0
        DO KV=1,NVAL
          IF(NVEC(KV).GT.ILAN)THEN
            KVV=KVV+1
            IROFF=KSIZE+KVV
            J0=KSIZE+IROFF*(IROFF-1)/2
            KZZ=0
            DO KZ=1,KV
              IF(NVEC(KZ).GT.ILAN)THEN
                KZZ=KZZ+1
                CALL DOT(V(1,KV),Z(1,KZ),P,ISIZE)
                XA(J0+KZZ)=P
                XACOPY(J0+KZZ)=P
              END IF
            END DO
          END IF
        END DO
C**GET CURRENT MATRIX SIZE FROM PREVIOUS SIZE
        KSIZE=KSIZE+NVAL
        DO K=1,NVAL
          IF(NVEC(K).LE.ILAN)KSIZE=KSIZE-1
        END DO
C**SET UP SQUARE MATRIX FOR QL
        IF(.NOT.LGIV)THEN
          CALL QLCOPY(XACOPY,KSIZE,XK,NVAL*NCYCLE)
        END IF
C**DIAGONALISE MATRIX ORDER KSIZE (NO OVERLAP)
        CALL DIAG(XA,XK,NVAL*NCYCLE,KSIZE,-1,SUP4,EVAL,WK,
     1  NVAL,NVAL,IP,ISIZMX,NMODE,XK,XK,IASSIG,ISIZMX,J21,NS)
C**************************************************************
C**************************************************************
C**FORM NEW FUNCTION X.....LOWEST-ENERGY EIGENFUNCTIONS K (NORMALISED)
C**************************************************************
C**************************************************************
        DO K=1,NVAL
CNCH      IF(ELAST(K).NE.0.D0)THEN
            DO I=1,ISIZE
              X(I,K)=0
            END DO
CNCH      END IF
        END DO
        JOFF=0
C*******************************DISC
        REWIND 52
C*******************************DISC
        DO J=1,ILAN
C*******************************DISC
          READ(52)V
C*******************************DISC
          DO K=1,NVAL
CNCH        IF(ELAST(K).NE.0.D0)THEN
              KOFF=JOFF
              DO N=1,NVAL
                IF(NVEC(N).GT.J)THEN
                  KOFF=KOFF+1
                  DO I=1,ISIZE
                    X(I,K)=X(I,K)+XK(KOFF,K)*V(I,N)
                  END DO
                END IF
              END DO
CNCH        END IF
          END DO
          JOFF=KOFF
        END DO
        DO K=1,NVAL
          IF(ELAST(K).EQ.0.D0)THEN
            DO I=1,ISIZE
              WRK(I,K)=X(I,K)
              X(I,K)=WRK(I,K)
            END DO
          END IF
        END DO
C**FORM NEW FUNCTION S
        DO K=1,NVAL
CNCH      IF(ELAST(K).NE.0.D0)ENERGY(K)=EVAL(K)*WAVENM
          ENERGY(K)=EVAL(K)*WAVENM
          DO I=1,ISIZE
            S(I,K)=0
          END DO
        END DO
        JOFF=0
C*******************************DISC
        REWIND 53
C*******************************DISC
        DO J=1,ILAN
C*******************************DISC
          READ(53)Z
C*******************************DISC
          DO K=1,NVAL
            KOFF=JOFF
            DO N=1,NVAL
              IF(NVEC(N).GT.J)THEN
                KOFF=KOFF+1
                DO I=1,ISIZE
                  S(I,K)=S(I,K)+XK(KOFF,K)*Z(I,N)
                END DO
              END IF
            END DO
          END DO
          JOFF=KOFF
        END DO
        IF(ILAN.EQ.NCYCLE)GO TO 1000
        KGOT=0
        DO K=1,NVAL
          IF(ELAST(K).NE.0.D0)THEN
            IF(DABS(ENERGY(K)-ELAST(K)).LT.TOLLAN)THEN
              ELAST(K)=0.D0
C**NVEC(K) IS FIRST CYCLE IN WHICH FUNCTION 'K' IS NOT INCLUDED
              NVEC(K)=ILAN+1
            ELSE
              ELAST(K)=ENERGY(K)
            END IF
          END IF
          IF(ELAST(K).NE.0.D0)KGOT=1
        END DO
        IF(KGOT.EQ.0)GO TO 2000
C**************************************************************
C**************************************************************
C**FORM NEXT BASIS FUNCTION V(ILAN+1) (UNNORMALISED)
C**************************************************************
C**************************************************************
        KGOT=0
        DO K=1,NVAL
          IF(ELAST(K).NE.0.D0)THEN
            DO I=1,ISIZE
              V(I,K)=-(S(I,K)-EVAL(K)*X(I,K))/
     1        (Y(I)-EVAL(K))
            END DO
C**ORTHONORMALISE V(ILAN+1) TO V(ILAN)
            CALL SCHMDB(V,K,WRK,NVAL,ISIZE,ILAN)
            IGOT=0
            DO I=1,ISIZE
              IF(V(I,K).NE.0)IGOT=1
            END DO
C******************************************
            IF(IGOT.EQ.0)THEN
              ELAST(K)=0.D0
              NVEC(K)=ILAN+1
            END IF
CSC         IF(IGOT.EQ.0)THEN
CSC           ELAST(K)=0.D0
CSC           NVEC(K)=ILAN+1
CSC         END IF
C******************************************
          ELSE
            DO I=1,ISIZE
              V(I,K)=0
            END DO
          END IF
          IF(ELAST(K).NE.0.D0)KGOT=1
        END DO
        IF(KGOT.EQ.0)GO TO 2000
C**SCHMIDT POSITIONS DISC 52 FOR NEW V(K)
C**OTHERWISE DISC 52 ALREADY POSITIONED FORMING NEW X
C*******************************DISC
        WRITE(52)V
C*******************************DISC
C**RE-LOAD MATRIX ORDER ILAN IN PREPARATION FOR MATRIX ORDER ILAN+1
        DO I=1,KSIZE*(KSIZE+1)/2
          XA(I)=XACOPY(I)
        END DO
      END DO
C**************************************************************
C**************************************************************
C**END LANCZOS ITERATIONS
C**************************************************************
C**************************************************************
1000  CONTINUE
      WRITE(IOUT,105)
      N=0
      DO K=1,NVAL
        IF(ELAST(K).NE.0)THEN
          N=N+1
          NVEC(N)=K
        END IF
      END DO
      WRITE(IOUT,115)(NVEC(K),K=1,N)
      GO TO 3000
2000  CONTINUE
      WRITE(IOUT,110)
3000  CONTINUE
C**************************************************************
C**************************************************************
C**ASSIGN LEVELS
C**************************************************************
C**************************************************************
      MT=1
      IF(EVL.NE.0.0D0)GO TO 4
      MT=2
4     CONTINUE
      IF(EVL.EQ.0)EVL=ENERGY(1)
      DO 7 I=1,NVAL
      ELAST(I)=ENERGY(I)/WAVENM
      ENERGY(I)=ENERGY(I)-EVL
7     CONTINUE
      XMAX(1)=2
      IF(ICID.EQ.0)THEN
        NSIZE=ISIZE
        KKC=KC
        KKA=KA
      ELSE
CC      J21=2*JTHIS+1
        NSIZE=NVALV
      END IF
      ICYCL=1
      IF(JPRINT.LT.-1.AND.ICID.EQ.0)ICYCL=2
      IF(JPRINT.LT.-2.AND.ICID.EQ.0)ICYCL=3
      DO 10 I=1,NVAL
      IF(ENERGY(I).GT.CUT)GO TO 1001
      DO 101 ICYC=1,ICYCL
      NSX=NS
      IOFF=NSX-NVSYM
      IF(IOFF.GT.0)NSX=IOFF
      IF(ICID.NE.0)THEN
        KKC=JTHIS
        KKA=0
        IF(TRIAT.AND.NS.GT.NVSYM)THEN
          KKC=JTHIS-1
          KKA=1
        END IF
      END IF
      JX=I+1-MT
      if (icyc.lt.3) then
c  find two max ci coeffs
        XMAX(ICYC+1)=0.D0
        DO 9 J=1,ISIZE
        IF(X(J,I).EQ.XMAX(ICYC))GO TO 9
        IF(ABS(X(J,I)).LE.ABS(XMAX(ICYC+1)))GO TO 9
        XMAX(ICYC+1)=X(J,I)
        IMAX=J
9       CONTINUE
      else
c  get third biggest coeff
        xmax(4)=0.0
        do 90 j=1,isize
        if((x(j,i).eq.xmax(2)).or.(x(j,i).eq.xmax(3))) go to 90
        if (abs(x(j,i)).lt.abs(xmax(4))) go to 90
        xmax(4)=x(j,i)
        imax=j
90      continue
      endif
      IIA=0
12    CONTINUE
      IF(NSIZE.GE.IMAX)GO TO 13
      IMAX=IMAX-NSIZE
      IIA=IIA+1
      IF(TRIAT)THEN
        KKA=KKA+1
        IF(MOD(IIA,2).EQ.0)KKC=KKC-2
      ELSE
        IA=MOD(IIA,2)
        KKA=KKA+IA
        KKC=JTHIS-KKA
        KKC=KKC+MOD(IIA+1,2)
      END IF
C**UPDATE VIBRATIONAL SYMMETRY FOR NEXT TIME IF REQUIRED
      IF((.NOT.TRIAT).AND.NVSYM.EQ.4)THEN
        IF(MOD(IIA,2).NE.0)THEN
          NSX=NSX+2
        ELSE
          IF(MOD(IIA,4).NE.0)THEN
            NSX=NSX+(-1)**NS
          ELSE
            NSX=NSX-(-1)**NS
          END IF
        END IF
        IF(NSX.GT.NVSYM)THEN
          INCR=NSX-NVSYM
          NSX=INCR
        END IF
        IF(NSX.EQ.0)NSX=NVSYM
      ELSE
        NSX=NSX+1
        IF(NSX.GT.NVSYM)NSX=1
      END IF
      GO TO 12
13    CONTINUE
      IF(ICID.EQ.0)THEN
C**STORE CI ASSIGNMENTS
        KROT=KEL
        IASSIG(I,KROT,ICYC,NSX)=IMAX
        JJMAX(ICYC)=IMAX
      ELSE
C**RECALL CI ASSIGNMENTS
        IOFF=0
        IF(NSX.GT.1)THEN
          DO K=1,NSX-1
            IOFF=IOFF+NTOT(K)
          END DO
        END IF
        IF(TRIAT)THEN
          KROT=2*KKA
          KOFF=1
          IF(NS.GT.NVSYM)KOFF=0
          JROT=KROT+KOFF
        ELSE
          KROT=2*KKA
          IF(KKA.EQ.0.OR.KKA+KKC.NE.JTHIS)KROT=KROT+1
          JROT=KROT
        END IF
        KKEL=(JROT-1)/KSTEP+1
        JJMAX(1)=IOFF+IASSIG(IMAX,KKEL,1,NSX)
        IF(JPRINT.LT.-1)JJMAX(2)=IOFF+IASSIG(IMAX,KKEL,2,NSX)
        IF(JPRINT.LT.-2)JJMAX(3)=IOFF+IASSIG(IMAX,KKEL,3,NSX)
      END IF
101   CONTINUE
      CALL PRCI(ISTAT,NSTAT,NMODE,JJMAX(1),JJMAX(2),JJMAX(3),
     1XMAX(2),XMAX(3),XMAX(4),JX,ENERGY(I),EVL,JTHIS,KKA,KKC,45)
1001  CONTINUE
10    CONTINUE
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE SCHMDB(V,KV,WRK,NVAL,ISIZE,NV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(ISIZE,NVAL),WRK(ISIZE,NVAL)
      COMMON/FILASS/IOUT,INP
C**SCHMIDT ORTHOGONALISE
      N=NV+1
C**NEW FUNCTION IS V(ISIZE,KV,N)
C**NORMALISE NEW FUNCTION
      CALL DOT(V(1,KV),V(1,KV),S,ISIZE)
      SQ=SQRT(S)
C**BEWARE!!! IS 10**-6 CORRECT?
C**WE WILL MAKE THE PROGRAM PRINT...
C**EACH TIME THIS IS OBEYED
      IF(SQ.GT.1.D-6)THEN
        S=1/SQ
      ELSE
        WRITE(IOUT,*)'ZERO VECTOR NO. ',KV
        S=0
C**POSITION DISC 52 FOR V
        REWIND 52
        DO I=1,N
C************************DISC
          IF(I.NE.N)READ(52)WRK
C************************DISC
        END DO
      END IF
      DO I=1,ISIZE
        V(I,KV)=V(I,KV)*S
      END DO
      IF(S.EQ.0)RETURN
C**IF KV=1, OLD FUNCTIONS ARE V(ISIZE,NVAL,NV)
C**IF KV>1, OLD FUNCTIONS ALSO INCLUDE V(ISIZE,KV-1,N)
C************************DISC
      REWIND 52
C************************DISC
      DO I=1,N
C************************DISC
        IF(I.NE.N)READ(52)WRK
C************************DISC
        NXVAL=NVAL
        IF(I.EQ.N)NXVAL=KV
        DO K=1,NXVAL
          IF(I.EQ.N)THEN
            DO J=1,ISIZE
              WRK(J,K)=V(J,K)
            END DO
          END IF
          CALL DOT(V(1,KV),WRK(1,K),S,ISIZE)
          IF(I.EQ.N.AND.K.EQ.KV)THEN
            SQ=SQRT(S)
C**BEWARE!!! IS 10**-6 CORRECT?
C**WE WILL MAKE THE PROGRAM PRINT...
C**EACH TIME THIS IS OBEYED
            IF(SQ.GT.1.D-6)THEN
              S=1.D0-1.D0/SQ
            ELSE
              WRITE(IOUT,*)'REJECTED VECTOR NO. ',K
              CALL FLUSH(IOUT)
              DO J=1,ISIZE
                V(J,KV)=0
                WRK(J,KV)=0
              END DO
            END IF
          END IF
          DO J=1,ISIZE
            V(J,KV)=V(J,KV)-WRK(J,K)*S
          END DO
        END DO
      END DO
C**THE ABOVE IS SCHMIDT OF W.M.
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE MX(X,V,WK,WKR,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 WK(N)
      REAL*4 WKR(N)
      DIMENSION X(N),V(N)
      COMMON/COUPLE/ICOUPL,JCOUPL
C**LOOP ROUND ROWS OF 'M'
      REWIND 20
      ISIZE=N+1
      DO I=1,N
        V(I)=0
      END DO
      DO I=1,N
        ISIZE=ISIZE-1
C**READ Ith ROW FROM DISK INTO WK
        IF(JCOUPL.GT.0)THEN
          READ(20)(WK(J),J=1,ISIZE)
        ELSE
          READ(20)(WKR(J),J=1,ISIZE)
        END IF
C**FORM M.X
        IF(JCOUPL.GT.0)THEN
          DO J=1,ISIZE
            II=N-ISIZE+J
            V(II)=V(II)+WK(J)*X(I)
          END DO
          DO J=1,ISIZE
            II=N-ISIZE+J
            V(I)=V(I)+WK(J)*X(II)
          END DO
          V(I)=V(I)-WK(1)*X(I)
        ELSE
          DO J=1,ISIZE
            II=N-ISIZE+J
            V(II)=V(II)+WKR(J)*X(I)
          END DO
          DO J=1,ISIZE
            II=N-ISIZE+J
            V(I)=V(I)+WKR(J)*X(II)
          END DO
          V(I)=V(I)-WKR(1)*X(I)
        END IF
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DOT(X,V,W,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),V(N)
      W=0
      DO J=1,N
        W=W+X(J)*V(J)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE LANOUT(XK,WK,WKR,Y,ZK,ZKR,N,ISTART,IEND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 WK(N),ZK(N)
      REAL*4 WKR(N),ZKR(N)
      DIMENSION XK(1),Y(N)
      COMMON/COUPLE/ICOUPL,JCOUPL
      IF(ISTART.EQ.1)REWIND 20
      L0=-ISTART+1
      DO I=ISTART,IEND
        J0=L0
        Y(I)=XK(J0+I)
        I0=I-1
C**ROWS TO I-1
        IF(I0.NE.0)THEN
          IF(ISTART.EQ.1)THEN
C**ALREADY IN XK
            KK0=N
            JJ0=I
            DO J=1,I0
              IF(JCOUPL.GT.0)THEN
                WK(J)=XK(JJ0)
              ELSE
                WKR(J)=XK(JJ0)
              END IF
              JJ0=JJ0+KK0-1
              KK0=KK0-1
            END DO
          ELSE
C**GET FROM DISC
C           REWIND 20
C           IF(JCOUPL.GT.0)THEN
C             DO J=1,I0
C               READ(20)ZK
C               WK(J)=ZK(I)
C             END DO
C           ELSE
C             DO J=1,I0
C               READ(20)ZKR
C               WKR(J)=ZKR(I)
C             END DO
C           END IF
          END IF
        END IF
C**REMAINING COLUMNS
        IF(JCOUPL.GT.0)THEN
          K=0
          DO J=I,N
            K=K+1
            WK(K)=XK(J0+J)
          END DO
          WRITE(20)(WK(J),J=1,K)
        ELSE
          K=0
          DO J=I,N
            K=K+1
            WKR(K)=XK(J0+J)
          END DO
          WRITE(20)(WKR(J),J=1,K)
        END IF
        L0=L0+N-I
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE QLCOPY(XACOPY,KSIZE,XK,NXK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XACOPY(1),XK(NXK,1)
      L=0
      DO I=1,KSIZE
        DO J=1,I
          L=L+1
          XK(J,I)=XACOPY(L)
          XK(I,J)=XK(J,I)
        END DO
      END DO
      RETURN
      END
