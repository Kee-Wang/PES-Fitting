C***********************************************************
C***********************************************************
C**MEMO
C***********************************************************
C***********************************************************
      PARAMETER(MAXSIZ=20000000,MSEG=100)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FILASS/IOUT,INP
      COMMON/CWORK/W(MAXSIZ)
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/TRANSF/LTRAN
C**TEMPORARY
C     OPEN(20,FILE='/scratch/STUART/lanczos20',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(21,FILE='/scratch2/STUART/vibints21',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(22,FILE='/scratch2/STUART/vibints22',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(23,FILE='/scratch2/STUART/vibints23',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(24,FILE='/scratch2/STUART/vibints24',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(31,FILE='/ray/stuart/rotints31',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(32,FILE='/ray/stuart/rotints32',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(33,FILE='/ray/stuart/rotints33',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(34,FILE='/ray/stuart/rotints34',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(35,FILE='/ray/stuart/rotints35',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(36,FILE='/ray/stuart/rotints36',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(37,FILE='/ray/stuart/rotints37',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(38,FILE='/ray/stuart/rotints38',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(39,FILE='/ray/stuart/rotints39',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(40,FILE='/ray/stuart/SCFenergy40',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(45,FILE='/ray/stuart/CIenergy45',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(50,FILE='/ray/stuart/scratch50',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(52,FILE='/ray/stuart/lanczos52',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(53,FILE='/ray/stuart/lanczos53',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(54,FILE='/ray/stuart/lanczos54',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(55,FILE='/ray/stuart/lanczos55',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(58,FILE='/ray/stuart/scratch58',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(59,FILE='/ray/stuart/scratch59',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(60,FILE='/ray/stuart/dump60',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(61,FILE='/ray/stuart/rotation61',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(62,FILE='/ray/stuart/rotation62',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(63,FILE='/ray/stuart/rotation63',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(64,FILE='/ray/stuart/rotation64',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(71,FILE='/scratch/STUART/potential71',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(72,FILE='/scratch/STUART/potential72',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(73,FILE='/scratch/STUART/potential73',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(74,FILE='/scratch/STUART/potential74',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(81,FILE='/scratch1/STUART/coriolis81',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(82,FILE='/scratch1/STUART/coriolis82',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(83,FILE='/scratch1/STUART/coriolis83',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(84,FILE='/scratch1/STUART/coriolis84',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(91,FILE='/ray/stuart/coriolis91',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(92,FILE='/ray/stuart/coriolis92',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(93,FILE='/ray/stuart/coriolis93',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C     OPEN(94,FILE='/ray/stuart/coriolis94',FORM='UNFORMATTED',
C    1status='UNKNOWN')
C**TEMPORARY
      LTRAN=MAXSIZ
      INP=1
      IOUT=2
C     CALL TIMIT(2)
      NADD=MSEG
      DO 123 N=1,NADD
      KADD(N)=0
123   LADD(N)=0
      DO 9999 N=1,MAXSIZ
9999  W(N)=0
      NSEG=0
      KFREE=MAXSIZ
      LFREE=1
      KINF=MAXSIZ
      CALL VSCF(W)
      WRITE(IOUT,66) MAXSIZ-KINF,MAXSIZ
 66   FORMAT('USED MEMORY:',I10,3X,'OF',I10)
      CALL TIMIT(4)
      STOP 'END OF MULTIMODE'
      END
C***********************************************************
      SUBROUTINE MEMO(MM,LA1,KA1,LA2,KA2,LA3,KA3,LA4,KA4,LA5,KA5)
      PARAMETER(MSEG=100)
      COMMON/FILASS/IOUT,INP
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
100   FORMAT(1X,'TOO MANY ENTRIES TO MEMO',/)
101   FORMAT(/,1X,'NO. ARRAYS MUST BE NON-ZERO',/)
102   FORMAT(/,1X,'ERROR IN MEMO CALLING SEQUENCE',/)
      M=MM
      IF(M.EQ.0)THEN
      WRITE(IOUT,101)
      WRITE(IOUT,102)
      STOP 'NO. ARRAYS ZERO'
      END IF
      IF(IABS(M).GT.5)THEN
      WRITE(IOUT,100)
      WRITE(IOUT,102)
      STOP 'TOO MANY ARRAYS'
      END IF
      MADD=1
 1    IF(IABS(MM).GE.1) CALL MEM1(M,LA1,KA1)
      MADD=2
 2    IF(IABS(MM).GE.2) CALL MEM1(M,LA2,KA2)
      MADD=3
 3    IF(IABS(MM).GE.3) CALL MEM1(M,LA3,KA3)
      MADD=4
 4    IF(IABS(MM).GE.4) CALL MEM1(M,LA4,KA4)
      MADD=5
 5    IF(IABS(MM).GE.5) CALL MEM1(M,LA5,KA5)
      KINF=MIN0(KINF,KFREE)
      RETURN
      END
C***************************************************
C***************************************************
      SUBROUTINE MEM1(M,LA,KA)
C.....RESERVES OR FREES MEMORY OF LENGTH KA AT ADDRESS LA
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(MSEG=100)
      COMMON/FILASS/IOUT,INP
      COMMON/CMEMO/NADD,NSEG,KFREE,LFREE,KINF,MADD
      COMMON/CSIZE/KADD(MSEG)
      COMMON/CADDR/LADD(MSEG)
      COMMON/CWORK/W(1)
500   FORMAT(/,1X,'ARRAY NUMBER ',I3,' ADDRESS ',I8,' SIZE ',I8,
     1' ALREADY EXISTS',/)
501   FORMAT(/,1X,'ARRAY TO BE ADDED MUST HAVE POSITIVE LENGTH',/)
502   FORMAT(/,1X,'MAXIMUM NO. ARRAYS REACHED.....',
     1'ADJUST MSEG PARAMETER (C.F. ALSO LADD,KADD)',/)
503   FORMAT(/,1X,'MEMORY EXHAUSTED.....',/,
     1         1X,I10,' AVAILABLE',4X,I10,' REQUESTED',4X,
     2         'REQUIRES EXTRA ',I10,/,1X,'ADJUST MAXSIZ PARAMETER',/)

504   FORMAT(/,1X,'ERROR IN MEMO CALLING SEQUENCE',/)
505   FORMAT(/,1X,'ARRAY ADDRESS ',I8,' SIZE ',I8,' NOT FOUND',/)
506   FORMAT(/,1X,'ARRAY TO BE ADDED MUST HAVE ZERO ADDRESS',/)
507   FORMAT(/,1X,'ARRAY TO BE DELETED HAS NON-POSITIVE ADDRESS ',I8,/)
508   FORMAT(/,1X,'ARRAY AT ADDRESS ',I8,' SIZE ',I8,' DOES NOT EXIST',
     1/)
509   FORMAT(/,1X,'ARRAY NUMBER ',I3,' ADDRESS ',I8,' LENGTH ',I8,
     1' HAS WRONG SIZE ',I8,/)
510   FORMAT(/,1X,'ARRAY ADDRESS NAME NOT IN COMMON/CADDR/',/)
511   FORMAT(/,1X,'ERROR IN MEMO SEQUENCE NUMBER ',I2,/)
      IF(M.LT.0) GOTO 100
C....................................CHECK FOR AVAILABLE MEMORY
      IF(LA.NE.0)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,506)
      NUM=0
      DO 1 I=1,NADD
      IF(LADD(I).NE.LA)GO TO 1
      NUM=I
1     CONTINUE
      IF(NUM.EQ.0)THEN
      WRITE(IOUT,505)LA,KA
      WRITE(IOUT,504)
      ELSE
      WRITE(IOUT,500)NUM,LA,KA
      END IF
      STOP 'ARRAY ADDRESS NON-ZERO'
      END IF
C***********
      IF(KA.LE.0)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,501)
      STOP 'ARRAY LENGTH NON-POSITIVE'
      END IF
C***********
      IF(KA.GT.KFREE)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,503)KFREE,KA,KA-KFREE
      STOP 'MEMORY EXHAUSTED'
      END IF
C***********
      IF(NSEG.EQ.MSEG)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,502)
      STOP 'TOO MANY ARRAYS'
      END IF
C***********
      NSEG=NSEG+1
      LA=LFREE
      NUM=0
      DO 2 I=1,NADD
      IF(LADD(I).NE.LA)GO TO 2
      NUM=I
2     CONTINUE
      IF(NUM.EQ.0)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,510)
      WRITE(IOUT,504)
      STOP 'ARRAY NOT FOUND'
      ELSE
      KADD(NUM)=KA
      END IF
      KFREE=KFREE-KA
      LFREE=LFREE+KA
      RETURN
C....................................FIND SEGMENT TO BE FREED
100   CONTINUE
      IF(LA.LE.0)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,507)LA
      STOP 'ARRAY ADDRESS NOT POSITIVE'
      END IF
C***********
      NUM=0
      DO 3 I=1,NADD
      IF(LADD(I).NE.LA)GO TO 3
      NUM=I
3     CONTINUE
      IF(NUM.EQ.0)THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,508)LA,KA
      WRITE(IOUT,510)
      WRITE(IOUT,504)
      STOP 'ARRAY NOT FOUND'
      ELSE
C***********
      IF(KA.NE.KADD(NUM))THEN
      WRITE(IOUT,511)MADD
      WRITE(IOUT,509)NUM,LA,KADD(NUM),KA
      WRITE(IOUT,504)
      STOP 'ARRAY HAS WRONG SIZE'
      END IF
      END IF
C***********
      KTOT=0
      DO 4 I=1,NADD
      IF(LADD(I).GT.LA)KTOT=KTOT+KADD(I)
4     CONTINUE
      IF(KTOT.GT.0)THEN
      DO 5 I=1,KTOT
      W(LA-1+I)=W(LA-1+I+KA)
5     CONTINUE
      END IF
      DO 6 I=1,NADD
      IF(LADD(I).GT.LA)LADD(I)=LADD(I)-KA
6     CONTINUE
      LA=0
      LFREE=LFREE-KA
      KFREE=KFREE+KA
      NSEG=NSEG-1
      RETURN
      END
