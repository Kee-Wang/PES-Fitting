C**************************************************************
C**************************************************************
      SUBROUTINE USERIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FILASS/IOUT,INP
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/DUMP/JDUMP(10),IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(148),NTERM,
     1KTERM(148,6),ITERM(148,6),JTERM(148,6),NDUMMY,
     2CDIP(100,3),NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3)
C     COMMON/ALPHA/ALPHA1,ALPHA2,ALPHA3
      COMMON/AXES/MX(3),MXDIP(3)
C**RVIB3 POTPAR
      COMMON/MORSE/ALPHA1(3),ALPHA2(3),ALX1(5),ALX2(5)
      COMMON/CUBIC/BLX1(5),BLX2(5),BLX3(5)
      COMMON/EXPAND/IEXPND(5)
      COMMON/VPOT1/C1(400,5),INDEX1(400,3,5),RO1(3,5),IGR1(3,5),
     1ISY1(3,5),IC1(3,5),NPA1(5),RP11(3,5),RP21(3,5),RP31(3,5),
     2ICTYP1(3,5),MOLSM1(5),CPL1(3,3,5)
C**RVIB3 POTPAR
200   FORMAT(//,1X,'EQUILIBRIUM GEOMETRY',/)
205   FORMAT(1X,'R1 = ',F10.6,'  R2 = ',F10.6,'  R3 = ',F10.6,/)
210   FORMAT(1X,'TH1 = ',F10.6,'  TH2 = ',F10.6,'  PHI = ',F10.6,/)
215   FORMAT(/,1X,'POTENTIAL TERMS: ITERM, JTERM, CTERM',/)
220   FORMAT(1X,6I2,1X,6I2,1X,D20.12)
225   FORMAT(/,1X,'PA X-COMPONENT',1X,'SYMMETRY ',I1,/)
230   FORMAT(/,1X,'PA Y-COMPONENT',1X,'SYMMETRY ',I1,/)
235   FORMAT(/,1X,'PA Z-COMPONENT',1X,'SYMMETRY ',I1,/)
240   FORMAT(/,1X,'DIPOLE TERMS: ITERM, JTERM, CTERM',/)
C**RVIB3 POTPAR
100   FORMAT(15I4)
101   FORMAT(3F15.10)
 102   FORMAT(3(3I2,F15.10))
C102   FORMAT(3(3I2,E15.8))
500   FORMAT(/,1X,I4,' COEFFICIENTS',/)
501   FORMAT(/,1X,' REFERENCE GEOMETRY:   ',3F15.10)
502   FORMAT(/,1X,' MORSE ALPHAS FOR BONDS 1 AND 2: ',2F15.10)
503   FORMAT(/,1X,' CUBIC BEND COEFFICIENTS: ',3F15.10)
504   FORMAT(/,1X,22(1H*),/)
505   FORMAT(/,1X,'FORBIDDEN POTENTIAL - ICTYP1(3) = ',I3)
506   FORMAT(/,1X,'POTENTIAL TERMS',/)
C***********************************************************
      READ(INP,*)
      IX=1
      READ(INP,100)NPA1(IX),IGR1(1,IX),IGR1(2,IX),IGR1(3,IX),ISY1(1,IX),
     1ISY1(2,IX),ISY1(3,IX),IC1(1,IX),IC1(2,IX),IC1(3,IX),
     2ICTYP1(1,IX),ICTYP1(2,IX),ICTYP1(3,IX),MOLSM1(IX)
      IF(ICTYP1(3,IX).EQ.4.OR.ICTYP1(3,IX).EQ.5)THEN
        WRITE(IOUT,505)ICTYP1(3,IX)
        STOP 'FORBIDDEN POTENTIAL'
      END IF
      IEXPND(IX)=ICTYP1(1,IX)+ICTYP1(3,IX)
      ICTYP1(1,IX)=IABS(ICTYP1(1,IX))
      ICTYP1(2,IX)=IABS(ICTYP1(2,IX))
      ICTYP1(3,IX)=IABS(ICTYP1(3,IX))
      DO I=1,NPA1(IX)
      READ(INP,*)INDEX1(I,1,IX),INDEX1(I,2,IX),INDEX1(I,3,IX),
     1C1(I,IX)
      END DO
      READ(INP,101) RO1(1,IX),RO1(2,IX),RO1(3,IX)
C     RO1(1,IX)=RO1(1,IX)
C     RO1(2,IX)=RO1(2,IX)
C     RO1(3,IX)=RO1(3,IX)
      RO1(3,IX)=RO1(3,IX)/RAD
      READ(INP,101) RP11(1,IX),RP11(2,IX),RP11(3,IX)
      READ(INP,101) RP21(1,IX),RP21(2,IX),RP21(3,IX)
      READ(INP,101) RP31(1,IX),RP31(2,IX),RP31(3,IX)
      READ(INP,101) CPL1(1,1,IX),CPL1(1,2,IX),CPL1(1,3,IX)
      READ(INP,101) CPL1(2,1,IX),CPL1(2,2,IX),CPL1(2,3,IX)
      READ(INP,101) CPL1(3,1,IX),CPL1(3,2,IX),CPL1(3,3,IX)
C***********************************************SURFIT-TYPE INPUT BOX
      WRITE(IOUT,506)
      WRITE(IOUT,500) NPA1(IX)
      WRITE(IOUT,102)(INDEX1(I,1,IX),INDEX1(I,2,IX),INDEX1(I,3,IX),
     1C1(I,IX),I=1,NPA1(IX))
      WRITE(IOUT,501) RO1(1,IX),RO1(2,IX),RO1(3,IX)*RAD
      IF(ICTYP1(1,IX).EQ.2.OR.ICTYP1(2,IX).EQ.2)THEN
        WRITE(IOUT,502)RP11(1,IX),RP11(2,IX)
      END IF
      IF(ICTYP1(3,IX).EQ.3)
     1WRITE(IOUT,503)RP11(3,IX),RP21(3,IX),RP31(3,IX)
      WRITE(IOUT,504)
      IF(IABS(IEXPND(IX)).EQ.2.OR.IABS(IEXPND(IX)).EQ.5.OR.
     1IABS(IEXPND(IX)).EQ.8.OR.IABS(IEXPND(IX)).EQ.11)THEN
        ALX1(IX)=RP11(1,IX)
        ALX2(IX)=RP11(2,IX)
      END IF
      IF(IABS(IEXPND(IX)).GT.2.AND.IABS(IEXPND(IX)).LT.6)THEN
        BLX1(IX)=RP11(3,IX)
        BLX2(IX)=RP21(3,IX)
        BLX3(IX)=RP31(3,IX)
      END IF
C**RVIB3 POTPAR
      IF(LDUMP.EQ.0)RETURN
C*************************
      DO L=2,3
        READ(INP,*)
      READ(INP,100)MXDIP(L),NPA1(L),IGR1(1,L),IGR1(2,L),IGR1(3,L),
     1ISY1(1,L),ISY1(2,L),ISY1(3,L),IC1(1,L),IC1(2,L),IC1(3,L),
     2ICTYP1(1,L),ICTYP1(2,L),ICTYP1(3,L),MOLSM1(L)
      IF(ICTYP1(3,L).EQ.4.OR.ICTYP1(3,L).EQ.5)THEN
        WRITE(IOUT,505)ICTYP1(3,L)
        STOP 'FORBIDDEN POTENTIAL'
      END IF
      IEXPND(L)=ICTYP1(1,L)+ICTYP1(3,L)
      ICTYP1(1,L)=IABS(ICTYP1(1,L))
      ICTYP1(2,L)=IABS(ICTYP1(2,L))
      ICTYP1(3,L)=IABS(ICTYP1(3,L))
      READ(INP,102)(INDEX1(I,1,L),INDEX1(I,2,L),INDEX1(I,3,L),
     1C1(I,L),I=1,NPA1(L))
      READ(INP,101) RO1(1,L),RO1(2,L),RO1(3,L)
      READ(INP,101) RP11(1,L),RP11(2,L),RP11(3,L)
      READ(INP,101) RP21(1,L),RP21(2,L),RP21(3,L)
      READ(INP,101) RP31(1,L),RP31(2,L),RP31(3,L)
      READ(INP,101) CPL1(1,1,L),CPL1(1,2,L),CPL1(1,3,L)
      READ(INP,101) CPL1(2,1,L),CPL1(2,2,L),CPL1(2,3,L)
      READ(INP,101) CPL1(3,1,L),CPL1(3,2,L),CPL1(3,3,L)
C***********************************************SURFIT-TYPE INPUT BOX
        IF(L.EQ.1)WRITE(IOUT,225)MXDIP(1)
        IF(L.EQ.2)WRITE(IOUT,230)MXDIP(2)
        IF(L.EQ.3)WRITE(IOUT,235)MXDIP(3)
      WRITE(IOUT,500) NPA1(L)
      WRITE(IOUT,102)(INDEX1(I,1,L),INDEX1(I,2,L),INDEX1(I,3,L),
     1C1(I,L),I=1,NPA1(L))
      WRITE(IOUT,501) RO1(1,L),RO1(2,L),RO1(3,L)*RAD
      IF(ICTYP1(1,L).EQ.2.OR.ICTYP1(2,L).EQ.2)THEN
        WRITE(IOUT,502)RP11(1,L),RP11(2,L)
      END IF
      IF(ICTYP1(3,L).EQ.3)
     1WRITE(IOUT,503)RP11(3,L),RP21(3,L),RP31(3,L)
      WRITE(IOUT,504)
      IF(IABS(IEXPND(L)).EQ.2.OR.IABS(IEXPND(L)).EQ.5.OR.
     1IABS(IEXPND(L)).EQ.8)THEN
        ALX1(L)=RP11(1,L)
        ALX2(L)=RP11(2,L)
      END IF
      IF(IABS(IEXPND(L)).GT.2.AND.IABS(IEXPND(L)).LT.6)THEN
        BLX1(L)=RP11(3,L)
        BLX2(L)=RP21(3,L)
        BLX3(L)=RP31(3,L)
      END IF

      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPOT(V,NATOM,XX,RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION QQ(6)
C     COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(148),NTERM,
C    1KTERM(148,6),ITERM(148,6),JTERM(148,6),NDUMMY,
C    2CDIP(100,3),NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3)
C     COMMON/ALPHA/ALPHA1,ALPHA2,ALPHA3
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/ECKCNT/ICNT,INTC
      COMMON/FILASS/IOUT,INP
C**RVIB3 POTPAR
      COMMON/MORSE/ALPHA1(3),ALPHA2(3),ALX1(5),ALX2(5)
      COMMON/CUBIC/BLX1(5),BLX2(5),BLX3(5)
      COMMON/EXPAND/IEXPND(5)
      COMMON/VPOT1/C1(400,5),INDEX1(400,3,5),RO1(3,5),IGR1(3,5),
     1ISY1(3,5),IC1(3,5),NPA1(5),RP11(3,5),RP21(3,5),RP31(3,5),
     2ICTYP1(3,5),MOLSM1(5),CPL1(3,3,5)
C**RVIB3 POTPAR
200   FORMAT(1X,6F12.6)
C**RVIB3 POTPAR
500   FORMAT(1X,E14.7,3F13.7)
C**RVIB3 POTPAR
      CALL BONDS(NATOM,RR,XX)
      R1=RR(1,3)
      R2=RR(1,2)
      R3=RR(2,3)
      CAL=(R1*R1+R2*R2-R3*R3)/(2*R1*R2)
      IF(CAL.GE.1)CAL=1
      IF(CAL.LE.-1)CAL=-1
      S1=R1
      S2=R2
      AL=DACOS(CAL)
C**RVIB3 POTPAR
      NPA=NPA1(1)
      IND=1
      V=FCT(S1,S2,AL,NPA,C1(1,IND),INDEX1(1,1,IND),
     1RO1(1,IND),RO1(2,IND),RO1(3,IND),MOLSM1(IND),CPL1(1,1,IND),IND)
      IF(AL*RAD.LT.110)V1=1.0
C**RVIB3 POTPAR
      RETURN
      END
C**RVIB3 POTPAR
C***********************************************************
C***********************************************************
      FUNCTION FCT(R1,R2,R3,NPA,C,INDEX,R10,R20,R30,MOLSYM,CPL,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/EXPAND/IEXPND(5)
      COMMON/FILASS/IOUT,INP,IP
      DIMENSION C(400),INDEX(400,3),IGR(3),ISY(3),IC(3)
      DIMENSION Q(3),ICTYP(3),RP1(3),RP2(3),RP3(3),CPL(3,3)
      CALL QCOORD(Q,R1,R2,R3,R10,R20,R30,MOLSYM,CPL,IND)
      FCT=0.D0
      DO I=1,NPA
        F=1.D0
        DO J=1,3
          IF (INDEX(I,J).NE.0) THEN
            F=F*Q(J)**INDEX(I,J)
          END IF
        END DO
        FCT=FCT+C(I)*F
      END DO
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE QCOORD(Q,R1,R2,R3,R10,R20,R30,MOLSYM,CPL,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RC(3),ROC(3),RO(3),CPL(3,3),Q(3)
      COMMON/FILASS/IOUT,INP,IP
      COMMON/MORSE/ALPHA1(3),ALPHA2(3),A1(5),A2(5)
      COMMON/CUBIC/B1(5),B2(5),B3(5)
      COMMON/EXPAND/IEXPN(5)
      IF(MOLSYM.GT.0)THEN
        RO(1)=R10
        RO(2)=R20
        RO(3)=R30
        RC(1)=R1
        RC(2)=R2
        RC(3)=R3
        CALL COUPLL(RC,RC,CPL)
        CALL COUPLL(RO,ROC,CPL)
        Q(1)=RC(1)-ROC(1)
        Q(2)=RC(2)-ROC(2)
        Q(3)=RC(3)-ROC(3)
      ELSE
        Q(1)=R1-R10
        Q(2)=R2-R20
        Q(3)=R3-R30
      END IF
      IF(IEXPN(IND).EQ.0)RETURN
      IF(IABS(IEXPN(IND)).EQ.1.OR.IABS(IEXPN(IND)).EQ.4.OR.
     1IABS(IEXPN(IND)).EQ.7)THEN
        Q(1)=Q(1)/R1
        Q(2)=Q(2)/R2
      END IF
      IF(IABS(IEXPN(IND)).EQ.2.OR.IABS(IEXPN(IND)).EQ.5.OR.
     1IABS(IEXPN(IND)).EQ.8.OR.IABS(IEXPN(IND)).EQ.11)THEN
        Q(1)=(1.D0-EXP(-A1(IND)*Q(1)))
        Q(2)=(1.D0-EXP(-A2(IND)*Q(2)))
      END IF
      IF (IABS(IEXPN(IND)).GT.2.AND.IABS(IEXPN(IND)).LT.6)THEN
        Q(3)=Q(3)*(B1(IND)+Q(3)*(B2(IND)+Q(3)*B3(IND)))
      END IF
      IF (IABS(IEXPN(IND)).GT.8)THEN
C       IF(Q(3).GT.RAD)Q(3)=RAD
C       IF(Q(3).LT.0)Q(3)=0
        Q(3)=COS(Q(3))
      END IF
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE COUPLL(R,S,CPL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(3),S(3),T(3),CPL(3,3)
      COMMON/FILASS/IOUT,INP,IP
      DO I=1,3
        T(I)=0.D0
        DO J=1,3
          T(I)=T(I)+CPL(I,J)*R(J)
        END DO
      END DO
      DO I=1,3
        S(I)=T(I)
      END DO
      RETURN
      END
C***********************************************************
C***********************************************************
C**RVIB3 POTPAR
C*****************************
C*****************************
      SUBROUTINE GETQPT
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETAPT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDIP(V,NATOM,XX,RR,LDIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM)
      DIMENSION QQ(6)
C     COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(148),NTERM,
C    1KTERM(148,6),ITERM(148,6),JTERM(148,6),NDUMMY,
C    2CDIP(100,3),NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3)
C     COMMON/ALPHA/ALPHA1,ALPHA2,ALPHA3
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/ECKCNT/ICNT,INTC
      COMMON/FILASS/IOUT,INP
C**RVIB3 POTPAR
      COMMON/MORSE/ALPHA1(3),ALPHA2(3),ALX1(5),ALX2(5)
      COMMON/CUBIC/BLX1(5),BLX2(5),BLX3(5)
      COMMON/EXPAND/IEXPND(5)
      COMMON/VPOT1/C1(400,5),INDEX1(400,3,5),RO1(3,5),IGR1(3,5),
     1ISY1(3,5),IC1(3,5),NPA1(5),RP11(3,5),RP21(3,5),RP31(3,5),
     2ICTYP1(3,5),MOLSM1(5),CPL1(3,3,5)
C**RVIB3 POTPAR
200   FORMAT(1X,6F12.6)
C**RVIB3 POTPAR
500   FORMAT(1X,E14.7,3F13.7)
C**RVIB3 POTPAR
      CALL BONDS(NATOM,RR,XX)
      R1=RR(1,3)
      R2=RR(1,2)
      R3=RR(2,3)
      CAL=(R1*R1+R2*R2-R3*R3)/(2*R1*R2)
      IF(CAL.GE.1)CAL=1
      IF(CAL.LE.-1)CAL=-1
      S1=R1
      S2=R2
      AL=DACOS(CAL)
C**RVIB3 POTPAR
      NPA=NPA1(LDIP)
      IND=LDIP
      V=FCT(S1,S2,AL,NPA,C1(1,IND),INDEX1(1,1,IND),
     1RO1(1,IND),RO1(2,IND),RO1(3,IND),MOLSM1(IND),CPL1(1,1,IND),IND)
C**RVIB3 POTPAR
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETQDT
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE MINPOT
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETPRP(VDP,NATOM,XX,RR,QQ,NMODE,IDIP,TAU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM),QQ(NMODE)
      CALL BONDS(NATOM,RR,XX)
      R1=RR(1,3)
      R2=RR(1,4)
      R3=RR(1,2)
      R4=RR(2,3)
      R5=RR(2,4)
      R6=RR(3,4)
      TH1=DACOS((R1*R1+R3*R3-R4*R4)/(2*R1*R3))
      TH2=DACOS((R2*R2+R3*R3-R5*R5)/(2*R2*R3))
      ALPHA=DACOS((R1*R1+R2*R2-R6*R6)/(2*R1*R2))
      COSPHI=((COS(ALPHA)-COS(TH1)*COS(TH2))/(SIN(TH1)*SIN(TH2)))
      IF(COSPHI.GE.1)COSPHI=1
      IF(COSPHI.LE.-1)COSPHI=-1
      PHI=DACOS(COSPHI)
      IF(IDIP.EQ.1)VDP=R1
      IF(IDIP.EQ.2)VDP=COS(ALPHA)
      RETURN
      END
C*****************************
      SUBROUTINE GETPOL
       RETURN
      END
C*****************************

