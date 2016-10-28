C****************************************************************
C****************************************************************
C**ROT
C****************************************************************
C****************************************************************
      SUBROUTINE ROTC(XX,XM,NATOM,A,B,C,BBAR,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),XM(NATOM),WR(3),E(3)
      COMMON/MOMI/XK(3,3),XMU(3,3)
      IF(IND.NE.0)THEN
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
      END IF
C**DIAGONALISE MOMENT OF INERTIA MATRIX
      CALL DIAG(XK,XK,3,3,-1,WR,E,WR,3,3,XK,3,3,XK,XK,XK,IDUM,IDUM,
     1IDUM)
C**E ARE IN ASCENDING ORDER
      A=1/(2*E(1))
      B=1/(2*E(2))
      C=1/(2*E(3))
      BBAR=(B+C)/2
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI0(NATOM,X0,XM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X0(NATOM,3),XM(NATOM)
      COMMON/FILASS/IOUT,INP
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
100   FORMAT(//,1X,'ROTATIONAL CONSTANTS AT EQUILIBRIUM ',3F10.5,///)
      CALL ROTC(X0,XM,NATOM,AROT,BROT,CROT,BBAR,1)
      WRITE(IOUT,100)AROT*WAVENM,BROT*WAVENM,CROT*WAVENM
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI1(XQ,MM,NMODE,NATOM,QQ,XZ,AB,B,AA,BB,
     1XX,X0,XL,XM,MODE,V,VR,S,J21,V0,W0,E0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM)
      REAL*4 VR(J21,MM)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ(MM),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/AXES/MXYZ(3)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
      MX=MXYZ(1)
      MY=MXYZ(2)
      MZ=MXYZ(3)
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE)=1
      DO M=1,MM
        QQ(MODE)=XQ(M)
        CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
        IF(KMAX.GT.0.AND.ICI.LT.0)THEN
          IF(JCOUPL.GT.0)THEN
            DO J=1,J21
              V(J,M)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
     1        S(J,3,J)*XMU(MZ,MZ))/2
            END DO
          ELSE
            DO J=1,J21
              VR(J,M)=(S(J,1,J)*XMU(MX,MX)+S(J,2,J)*XMU(MY,MY)+
     1        S(J,3,J)*XMU(MZ,MZ))/2
            END DO
          END IF
        ELSE
          CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
          DO I=1,J21
            DO J=1,J21
              V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
            END DO
          END DO
          CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,V0,
     1    V0,IDUM,IDUM,IDUM)
          IF(JCOUPL.GT.0)THEN
            DO J=1,J21
              V(J,M)=E0(J)
            END DO
          ELSE
            DO J=1,J21
              VR(J,M)=E0(J)
            END DO
          END IF
        END IF
      END DO
      IF(JCOUPL.GT.0)THEN
        WRITE(61)V
      ELSE
        WRITE(61)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI2(XQ1,XQ2,MM1,MM2,NMODE,NATOM,QQ,XZ,AB,B,AA,BB,
     1XX,X0,XL,XM,MODE1,MODE2,V,VR,S,J21,V0,W0,E0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM2,MM1)
      REAL*4 VR(J21,MM2,MM1)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/AXES/MXYZ(3)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
      MX=MXYZ(1)
      MY=MXYZ(2)
      MZ=MXYZ(3)
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE1)=1
      MULT(MODE2)=1
      DO M1=1,MM1
        QQ(MODE1)=XQ1(M1)
        DO M2=1,MM2
          QQ(MODE2)=XQ2(M2)
          CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
          IF(KMAX.GT.0.AND.ICI.LT.0)THEN
            IF(JCOUPL.GT.0)THEN
              DO J=1,J21
                V(J,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1          S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
              END DO
            ELSE
              DO J=1,J21
                VR(J,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1          S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
              END DO
            END IF
          ELSE
            CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
            DO I=1,J21
              DO J=1,J21
                V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
              END DO
            END DO
            CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,V0,
     1      V0,V0,IDUM,IDUM,IDUM)
            IF(JCOUPL.GT.0)THEN
              DO J=1,J21
                V(J,M2,M1)=E0(J)
              END DO
            ELSE
              DO J=1,J21
                VR(J,M2,M1)=E0(J)
              END DO
            END IF
          END IF
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        WRITE(62)V
      ELSE
        WRITE(62)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI3(XQ1,XQ2,XQ3,MM1,MM2,MM3,NMODE,NATOM,QQ,
     1XZ,AB,B,AA,BB,XX,X0,XL,XM,MODE1,MODE2,MODE3,V,VR,S,J21,V0,W0,E0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM3,MM2,MM1)
      REAL*4 VR(J21,MM3,MM2,MM1)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/FILASS/IOUT,INP
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/AXES/MXYZ(3)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
      MX=MXYZ(1)
      MY=MXYZ(2)
      MZ=MXYZ(3)
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE1)=1
      MULT(MODE2)=1
      MULT(MODE3)=1
      DO M1=1,MM1
        QQ(MODE1)=XQ1(M1)
        DO M2=1,MM2
          QQ(MODE2)=XQ2(M2)
          DO M3=1,MM3
            QQ(MODE3)=XQ3(M3)
            CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
            IF(KMAX.GT.0.AND.ICI.LT.0)THEN
              IF(JCOUPL.GT.0)THEN
                DO J=1,J21
                  V(J,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1            S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
                END DO
              ELSE
                DO J=1,J21
                  VR(J,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1            S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
                END DO
              END IF
            ELSE
              CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
              DO I=1,J21
                DO J=1,J21
                  V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
                END DO
              END DO
              CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,
     1        V0,V0,V0,IDUM,IDUM,IDUM)
              IF(JCOUPL.GT.0)THEN
                DO J=1,J21
                  V(J,M3,M2,M1)=E0(J)
                END DO
              ELSE
                DO J=1,J21
                  VR(J,M3,M2,M1)=E0(J)
                END DO
              END IF
            END IF
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        WRITE(63)V
      ELSE
        WRITE(63)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETMI4(XQ1,XQ2,XQ3,XQ4,MM1,MM2,MM3,MM4,NMODE,NATOM,QQ,
     1XZ,AB,B,AA,BB,XX,X0,XL,XM,MODE1,MODE2,MODE3,MODE4,V,VR,S,J21,V0,
     2W0,E0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 V(J21,MM4,MM3,MM2,MM1)
      REAL*4 VR(J21,MM4,MM3,MM2,MM1)
      DIMENSION V0(J21,J21),W0(J21),E0(J21)
      DIMENSION XZ(NMODE,NMODE,3),AB(NMODE,3),B(NMODE,NMODE)
      DIMENSION AA(NMODE,3,3),BB(NMODE)
      DIMENSION S(J21,9,J21)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4),QQ(NMODE)
      DIMENSION XX(NATOM,3),X0(NATOM,3),XL(NATOM,NMODE,3),XM(NATOM)
      COMMON/ROTS/JMAX,KMAX
      COMMON/CIDIAG/ICID,ICI,JCI
      COMMON/AXES/MXYZ(3)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/MOMI/XK(3,3),XMU(3,3)
C**TEMPORARY (DIMENSIONS)
      COMMON/MULT/MULT(1000)
      MX=MXYZ(1)
      MY=MXYZ(2)
      MZ=MXYZ(3)
      DO K=1,NMODE
        QQ(K)=0
        MULT(K)=0
      END DO
      MULT(MODE1)=1
      MULT(MODE2)=1
      MULT(MODE3)=1
      MULT(MODE4)=1
      DO M1=1,MM1
        QQ(MODE1)=XQ1(M1)
        DO M2=1,MM2
          QQ(MODE2)=XQ2(M2)
          DO M3=1,MM3
            QQ(MODE3)=XQ3(M3)
            DO M4=1,MM4
              QQ(MODE4)=XQ4(M4)
              CALL CORIOL(NMODE,NATOM,QQ,XZ,AB,B,AA,BB,XX,X0,XL,XM,ZZ)
              IF(KMAX.GT.0.AND.ICI.LT.0)THEN
                IF(JCOUPL.GT.0)THEN
                  DO J=1,J21
                    V(J,M4,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1              S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
                  END DO
                ELSE
                  DO J=1,J21
                    VR(J,M4,M3,M2,M1)=(S(J,1,J)*XMU(MX,MX)+
     1              S(J,2,J)*XMU(MY,MY)+S(J,3,J)*XMU(MZ,MZ))/2
                  END DO
                END IF
              ELSE
                CALL ROTC(XX,XM,NATOM,AROT,BROT,CROT,BBAR,0)
                DO I=1,J21
                  DO J=1,J21
                    V0(J,I)=S(J,3,I)*AROT+S(J,1,I)*BROT+S(J,2,I)*CROT
                  END DO
                END DO
                CALL DIAG(V0,V0,J21,J21,-1,W0,E0,W0,J21,J21,V0,J21,J21,
     1          V0,V0,V0,IDUM,IDUM,IDUM)
                IF(JCOUPL.GT.0)THEN
                  DO J=1,J21
                    V(J,M4,M3,M2,M1)=E0(J)
                  END DO
                ELSE
                  DO J=1,J21
                    VR(J,M4,M3,M2,M1)=E0(J)
                  END DO
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
      IF(JCOUPL.GT.0)THEN
        WRITE(64)V
      ELSE
        WRITE(64)VR
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI1(NMODE,MODE,H,XQ,XA,NN,MM,IP,ISIZE,
     2VM,VMR,J21,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM,6)
      REAL*4 VMR(MM,6)
      DIMENSION H(NN,MM,3),XQ(MM),XA(1)
      DIMENSION IP(ISIZE,NMODE)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      FACT=0
      IF(ICOUPL.EQ.1)FACT=1/FACTOR
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.1)IFACT=1
        IF(ICOUPL.EQ.2)IFACT=-(NMODE-2)
        IF(ICOUPL.EQ.3)THEN
          IFACT=(NMODE-3)*(NMODE-1)
          DO I=2,NMODE-2
            IFACT=IFACT-I
          END DO
        END IF
        IF(ICOUPL.EQ.4)THEN
          IFACT=-(NMODE-4)*(NMODE-3)*(NMODE-1)
          DO I=1,NMODE-4
            IFACT=IFACT-I*(NMODE-3-I)
          END DO
          DO I=1,NMODE-4
            IFACT=IFACT+(NMODE-2)*I
          END DO
          DO I=2,NMODE-2
            IFACT=IFACT+(NMODE-4)*I
          END DO
        END IF
      END IF
      FACTRC=FACT
      IF(JCOUPL.GT.0)THEN
        READ(91)VM
      ELSE
        READ(91)VMR
      END IF
      DO M=1,MM
        IF(JCOUPL.GT.0)THEN
          TERM=VM(M,IABC)*FACTRC
        ELSE
          TERM=VMR(M,IABC)*FACTRC
        END IF
C**ISIZE IS NO. UNIQUE INTEGRALS (1-DIM)
        DO IRHS=1,ISIZE
          NR=IP(IRHS,MODE)
          X=TERM*H(NR,M,1)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL=IP(ILHS,MODE)
            Y=H(NL,M,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*X
          END DO
        END DO
      END DO
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI1(NMODE,MODE,H,XQ,XA,XK,NN,MM,IPL,IPR,ISIZMX,
     1ISIZEL,ISIZER,IP1,ISIZE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H(NN,MM,3),XQ(MM),XA(ISIZEL,ISIZER),XK(1)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP1(ISIZE1,1)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM1A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI1'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR=IPR(IRHS,MODE)
C**FIND RHS INDEX (TRIVIAL CASE)
        DO IR=1,ISIZE1
          IF(NR.EQ.IP1(IR,1))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL=IPL(ILHS,MODE)
C**FIND LHS INDEX (TRIVIAL CASE)
          DO IL=1,ISIZE1
            IF(NL.EQ.IP1(IL,1))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM1A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM1A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI2(NMODE,MODE1,MODE2,H1,XQ1,H2,XQ2,NN1,MM1,NN2,MM2,
     1XA,IP,ISIZE,TEMP,JCI,VM,VMR,J21,IABC,IP1,ISIZE1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM2,MM1,12)
      REAL*4 VMR(MM2,MM1,12)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),XA(1)
      DIMENSION IP(ISIZE,NMODE),IP1(ISIZE1,1)
      DIMENSION TEMP(JCI,JCI,3)
      DIMENSION XQ1(MM1),XQ2(MM2)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      FACT=0
      IF(ICOUPL.EQ.2)FACT=1/FACTOR
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.2)IFACT=1
        IF(ICOUPL.EQ.3)IFACT=-(NMODE-3)
        IF(ICOUPL.EQ.4)THEN
          IFACT=(NMODE-4)*(NMODE-3)
          DO I=1,NMODE-4
            IFACT=IFACT-I
          END DO
        END IF
      END IF
      FACTRC=FACT
      IF(JCOUPL.GT.0)THEN
        READ(92)VM
      ELSE
        READ(92)VMR
      END IF
      IOFF=(IABC-7)*2
      DO M1=1,MM1
        DO K=1,2
          DO ILHS2=1,JCI
            DO IRHS2=1,JCI
              TEMP(ILHS2,IRHS2,K)=0
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
          X2=0
          IF(JCOUPL.GT.0)THEN
            DO M2=1,MM2
              DO IX1=1,ISIZE1
                IRHS2=IP1(IX1,1)
                YY1=H2(IRHS2,M2,1)*FACTRC
                X1=VM(M2,M1,IABC)*YY1
                DO IY1=1,ISIZE1
                  ILHS2=IP1(IY1,1)
                  Y0=H2(ILHS2,M2,1)
                  TEMP(ILHS2,IRHS2,1)=TEMP(ILHS2,IRHS2,1)+Y0*X1
                  TEMP(ILHS2,IRHS2,2)=TEMP(ILHS2,IRHS2,2)+Y0*X2
                END DO
              END DO
            END DO
          ELSE
            DO M2=1,MM2
              DO IX1=1,ISIZE1
                IRHS2=IP1(IX1,1)
                YY1=H2(IRHS2,M2,1)*FACTRC
                X1=VMR(M2,M1,IABC)*YY1
                DO IY1=1,ISIZE1
                  ILHS2=IP1(IY1,1)
                  Y0=H2(ILHS2,M2,1)
                  TEMP(ILHS2,IRHS2,1)=TEMP(ILHS2,IRHS2,1)+Y0*X1
                  TEMP(ILHS2,IRHS2,2)=TEMP(ILHS2,IRHS2,2)+Y0*X2
                END DO
              END DO
            END DO
          END IF
        ELSE
          IF(JCOUPL.GT.0)THEN
            DO M2=1,MM2
              DO IX1=1,ISIZE1
                IRHS2=IP1(IX1,1)
                Y1=H2(IRHS2,M2,1)*IFACT
                Y2=H2(IRHS2,M2,2)*IFACT
                X1=VM(M2,M1,8+IOFF)*Y2
                X2=VM(M2,M1,7+IOFF)*Y1
                DO IY1=1,ISIZE1
                  ILHS2=IP1(IY1,1)
                  Y0=H2(ILHS2,M2,1)
                  TEMP(ILHS2,IRHS2,1)=TEMP(ILHS2,IRHS2,1)+Y0*X1
                  TEMP(ILHS2,IRHS2,2)=TEMP(ILHS2,IRHS2,2)+Y0*X2
                END DO
              END DO
            END DO
          ELSE
            DO M2=1,MM2
              DO IX1=1,ISIZE1
                IRHS2=IP1(IX1,1)
                Y1=H2(IRHS2,M2,1)*IFACT
                Y2=H2(IRHS2,M2,2)*IFACT
                X1=VMR(M2,M1,8+IOFF)*Y2
                X2=VMR(M2,M1,7+IOFF)*Y1
                DO IY1=1,ISIZE1
                  ILHS2=IP1(IY1,1)
                  Y0=H2(ILHS2,M2,1)
                  TEMP(ILHS2,IRHS2,1)=TEMP(ILHS2,IRHS2,1)+Y0*X1
                  TEMP(ILHS2,IRHS2,2)=TEMP(ILHS2,IRHS2,2)+Y0*X2
                END DO
              END DO
            END DO
          END IF
        END IF
C**ISIZE IS NO. UNIQUE INTEGRALS (2-DIM)
        DO IRHS=1,ISIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          X1=H1(NR1,M1,1)
          X2=H1(NR1,M1,2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            Y=H1(NL1,M1,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*(TEMP(NL2,NR2,1)*X1+
     1      TEMP(NL2,NR2,2)*X2)
          END DO
        END DO
      END DO
      IF(ITIM2A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM2A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI2(NMODE,MODE1,MODE2,H1,XQ1,H2,XQ2,NN1,MM1,NN2,MM2,
     1XA,XK,IPL,IPR,ISIZMX,ISIZEL,ISIZER,IP2,ISIZE2,JCI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3)
      DIMENSION XA(ISIZEL,ISIZER),XK(1),
     1IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP2(ISIZE2,2)
      DIMENSION XQ1(MM1),XQ2(MM2)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM2A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI2'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
C**FIND RHS INDEX
        DO IR=1,ISIZE2
          IF(NR1.EQ.IP2(IR,1).AND.NR2.EQ.IP2(IR,2))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.(IPR(IRHS,K).NE.
     1      IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
C**FIND LHS INDEX
          DO IL=1,ISIZE2
            IF(NL1.EQ.IP2(IL,1).AND.NL2.EQ.IP2(IL,2))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM2A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM2A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI3(NMODE,MODE1,MODE2,MODE3,H1,XQ1,H2,XQ2,
     1H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,XA,IP,ISIZE,TEMP,JCI,VM,
     3VMR,J21,IABC,IP2,ISIZE2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM3,MM2,MM1,15)
      REAL*4 VMR(MM3,MM2,MM1,15)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION IP(ISIZE,NMODE),IP2(ISIZE2,2)
      DIMENSION XA(1),TEMP(JCI,JCI,JCI,JCI,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      FACT=0
      IF(ICOUPL.EQ.3)FACT=1/FACTOR
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.3)IFACT=1
        IF(ICOUPL.EQ.4)IFACT=-(NMODE-4)
      END IF
      FACTRC=FACT
      IF(JCOUPL.GT.0)THEN
        READ(93)VM
      ELSE
        READ(93)VMR
      END IF
      IOFF=(IABC-7)*3
      DO M1=1,MM1
        DO K=1,2
          DO ILHS2=1,JCI
            DO ILHS3=1,JCI
              DO IRHS2=1,JCI
                DO IRHS3=1,JCI
                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,K)=0
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
          X2=0
          IF(JCOUPL.GT.0)THEN
            DO IX2=1,ISIZE2
              IRHS2=IP2(IX2,1)
              IRHS3=IP2(IX2,2)
              DO M2=1,MM2
                YY1=H2(IRHS2,M2,1)*FACTRC
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  X1=VM(M3,M2,M1,IABC)*YY1*Z1
                  DO IY2=1,ISIZE2
                    ILHS2=IP2(IY2,1)
                    ILHS3=IP2(IY2,2)
                    Y0=H2(ILHS2,M2,1)
                    Z0=H3(ILHS3,M3,1)
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                  END DO
                END DO
              END DO
            END DO
          ELSE
            DO IX2=1,ISIZE2
              IRHS2=IP2(IX2,1)
              IRHS3=IP2(IX2,2)
              DO M2=1,MM2
                YY1=H2(IRHS2,M2,1)*FACTRC
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  X1=VMR(M3,M2,M1,IABC)*YY1*Z1
                  DO IY2=1,ISIZE2
                    ILHS2=IP2(IY2,1)
                    ILHS3=IP2(IY2,2)
                    Y0=H2(ILHS2,M2,1)
                    Z0=H3(ILHS3,M3,1)
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                  END DO
                END DO
              END DO
            END DO
          END IF
        ELSE
          IF(JCOUPL.GT.0)THEN
            DO IX2=1,ISIZE2
              IRHS2=IP2(IX2,1)
              IRHS3=IP2(IX2,2)
              DO M2=1,MM2
                Y1=H2(IRHS2,M2,1)*IFACT
                Y2=H2(IRHS2,M2,2)*IFACT
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  Z2=H3(IRHS3,M3,2)
                  X1=VM(M3,M2,M1,8+IOFF)*Y2*Z1+VM(M3,M2,M1,9+IOFF)*
     1            Y1*Z2
                  X2=VM(M3,M2,M1,7+IOFF)*Y1*Z1
                  DO IY2=1,ISIZE2
                    ILHS2=IP2(IY2,1)
                    ILHS3=IP2(IY2,2)
                    Y0=H2(ILHS2,M2,1)
                    Z0=H3(ILHS3,M3,1)
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                  END DO
                END DO
              END DO
            END DO
          ELSE
            DO IX2=1,ISIZE2
              IRHS2=IP2(IX2,1)
              IRHS3=IP2(IX2,2)
              DO M2=1,MM2
                Y1=H2(IRHS2,M2,1)*IFACT
                Y2=H2(IRHS2,M2,2)*IFACT
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  Z2=H3(IRHS3,M3,2)
                  X1=VMR(M3,M2,M1,8+IOFF)*Y2*Z1+VMR(M3,M2,M1,9+IOFF)*
     1            Y1*Z2
                  X2=VMR(M3,M2,M1,7+IOFF)*Y1*Z1
                  DO IY2=1,ISIZE2
                    ILHS2=IP2(IY2,1)
                    ILHS3=IP2(IY2,2)
                    Y0=H2(ILHS2,M2,1)
                    Z0=H3(ILHS3,M3,1)
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                    TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1              TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                  END DO
                END DO
              END DO
            END DO
          END IF
        END IF
C**ISIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,ISIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          NR3=IP(IRHS,MODE3)
          X1=H1(NR1,M1,1)
          X2=H1(NR1,M1,2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            NL3=IP(ILHS,MODE3)
            Y=H1(NL1,M1,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*
     1      (TEMP(NL2,NL3,NR2,NR3,1)*X1+TEMP(NL2,NL3,NR2,NR3,2)*X2)
          END DO
        END DO
      END DO
      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI3(NMODE,MODE1,MODE2,MODE3,H1,XQ1,H2,XQ2,
     1H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,XA,XK,
     2IPL,IPR,ISIZMX,ISIZEL,ISIZER,IP3,ISIZE3,JCI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP3(ISIZE3,3)
      DIMENSION XA(ISIZEL,ISIZER),XK(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
        NR3=IPR(IRHS,MODE3)
C**FIND RHS INDEX
        DO IR=1,ISIZE3
          IF(NR1.EQ.IP3(IR,1).AND.NR2.EQ.IP3(IR,2).AND.
     1       NR3.EQ.IP3(IR,3))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.(
     1      IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
          NL3=IPL(ILHS,MODE3)
C**FIND LHS INDEX
          DO IL=1,ISIZE3
            IF(NL1.EQ.IP3(IL,1).AND.NL2.EQ.IP3(IL,2).AND.
     1         NL3.EQ.IP3(IL,3))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
          XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VXMI3(NMODE,MODE1,MODE2,MODE3,H1,XQ1,H2,XQ2,
     1H3,XQ3,NN1,MM1,NN2,MM2,NN3,MM3,XA,
     2IPL,IPR,ISIZMX,ISIZEL,ISIZER,TEMP,JCI,VM,VMR,J21,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM3,MM2,MM1,15)
      REAL*4 VMR(MM3,MM2,MM1,15)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE)
      DIMENSION XA(ISIZEL,ISIZER),TEMP(JCI,JCI,JCI,JCI,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM3A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VXMI3'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      FACT=0
      IF(ICOUPL.EQ.3)FACT=1/FACTOR
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.3)IFACT=1
        IF(ICOUPL.EQ.4)IFACT=-(NMODE-4)
      END IF
      FACTRC=FACT
      IF(JCOUPL.GT.0)THEN
        READ(93)VM
      ELSE
        READ(93)VMR
      END IF
      IOFF=(IABC-7)*3
      DO M1=1,MM1
        DO K=1,2
          DO ILHS2=1,JCI
            DO ILHS3=1,JCI
              DO IRHS2=1,JCI
                DO IRHS3=1,JCI
                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,K)=0
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
          X2=0
          IF(JCOUPL.GT.0)THEN
            DO M2=1,MM2
C             DO M3=1,MM3
                DO IRHS2=1,JCI
                  YY1=H2(IRHS2,M2,1)*FACTRC
              DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
                    X1=VM(M3,M2,M1,IABC)*YY1*Z1
                    DO ILHS2=1,JCI
                      Y0=H2(ILHS2,M2,1)
                      DO ILHS3=1,JCI
                        Z0=H3(ILHS3,M3,1)
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                      END DO
                    END DO
                  END DO
              END DO
                END DO
C             END DO
            END DO
          ELSE
            DO M2=1,MM2
C             DO M3=1,MM3
                DO IRHS2=1,JCI
                  YY1=H2(IRHS2,M2,1)*FACTRC
              DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
                    X1=VMR(M3,M2,M1,IABC)*YY1*Z1
                    DO ILHS2=1,JCI
                      Y0=H2(ILHS2,M2,1)
                      DO ILHS3=1,JCI
                        Z0=H3(ILHS3,M3,1)
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                      END DO
                    END DO
                  END DO
              END DO
                END DO
C             END DO
            END DO
          END IF
        ELSE
          IF(JCOUPL.GT.0)THEN
            DO M2=1,MM2
C             DO M3=1,MM3
                DO IRHS2=1,JCI
                  Y1=H2(IRHS2,M2,1)*IFACT
                  Y2=H2(IRHS2,M2,2)*IFACT
              DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
                    Z2=H3(IRHS3,M3,2)
                    X1=VM(M3,M2,M1,8+IOFF)*Y2*Z1+VM(M3,M2,M1,9+IOFF)*
     1              Y1*Z2
                    X2=VM(M3,M2,M1,7+IOFF)*Y1*Z1
                    DO ILHS2=1,JCI
                      Y0=H2(ILHS2,M2,1)
                      DO ILHS3=1,JCI
                        Z0=H3(ILHS3,M3,1)
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                      END DO
                    END DO
                  END DO
              END DO
                END DO
C             END DO
            END DO
          ELSE
            DO M2=1,MM2
C             DO M3=1,MM3
                DO IRHS2=1,JCI
                  Y1=H2(IRHS2,M2,1)*IFACT
                  Y2=H2(IRHS2,M2,2)*IFACT
              DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
                    Z2=H3(IRHS3,M3,2)
                    X1=VMR(M3,M2,M1,8+IOFF)*Y2*Z1+VMR(M3,M2,M1,9+IOFF)*
     1              Y1*Z2
                    X2=VMR(M3,M2,M1,7+IOFF)*Y1*Z1
                    DO ILHS2=1,JCI
                      Y0=H2(ILHS2,M2,1)
                      DO ILHS3=1,JCI
                        Z0=H3(ILHS3,M3,1)
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,1)+Y0*Z0*X1
                        TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)=
     1                  TEMP(ILHS2,ILHS3,IRHS2,IRHS3,2)+Y0*Z0*X2
                      END DO
                    END DO
                  END DO
              END DO
                END DO
C             END DO
            END DO
          END IF
        END IF
C**ISIZE IS NO. UNIQUE INTEGRALS (3-DIM)
        DO IRHS=1,ISIZER
          NR1=IPR(IRHS,MODE1)
          NR2=IPR(IRHS,MODE2)
          NR3=IPR(IRHS,MODE3)
          X1=H1(NR1,M1,1)
          X2=H1(NR1,M1,2)
          DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NMODE
              IF(IS.EQ.0)GO TO 1000
              IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.(
     1        IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
            END DO
1000        CONTINUE
            NL1=IPL(ILHS,MODE1)
            NL2=IPL(ILHS,MODE2)
            NL3=IPL(ILHS,MODE3)
            Y=H1(NL1,M1,1)
            XA(ILHS,IRHS)=XA(ILHS,IRHS)+Y*IS*
     1      (TEMP(NL2,NL3,NR2,NR3,1)*X1+TEMP(NL2,NL3,NR2,NR3,2)*X2)
          END DO
        END DO
      END DO
      IF(ITIM3A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM3A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE V0MI4(NMODE,MODE1,MODE2,MODE3,MODE4,H1,XQ1,H2,XQ2,
     1H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,XA,
     2IP,ISIZE,TEMP,JCI,VM,VMR,J21,IABC,IP3,ISIZE3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM4,MM3,MM2,MM1,18)
      REAL*4 VMR(MM4,MM3,MM2,MM1,18)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION IP(ISIZE,NMODE),IP3(ISIZE3,3)
      DIMENSION XA(1),TEMP(JCI,JCI,JCI,JCI,JCI,JCI,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating V0MI4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      FACT=0
      IF(ICOUPL.EQ.4)FACT=1/FACTOR
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.4)IFACT=1
      END IF
      FACTRC=FACT
      IF(JCOUPL.GT.0)THEN
        READ(94)VM
      ELSE
        READ(94)VMR
      END IF
      IOFF=(IABC-7)*4
      DO M1=1,MM1
        DO K=1,2
          DO ILHS2=1,JCI
            DO ILHS3=1,JCI
              DO ILHS4=1,JCI
                DO IRHS2=1,JCI
                  DO IRHS3=1,JCI
                    DO IRHS4=1,JCI
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,K)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
          X2=0
          IF(JCOUPL.GT.0)THEN
            DO IX3=1,ISIZE3
              IRHS2=IP3(IX3,1)
              IRHS3=IP3(IX3,2)
              IRHS4=IP3(IX3,3)
              DO M2=1,MM2
                YY1=H2(IRHS2,M2,1)*FACTRC
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  DO M4=1,MM4
                    W1=H4(IRHS4,M4,1)
                    X1=VM(M4,M3,M2,M1,IABC)*YY1*Z1*W1
                    DO IY3=1,ISIZE3
                      ILHS2=IP3(IY3,1)
                      ILHS3=IP3(IY3,2)
                      ILHS4=IP3(IY3,3)
                      Y0=H2(ILHS2,M2,1)
                      Z0=H3(ILHS3,M3,1)
                      W0=H4(ILHS4,M4,1)
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                +Y0*Z0*W0*X1
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                +Y0*Z0*W0*X2
                    END DO
                  END DO
                END DO
              END DO
            END DO
          ELSE
            DO IX3=1,ISIZE3
              IRHS2=IP3(IX3,1)
              IRHS3=IP3(IX3,2)
              IRHS4=IP3(IX3,3)
              DO M2=1,MM2
                YY1=H2(IRHS2,M2,1)*FACTRC
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  DO M4=1,MM4
                    W1=H4(IRHS4,M4,1)
                    X1=VMR(M4,M3,M2,M1,IABC)*YY1*Z1*W1
                    DO IY3=1,ISIZE3
                      ILHS2=IP3(IY3,1)
                      ILHS3=IP3(IY3,2)
                      ILHS4=IP3(IY3,3)
                      Y0=H2(ILHS2,M2,1)
                      Z0=H3(ILHS3,M3,1)
                      W0=H4(ILHS4,M4,1)
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                +Y0*Z0*W0*X1
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                +Y0*Z0*W0*X2
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END IF
        ELSE
          IF(JCOUPL.GT.0)THEN
            DO IX3=1,ISIZE3
              IRHS2=IP3(IX3,1)
              IRHS3=IP3(IX3,2)
              IRHS4=IP3(IX3,3)
              DO M2=1,MM2
                Y1=H2(IRHS2,M2,1)*IFACT
                Y2=H2(IRHS2,M2,2)*IFACT
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  Z2=H3(IRHS3,M3,2)
                  DO M4=1,MM4
                    W1=H4(IRHS4,M4,1)
                    W2=H4(IRHS4,M4,2)
                    X1=VM(M4,M3,M2,M1,8+IOFF)*Y2*Z1*W1+
     1              VM(M4,M3,M2,M1,9+IOFF)*Y1*Z2*W1+
     2              VM(M4,M3,M2,M1,10+IOFF)*Y1*Z1*W2
                    X2=VM(M4,M3,M2,M1,7+IOFF)*Y1*Z1*W1
                    DO IY3=1,ISIZE3
                      ILHS2=IP3(IY3,1)
                      ILHS3=IP3(IY3,2)
                      ILHS4=IP3(IY3,3)
                      Y0=H2(ILHS2,M2,1)
                      Z0=H3(ILHS3,M3,1)
                      W0=H4(ILHS4,M4,1)
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                +Y0*Z0*W0*X1
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                +Y0*Z0*W0*X2
                    END DO
                  END DO
                END DO
              END DO
            END DO
          ELSE
            DO IX3=1,ISIZE3
              IRHS2=IP3(IX3,1)
              IRHS3=IP3(IX3,2)
              IRHS4=IP3(IX3,3)
              DO M2=1,MM2
                Y1=H2(IRHS2,M2,1)*IFACT
                Y2=H2(IRHS2,M2,2)*IFACT
                DO M3=1,MM3
                  Z1=H3(IRHS3,M3,1)
                  Z2=H3(IRHS3,M3,2)
                  DO M4=1,MM4
                    W1=H4(IRHS4,M4,1)
                    W2=H4(IRHS4,M4,2)
                    X1=VMR(M4,M3,M2,M1,8+IOFF)*Y2*Z1*W1+
     1              VMR(M4,M3,M2,M1,9+IOFF)*Y1*Z2*W1+
     2              VMR(M4,M3,M2,M1,10+IOFF)*Y1*Z1*W2
                    X2=VMR(M4,M3,M2,M1,7+IOFF)*Y1*Z1*W1
                    DO IY3=1,ISIZE3
                      ILHS2=IP3(IY3,1)
                      ILHS3=IP3(IY3,2)
                      ILHS4=IP3(IY3,3)
                      Y0=H2(ILHS2,M2,1)
                      Z0=H3(ILHS3,M3,1)
                      W0=H4(ILHS4,M4,1)
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                +Y0*Z0*W0*X1
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                +Y0*Z0*W0*X2
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END IF
        END IF
C**ISIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,ISIZE
          NR1=IP(IRHS,MODE1)
          NR2=IP(IRHS,MODE2)
          NR3=IP(IRHS,MODE3)
          NR4=IP(IRHS,MODE4)
          X1=H1(NR1,M1,1)
          X2=H1(NR1,M1,2)
          J0=IRHS*(IRHS-1)/2
          DO ILHS=1,IRHS
            NL1=IP(ILHS,MODE1)
            NL2=IP(ILHS,MODE2)
            NL3=IP(ILHS,MODE3)
            NL4=IP(ILHS,MODE4)
            Y=H1(NL1,M1,1)
            XA(ILHS+J0)=XA(ILHS+J0)+Y*(
     1      TEMP(NL2,NL3,NL4,NR2,NR3,NR4,1)*X1+
     2      TEMP(NL2,NL3,NL4,NR2,NR3,NR4,2)*X2)
          END DO
        END DO
      END DO
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VMI4(NMODE,MODE1,MODE2,MODE3,MODE4,H1,XQ1,H2,XQ2,
     1H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,XA,XK,
     2IPL,IPR,ISIZMX,ISIZEL,ISIZER,IP4,ISIZE4,JCI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE),IP4(ISIZE4,4)
      DIMENSION XA(ISIZEL,ISIZER),XK(1)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VMI4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      DO IRHS=1,ISIZER
        NR1=IPR(IRHS,MODE1)
        NR2=IPR(IRHS,MODE2)
        NR3=IPR(IRHS,MODE3)
        NR4=IPR(IRHS,MODE4)
C**FIND RHS INDEX
        DO IR=1,ISIZE4
          IF(NR1.EQ.IP4(IR,1).AND.NR2.EQ.IP4(IR,2).AND.
     1       NR3.EQ.IP4(IR,3).AND.NR4.EQ.IP4(IR,4))GO TO 1000
        END DO
1000    CONTINUE
        DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
          IS=1
          DO K=1,NMODE
            IF(IS.EQ.0)GO TO 3000
            IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1      MODE4.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
          END DO
C**OVERLAP OF REMAINING STATES
          NL1=IPL(ILHS,MODE1)
          NL2=IPL(ILHS,MODE2)
          NL3=IPL(ILHS,MODE3)
          NL4=IPL(ILHS,MODE4)
C**FIND LHS INDEX
          DO IL=1,ISIZE4
            IF(NL1.EQ.IP4(IL,1).AND.NL2.EQ.IP4(IL,2).AND.
     1         NL3.EQ.IP4(IL,3).AND.NL4.EQ.IP4(IL,4))GO TO 2000
          END DO
2000      CONTINUE
C**GET MATRIX ELEMENT
          MR=IR
          ML=IL
          IF(IR.LT.IL)THEN
            MR=IL
            ML=IR
          END IF
          I=MR*(MR-1)/2+ML
          X=XK(I)
          IF(IL.GT.IR)X=X*IHERM
            XA(ILHS,IRHS)=XA(ILHS,IRHS)+X*IS
3000      CONTINUE
        END DO
      END DO
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE VXMI4(NMODE,MODE1,MODE2,MODE3,MODE4,H1,XQ1,H2,XQ2,
     1H3,XQ3,H4,XQ4,NN1,MM1,NN2,MM2,NN3,MM3,NN4,MM4,XA,
     2IPL,IPR,ISIZMX,ISIZEL,ISIZER,TEMP,JCI,VM,VMR,J21,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 VM(MM4,MM3,MM2,MM1,18)
      REAL*4 VMR(MM4,MM3,MM2,MM1,18)
      DIMENSION H1(NN1,MM1,3),H2(NN2,MM2,3),H3(NN3,MM3,3),H4(NN4,MM4,3)
      DIMENSION IPL(ISIZMX,NMODE),IPR(ISIZMX,NMODE)
      DIMENSION XA(ISIZEL,ISIZER),TEMP(JCI,JCI,JCI,JCI,JCI,JCI,3)
      DIMENSION XQ1(MM1),XQ2(MM2),XQ3(MM3),XQ4(MM4)
      COMMON/HERM/IHERM
      COMMON/COUPLE/ICOUPL,JCOUPL
      COMMON/WHICH/IWHICH
      COMMON/FACTOR/FACTOR
      COMMON/TIMER/ITIM1A,ITIM1B,ITIM2A,ITIM2B,ITIM3A,ITIM3B,
     1ITIM4A,ITIM4B,ITIM
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      IF(ITIM4A.EQ.0)THEN
        WRITE(IOUT,*)'Calculating VXMI4'
        CALL FLUSH(IOUT)
        CALL TIMIT(1)
      END IF
      FACT=0
      IF(ICOUPL.EQ.4)FACT=1/FACTOR
      IFACT=1
      IF(IWHICH.NE.0)THEN
        IF(ICOUPL.EQ.4)IFACT=1
      END IF
      FACTRC=FACT
      IF(JCOUPL.GT.0)THEN
        READ(94)VM
      ELSE
        READ(94)VMR
      END IF
      IOFF=(IABC-7)*4
      DO M1=1,MM1
        DO K=1,2
          DO ILHS2=1,JCI
            DO ILHS3=1,JCI
              DO ILHS4=1,JCI
                DO IRHS2=1,JCI
                  DO IRHS3=1,JCI
                    DO IRHS4=1,JCI
                      TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,K)=0
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
        IF(IABC.LT.7)THEN
          X2=0
          IF(JCOUPL.GT.0)THEN
          DO M2=1,MM2
C           DO M3=1,MM3
C             DO M4=1,MM4
                DO IRHS2=1,JCI
                  YY1=H2(IRHS2,M2,1)*FACTRC
            DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
              DO M4=1,MM4
                    DO IRHS4=1,JCI
                      W1=H4(IRHS4,M4,1)
                      X1=VM(M4,M3,M2,M1,IABC)*YY1*Z1*W1
                      DO ILHS2=1,JCI
                        Y0=H2(ILHS2,M2,1)
                        DO ILHS3=1,JCI
                          Z0=H3(ILHS3,M3,1)
                          DO ILHS4=1,JCI
                            W0=H4(ILHS4,M4,1)
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                     +Y0*Z0*W0*X1
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                     +Y0*Z0*W0*X2
                          END DO
                        END DO
                      END DO
                    END DO
              END DO
                  END DO
            END DO
                END DO
C             END DO
C           END DO
          END DO
          ELSE
          DO M2=1,MM2
C           DO M3=1,MM3
C             DO M4=1,MM4
                DO IRHS2=1,JCI
                  YY1=H2(IRHS2,M2,1)*FACTRC
            DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
              DO M4=1,MM4
                    DO IRHS4=1,JCI
                      W1=H4(IRHS4,M4,1)
                      X1=VMR(M4,M3,M2,M1,IABC)*YY1*Z1*W1
                      DO ILHS2=1,JCI
                        Y0=H2(ILHS2,M2,1)
                        DO ILHS3=1,JCI
                          Z0=H3(ILHS3,M3,1)
                          DO ILHS4=1,JCI
                            W0=H4(ILHS4,M4,1)
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                     +Y0*Z0*W0*X1
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                     +Y0*Z0*W0*X2
                          END DO
                        END DO
                      END DO
                    END DO
              END DO
                  END DO
            END DO
                END DO
C             END DO
C           END DO
          END DO
          END IF
        ELSE
          IF(JCOUPL.GT.0)THEN
          DO M2=1,MM2
C           DO M3=1,MM3
C             DO M4=1,MM4
                DO IRHS2=1,JCI
                  Y1=H2(IRHS2,M2,1)*IFACT
                  Y2=H2(IRHS2,M2,2)*IFACT
            DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
                    Z2=H3(IRHS3,M3,2)
              DO M4=1,MM4
                    DO IRHS4=1,JCI
                      W1=H4(IRHS4,M4,1)
                      W2=H4(IRHS4,M4,2)
                      X1=VM(M4,M3,M2,M1,8+IOFF)*Y2*Z1*W1+
     1                VM(M4,M3,M2,M1,9+IOFF)*Y1*Z2*W1+
     2                VM(M4,M3,M2,M1,10+IOFF)*Y1*Z1*W2
                      X2=VM(M4,M3,M2,M1,7+IOFF)*Y1*Z1*W1
                      DO ILHS2=1,JCI
                        Y0=H2(ILHS2,M2,1)
                        DO ILHS3=1,JCI
                          Z0=H3(ILHS3,M3,1)
                          DO ILHS4=1,JCI
                            W0=H4(ILHS4,M4,1)
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                     +Y0*Z0*W0*X1
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                     +Y0*Z0*W0*X2
                          END DO
                        END DO
                      END DO
                    END DO
              END DO
                  END DO
            END DO
                END DO
C             END DO
C           END DO
          END DO
          ELSE
          DO M2=1,MM2
C           DO M3=1,MM3
C             DO M4=1,MM4
                DO IRHS2=1,JCI
                  Y1=H2(IRHS2,M2,1)*IFACT
                  Y2=H2(IRHS2,M2,2)*IFACT
            DO M3=1,MM3
                  DO IRHS3=1,JCI
                    Z1=H3(IRHS3,M3,1)
                    Z2=H3(IRHS3,M3,2)
              DO M4=1,MM4
                    DO IRHS4=1,JCI
                      W1=H4(IRHS4,M4,1)
                      W2=H4(IRHS4,M4,2)
                      X1=VMR(M4,M3,M2,M1,8+IOFF)*Y2*Z1*W1+
     1                VMR(M4,M3,M2,M1,9+IOFF)*Y1*Z2*W1+
     2                VMR(M4,M3,M2,M1,10+IOFF)*Y1*Z1*W2
                      X2=VMR(M4,M3,M2,M1,7+IOFF)*Y1*Z1*W1
                      DO ILHS2=1,JCI
                        Y0=H2(ILHS2,M2,1)
                        DO ILHS3=1,JCI
                          Z0=H3(ILHS3,M3,1)
                          DO ILHS4=1,JCI
                            W0=H4(ILHS4,M4,1)
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,1)
     2                     +Y0*Z0*W0*X1
                            TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     1                     =TEMP(ILHS2,ILHS3,ILHS4,IRHS2,IRHS3,IRHS4,2)
     2                     +Y0*Z0*W0*X2
                          END DO
                        END DO
                      END DO
                    END DO
              END DO
                  END DO
            END DO
                END DO
C             END DO
C           END DO
          END DO
          END IF
        END IF
C**ISIZE IS NO. UNIQUE INTEGRALS (4-DIM)
        DO IRHS=1,ISIZER
          NR1=IPR(IRHS,MODE1)
          NR2=IPR(IRHS,MODE2)
          NR3=IPR(IRHS,MODE3)
          NR4=IPR(IRHS,MODE4)
          X1=H1(NR1,M1,1)
          X2=H1(NR1,M1,2)
          DO ILHS=1,ISIZEL
C**OVERLAP OF REMAINING STATES
            IS=1
            DO K=1,NMODE
              IF(IS.EQ.0)GO TO 1000
              IF(K.NE.MODE1.AND.K.NE.MODE2.AND.K.NE.MODE3.AND.K.NE.
     1        MODE4.AND.(IPR(IRHS,K).NE.IPL(ILHS,K)))IS=0
            END DO
1000        CONTINUE
            NL1=IPL(ILHS,MODE1)
            NL2=IPL(ILHS,MODE2)
            NL3=IPL(ILHS,MODE3)
            NL4=IPL(ILHS,MODE4)
            Y=H1(NL1,M1,1)
            XA(ILHS,IRHS)=XA(ILHS,IRHS)+IS*Y*(
     1      TEMP(NL2,NL3,NL4,NR2,NR3,NR4,1)*X1+
     2      TEMP(NL2,NL3,NL4,NR2,NR3,NR4,2)*X2)
          END DO
        END DO
      END DO
      IF(ITIM4A.EQ.0)THEN
        CALL TIMIT(3)
        ITIM4A=1
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DUMVM(XK,ISIZEL,ISIZER,IND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XK(ISIZEL,ISIZER)
      COMMON/FILASS/IOUT,INP
      DO IX=1,ISIZER
        WRITE(IND)(XK(IY,IX),IY=1,ISIZEL)
      END DO
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DIAGEL(ISIZE,KROT,XA,EVCI,NVAL,J21,JROT,NS,ISTART,
     1KSTART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ
      DIMENSION XA(1),EVCI(NVAL,J21,1)
      COMMON/LANCZO/LANCZ
      COMMON/MATRIX/IDUM,NVALR,KSTEP,KSIGN
      COMMON/FILASS/IOUT,INP
      IEL=(JROT-1)/KSTEP+1
      IF(LANCZ)THEN
        KSIZE=(KROT-KSTART)*NVAL
        JSIZE=ISIZE-KSIZE-(KSTART-1)*NVAL
        L0=KSIZE*JSIZE+KSIZE*(KSIZE+1)/2
        DO ILHS=1,NVAL
          XA(ILHS+L0)=EVCI(ILHS,IEL,NS)
          L0=L0+JSIZE-ILHS
        END DO
      ELSE
        IOFF=(KROT-1)*NVAL
        DO I=1,NVAL
          JOFF=IOFF+I-1
          J0=JOFF*(JOFF+1)/2+IOFF
          XA(J0+I)=EVCI(I,IEL,NS)
        END DO
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE OFFDEL(ISIZE,XA,CFSL,CFSR,ISIZMX,TEMP,ISIZEL,
     1ISIZER,NVAL,YK,S,J21,KROTL,KROTR,JROTL,JROTR,IELX,IERX,
     2NSL,NSR,ISTART,KSTART,IABC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LANCZ
      DIMENSION XA(1)
      DIMENSION CFSL(ISIZMX,NVAL,1),CFSR(ISIZMX,NVAL,1)
      DIMENSION TEMP(ISIZMX,NVAL)
      DIMENSION YK(ISIZMX,ISIZMX),S(J21,9,J21)
      COMMON/LANCZO/LANCZ
C**TEMPORARY (DIMENSIONS)
      COMMON/SYMM/NVSYM,NWSYM,NSYM(10),ISYM(10,100),NTOT(10)
      COMMON/FILASS/IOUT,INP
      COMMON/MATRIX/IDUM,NVALR,KSTEP,KSIGN
      IELR=(IERX-1)/KSTEP+1
      IELL=(IELX-1)/KSTEP+1
      ILOFF=NVAL*(KROTL-1)
      IROFF=NVAL*(KROTR-1)
      IF(LANCZ)THEN
        KL0=KROTL
        KSIZE=(KROTL-KSTART)*NVAL
        JSIZE=ISIZE-KSIZE-(KSTART-1)*NVAL
        L0=KSIZE*JSIZE+KSIZE*(KSIZE+1)/2
      END IF
      IND=30+IABC
      IROT=IABC
C**READ BASIC INTEGRALS
      REWIND IND
      MSL=NSL
      IF(NSR.LT.NSL)MSL=NSR
      MSR=NSR
      IF(NSR.LT.NSL)MSR=NSL
      DO KSL=1,MSL
        KSIZEL=NTOT(KSL)
        DO KSR=KSL,NVSYM
          KSIZER=NTOT(KSR)
          IF(NSL.LE.NSR)THEN
            DO IX=1,KSIZER
              READ(IND)(YK(IY,IX),IY=1,KSIZEL)
            END DO
          ELSE
            DO IX=1,KSIZER
              READ(IND)(YK(IX,IY),IY=1,KSIZEL)
            END DO
          END IF
          IF(KSL.EQ.MSL.AND.KSR.EQ.MSR)GO TO 1000
        END DO
      END DO
1000  CONTINUE
C**RHS
C**CALL MATRIX MULT. ROUTINE MXMA TO SET UP TEMP(IY,I2)
C     CALL MXMA(YK(1,1),1,ISIZMX,CFSR(1,1,IELR),1,ISIZMX,
C    &      TEMP(1,1),1,ISIZMX,ISIZEL,ISIZER,NVAL)
      CALL DGEMM('N','N',ISIZEL,NVAL,ISIZER,1.0D0,YK(1,1),ISIZMX,
     &       CFSR(1,1,IELR),ISIZMX,0.0D0,TEMP,ISIZMX)
      DO IX=1,NVAL
        DO IY=1,NVAL
          YK(IY,IX)=0
        END DO
      END DO
C**LHS
C**CALL MXMA TO MULT. TEMP() BY LHS CFS
C     CALL MXMB(CFSL(1,1,IELL),ISIZMX,1,TEMP(1,1),1,ISIZMX,
C    &      YK(1,1),1,ISIZMX,NVAL,ISIZEL,NVAL)
      CALL DGEMM('T','N',NVAL,NVAL,ISIZEL,1.0D0,CFSL(1,1,IELL),
     &  ISIZMX,TEMP,ISIZMX,1.0D0,YK(1,1),ISIZMX)
C**MULTIPLY OFF-DIAGONAL BLOCK BY RELEVANT S-MATRIX ELEMENT
      IF(KROTL.EQ.KROTR)THEN
        IF(LANCZ)THEN
          IY=0
          JY=0
          IX0=0
          DO ILHS=1,NVAL
            IY=IY+1
            KR0=KL0
            IX=IX0
            DO IRHS=ILHS,JSIZE
              JRHS=IRHS-ILHS+1+JY
              IX=IX+1
              IF(KL0.EQ.KROTL.AND.KR0.EQ.KROTR)THEN
                XA(IRHS+L0)=XA(IRHS+L0)+YK(IY,IX)*S(JROTL,IROT,JROTR)*
     1          KSIGN
              END IF
              IF(MOD(JRHS,NVAL).EQ.0)THEN
                KR0=KR0+1
                IX=IX0
              END IF
            END DO
            L0=L0+JSIZE-ILHS
            JY=JY+1
            IX0=IX0+1
          END DO
        ELSE
          DO IX=1,NVAL
            JOFF=IROFF+IX-1
            J0=JOFF*(JOFF+1)/2+ILOFF
            DO IY=1,IX
              XA(J0+IY)=XA(J0+IY)+YK(IY,IX)*S(JROTL,IROT,JROTR)*KSIGN
            END DO
          END DO
        END IF
      ELSE
        IF(LANCZ)THEN
          IY=0
          JY=0
          DO ILHS=1,NVAL
            IY=IY+1
            KR0=KL0
            IX=0
            DO IRHS=ILHS,JSIZE
              JRHS=IRHS-ILHS+1+JY
              IX=IX+1
              IF(KL0.EQ.KROTL.AND.KR0.EQ.KROTR)THEN
                XA(IRHS+L0)=XA(IRHS+L0)+YK(IY,IX)*S(JROTL,IROT,JROTR)*
     1          KSIGN
              END IF
              IF(MOD(JRHS,NVAL).EQ.0)THEN
                KR0=KR0+1
                IX=0
              END IF
            END DO
            L0=L0+JSIZE-ILHS
            JY=JY+1
          END DO
        ELSE
          DO IX=1,NVAL
            JOFF=IROFF+IX-1
            J0=JOFF*(JOFF+1)/2+ILOFF
            DO IY=1,NVAL
              XA(J0+IY)=XA(J0+IY)+YK(IY,IX)*S(JROTL,IROT,JROTR)*KSIGN
            END DO
          END DO
        END IF
      END IF
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE ROTEL(MAXJ2,S,SX,JM,N,M1,M2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SX(JM,9,JM),S(JM,9,JM)
      COMMON/PRINT/IPRINT
      COMMON/FILASS/IOUT,INP
      DO K=1,JM
        DO M=1,9
          DO I=1,JM
            S(I,M,K)=0
          END DO
        END DO
      END DO
C**ON ENTRY MAXJ2 = 2*J
C**         JM = 2*J + 1
      JK=MAXJ2+1
      J=MAXJ2/2
      DO I2=1,JM
        DO K=1,9
          DO I1=1,JM
            SX(I1,K,I2)=0.D0
          END DO
        END DO
      END DO
C*********************************************
C     L1=L2
      L=-J-1
      DO I=1,JK
        L=L+1
C**NB....FUNCTION FOR MZ IS EXP(+I.K.GAMMA),NB.....NOT BRINK & SATCHLER
C**Jx**2
        SX(I,1,I)=0.5D0*(J*(J+1)-L*L)
C**Jy**2
        SX(I,2,I)=SX(I,1,I)
C**Jz**2
        SX(I,3,I)=L*L
C**Jz TAKE ACCOUNT OF 'i'
        SX(I,9,I)=L
      END DO
C*********************************************
C     L1=L2+1,  L1=L2-1
      L1=-J-1
      JK1=JK-1
      DO I1=1,JK1
        I2=I1+1
        L1=L1+1
        L2=L1+1
C**Jx ALWAYS NEGATIVE
        SX(I1,7,I2)=-0.5D0*SQRT((J+L2)*(J-L2+1.D0))
        SX(I2,7,I1)=+0.5D0*SQRT((J-L1)*(J+L1+1.D0))
C**Jy INITIALLY POSITIVE
        SX(I1,8,I2)=-SX(I1,7,I2)
        SX(I2,8,I1)=SX(I2,7,I1)
C**JxJz + JzJx INITIALLY LIKE POSITIVE Jx
        SX(I1,4,I2)=-(2*L2-1)*SX(I1,7,I2)
        SX(I2,4,I1)=-(2*L1+1)*SX(I2,7,I1)
C**JyJz + JzJy ALWAYS LIKE POSITIVE Jy
        SX(I1,5,I2)=(2*L2-1)*SX(I1,8,I2)
        SX(I2,5,I1)=(2*L1+1)*SX(I2,8,I1)
      END DO
C*********************************************
C     L1=L2+2,   L1=L2-2
      L1=-J-1
      JK2=JK-2
      IF(JK2.GT.0)THEN
        DO I1=1,JK2
          I2=I1+2
          L1=L1+1
          L2=L1+2
C**Jx**2
          SX(I1,1,I2)=-0.25D0*SQRT((J-L2+1.D0)*(J-L2+2.D0)*(J+L2)*
     1    (J+L2-1.D0))
          SX(I2,1,I1)=-0.25D0*SQRT((J+L1+1.D0)*(J+L1+2.D0)*(J-L1)*
     1    (J-L1-1.D0))
C**Jy**2
          SX(I1,2,I2)=-SX(I1,1,I2)
          SX(I2,2,I1)=-SX(I2,1,I1)
C**JxJy + JyJx INITIALLY POSITIVE
          SX(I1,6,I2)=0.5D0*SQRT((J+L2-1.D0)*(J-L2+2.D0)*(J+L2)*
     1    (J-L2+1.D0))
          SX(I2,6,I1)=-0.5D0*SQRT((J-L1-1.D0)*(J+L1+2.D0)*(J-L1)*
     1    (J+L1+1.D0))
        END DO
      END IF
C***********************************************************
C***********************************************************
      DO I=1,9
C**Jy AND [JxJz]+ TAKE ACCOUNT OF 'i'
      IS48=1
      IF(I.EQ.4.OR.I.EQ.8)IS48=-1
C**[JxJy]+ TAKE ACCOUNT OF 'i'
      IS6=1
      IF(I.EQ.6)IS6=-1
C**TAKE ACCOUNT OF ANTI-COMMUTATION
      IS789=1
      IF(I.EQ.7.OR.I.EQ.8.OR.I.EQ.9)IS789=-1
C**INDEX FOR D(0,0)J
        J0=J+1
C*******************
        I1=1
        SQ2=SQRT(2.D0)
        I2=I1
C**<D(0,0)J/O/D(0,0)J>
        S(I1,I,I2)=SX(J0,I,J0)*M1
        IF(J.EQ.0)GO TO 5000
C**<D(0,0)J/O/D(0,K)J>
        DO L=1,J
          IF(L.GT.2)GO TO 2000
          I2=I2+1
          IS=(-1)**(N+L)
          J1=J0+L
          J2=J0-L
C**Jy AND [JxJz]+ TAKE ACCOUNT OF 'i'
C**[JxJy]+ TAKE ACCOUNT OF 'i'
          S(I1,I,I2)=IS48*IS6*(SX(J0,I,J1)+SX(J0,I,J2)*IS)/SQ2
          I2=I2+1
          S(I1,I,I2)=(SX(J0,I,J1)-SX(J0,I,J2)*IS)/SQ2
2000      CONTINUE
        END DO
        I1=2
        ISTART=I1
C**<D(0,K+)J/O/D(0,K+)J>
C**<D(0,K+)J/O/D(0,K-)J>
C**<D(0,K-)J/O/D(0,K+)J>
C**<D(0,K-)J/O/D(0,K-)J>
        DO L=1,J
          I2=I1
          IS=(-1)**(N+L)
          J1=J0+L
          J2=J0-L
          S(I1,I,I2)=(M1*(SX(J1,I,J1)+SX(J2,I,J2))+IS*M2*(SX(J1,I,J2)
     1    +SX(J2,I,J1)))/2.D0
          I2=I2+1
          S(I1,I,I2)=(SX(J1,I,J1)-SX(J2,I,J2)-IS*(SX(J1,I,J2)-
     1    SX(J2,I,J1)))/2.D0
          S(I2,I,I1)=(SX(J1,I,J1)-SX(J2,I,J2)+IS*(SX(J1,I,J2)-
     1    SX(J2,I,J1)))/2.D0
          I1=I1+1
          I2=I1
          S(I1,I,I2)=(M1*(SX(J1,I,J1)+SX(J2,I,J2))-IS*M2*(SX(J1,I,J2)
     1    +SX(J2,I,J1)))/2.D0
          I1=I1+1
        END DO
        I1=ISTART
        L1=1
3000    I2=I1+2
C**<D(0,K)J/O/D(0,K+1)J>
        L2=L1+1
        IF(L2.GT.J)GO TO 5000
        IS1=(-1)**(N+L1)
        IS2=(-1)**(N+L2)
        J1=J0+L1
        J2=J0-L1
        K1=J0+L2
        K2=J0-L2
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I1=I1+1
        I2=I2-1
C**Jy AND [JxJz]+ TAKE ACCOUNT OF 'i'
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))*IS48/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))/2.D0
C**<D(0,K)J/O/D(0,K+2)J>
        I1=I1-1
        I2=I2+1
        L2=L2+1
        IF(L2.GT.J)GO TO 4000
        IS2=(-1)**(N+L2)
        K1=J0+L2
        K2=J0-L2
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)+
     1  IS1*SX(J2,I,K1))/2.D0
        I1=I1+1
        I2=I2-1
C**[JxJy]+ TAKE ACCOUNT OF 'i'
        S(I1,I,I2)=(SX(J1,I,K1)-IS1*IS2*SX(J2,I,K2)+IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))*IS6/2.D0
        I2=I2+1
        S(I1,I,I2)=(SX(J1,I,K1)+IS1*IS2*SX(J2,I,K2)-IS2*SX(J1,I,K2)-
     1  IS1*SX(J2,I,K1))/2.D0
4000    CONTINUE
        I1=I1+1
        L1=L1+1
        GO TO 3000
5000    CONTINUE
C**TAKE ACCOUNT OF ANTI-COMMUTATION
        DO IX=1,JM
        DO JX=1,JM
          S(JX,I,IX)=S(JX,I,IX)*IS789
        END DO
        END DO
      END DO
      IF(IPRINT.GT.2)WRITE(IOUT,*)'ROTATION MATRIX'
      DO I=1,JM
      DO J=1,JM
      IF(IPRINT.GT.2)WRITE(IOUT,500)I,J,(S(I,K,J),K=1,9)
500   FORMAT(/,1X,2I3,6(5X,F12.8))
      END DO
      END DO
      RETURN
      END
