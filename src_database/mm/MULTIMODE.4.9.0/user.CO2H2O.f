C**************************************************************
C**************************************************************
      SUBROUTINE USERIN
      use pes_shell
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/FILASS/IOUT,INP
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/DUMP/JDUMP,IDUMP,KDUMP,MDUMP,LDUMP
      COMMON/ROTPAR/RX1,RX2,RX3,THX1,THX2,TAUX
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
      call pes_init() 
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETPOT(V,NATOM,XX,RR)
      use pes_shell
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),RR(NATOM,NATOM),RRR(3)
      DIMENSION QQ(6)
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/ECKCNT/ICNT,INTC
      COMMON/FILASS/IOUT,INP
      DIMENSION dd(3,NATOM)
      

      do i=1,NATOM
        dd(:,i)=XX(i,:)
      end do
       
       V=f(dd) + 0.004798548260072554 
     
      
 
        if (V .LE. -0.1) then
        open(unit=11,status='unknown',file='err.xyz')
        write (11,*) 6
        write (11,*)  V*219474.63
        write (11,*) 'O',XX(1,:)*0.52918
        write (11,*) 'O',XX(2,:)*0.52918
        write (11,*) 'O',XX(3,:)*0.52918
        write (11,*) 'H',XX(4,:)*0.52918
        write (11,*) 'H',XX(5,:)*0.52918
        write (11,*) 'C',XX(6,:)*0.52918
        end if
      
       CALL BONDS(NATOM,RR,XX)


       RETURN
       END

C****************************************************************
C****************************************************************
      SUBROUTINE GETQPT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETAPT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE GETDIP
      RETURN
      END
C*****************************
C*****************************
      SUBROUTINE GETQDT
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE MINPOT(TAU,NATOM,XR,RR,TAUE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XR(NATOM,3),RR(NATOM,NATOM)
      DIMENSION WRK(200),SOL(5),F(5)
      COMMON/ROTPAR/RX1,RX2,RX3,THX1,THX2,TAUX
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
      COMMON/FUNDAM/WAVENM,ATTOJ,BOHR,ELMASS,RAD
      COMMON/FILASS/IOUT,INP
C**RETURNS WITH CARTESIAN COORDINATES IN XR CORRESPONDING TO
C**MINIMUM ENERGY AT FIXED TAU
C*****************************************************
      EXTERNAL DERIV,MONIT
C*****************************************************
      IFAIL=1
      MFIT=5
C**TEMPORARY
C     STEP=10.D0
      STEP=5.D0
      XTOL=1.D-8
      FTOL=1.D-15
      MAXCAL=1000
      KPRINT=1000
C**TEMPORARY
C     WRITE(IOUT,*)'*******************************'
      IWXY=2*MFIT*(MFIT+MFIT)+2*MFIT+5*MFIT
      IF(IWXY.GT.200)THEN
        WRITE(IOUT,*)
        WRITE(IOUT,*)'WRK TOO SMALL'
        STOP 'WRK TOO SMALL IN MINPOT'
      END IF
C**SET CURRENT TAU
      TAUX=TAU+PHIE
      TAUE=TAUX
C**ANSWER NEAR LAST STRUCTURE, SO INITIAL GUESS SET TO IT
      SOL(1)=RX1
      SOL(2)=RX2
      SOL(3)=RX3
      SOL(4)=THX1
      SOL(5)=THX2
      CALL E04FBF(MFIT,MFIT,SOL,F,SUMSQ,FTOL,XTOL,STEP,WRK,200,
     1DERIV,MONIT,KPRINT,MAXCAL,IFAIL)
C     WRITE(IOUT,*)
C     WRITE(IOUT,*)'MINPOT IFAIL = ',IFAIL
      IF(IFAIL.NE.0.AND.IFAIL.NE.3)THEN
        STOP 'ERROR IN E04FBF'
      END IF
C**SAVE THIS STRUCTURE FOR NEXT TIME
      RX1=SOL(1)
      RX2=SOL(2)
      RX3=SOL(3)
      THX1=SOL(4)
      THX2=SOL(5)
C**TEMPORARY
C     WRITE(IOUT,*)'TAU: ',TAUX*RAD
C     WRITE(IOUT,*)'GEOM: ',RX1*BOHR,RX2*BOHR,RX3*BOHR,
C    1THX1*RAD,THX2*RAD
C**TEMPORARY
C**SET CARTESIAN COORDINATES FOR THIS STRUCTURE
      CT=DCOS(TAUX/2)
      IF(DABS(CT).LT.1.D-10)CT=0.D0
      ST=DSIN(TAUX/2)
      IF(DABS(ST).LT.1.D-10)ST=0.D0
      XR(1,1)=RX1*DSIN(THX1)*CT
      XR(1,2)=RX1*DSIN(THX1)*ST
      XR(1,3)=-RX3/2+RX1*DCOS(THX1)
      XR(2,1)=RX2*DSIN(THX2)*CT
      XR(2,2)=-RX2*DSIN(THX2)*ST
      XR(2,3)=RX3/2-RX2*DCOS(THX2)
      XR(3,1)=0.D0
      XR(3,2)=0.D0
      XR(3,3)=-RX3/2
      XR(4,1)=0.D0
      XR(4,2)=0.D0
      XR(4,3)=RX3/2
      CALL GETPOT(V,NATOM,XR,RR)
C     WRITE(IOUT,*)'V = ',V*WAVENM
      RETURN
      END
C****************************************************************
C****************************************************************
      SUBROUTINE DERIV(M,N,X,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),F(M)
      DIMENSION QQ(6),TERM(6)
      COMMON/ROTPAR/RX1,RX2,RX3,THX1,THX2,TAUX
      COMMON/POTPAR/RE1,RE2,RE3,THE1,THE2,PHIE,CTERM(200),NTERM,
     1KTERM(200,6),ITERM(200,6),JTERM(200,6),
     2NDIP(3),KDIP(100,6,3),IDIP(100,6,3),JDIP(100,6,3),CDIP(100,3)
      COMMON/FILASS/IOUT,INP
C**RESET INTERNAL COORDINATES
      R1=X(1)
      R2=X(2)
      R3=X(3)
      TH1=X(4)
      TH2=X(5)
      QQ(1)=(R1-RE1)/R1
      QQ(2)=(R2-RE2)/R2
      QQ(3)=(R3-RE3)/R3
      QQ(4)=TH1-THE1
      QQ(5)=TH2-THE2
      QQ(6)=TAUX
C**FIRST DERIVATIVES
      DO I=1,5
        F(I)=0.D0
      END DO
      DO I=1,NTERM
        DO ID=1,5
C**GET DERIVATIVE WRT COORDINATE ID
          FACT=0.D0
          DO K=1,6
            TERM(K)=1.D0
          END DO
C**MAXIMUM OF 6 MODES COUPLED
          DO J=1,6
            IF(ITERM(I,J).NE.0)THEN
              K=ITERM(I,J)
              L=JTERM(I,J)
              IF(K.EQ.ID)THEN
                TERM(K)=CTERM(I)
                FACT=L*QQ(K)**(L-1)
                IF(K.EQ.1)FACT=FACT*RE1/(R1*R1)
                IF(K.EQ.2)FACT=FACT*RE2/(R2*R2)
                IF(K.EQ.3)FACT=FACT*RE3/(R3*R3)
              ELSE
                TERM(K)=QQ(K)**L
                IF(K.EQ.6)TERM(K)=DCOS(L*TAUX)
              END IF
            END IF
          END DO
          Z=1.D0
          DO K=1,6
            Z=Z*TERM(K)
          END DO
          F(ID)=F(ID)+Z*FACT
        END DO
      END DO
CCCC  WRITE(IOUT,*)(F(I),I=1,5)
      RETURN
      END

