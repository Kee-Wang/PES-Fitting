C***********************************************************
C***********************************************************
C**USER-SUPPLIED
C***********************************************************
C***********************************************************
      SUBROUTINE ROTATE(NATOM,X0,OMEGA,NMODE,XL,WAVENM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X0(NATOM,3),OMEGA(NMODE),XL(NATOM,NMODE,3)
C*************************************************
C**ROUTINE TO ROTATE NORMAL COORDINATE VECTORS
C*************************************************
C**TEMPORARY
CCC   IF(DABS(X0(1,1))-DABS(X0(1,2)).LT.1.D-6)THEN
C**FIND E SYMMETRIES
        DO J=1,NMODE
          IF(J.LT.NMODE)THEN
            IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-3)THEN
C**ROTATION ANGLE ZERO
              IF(DABS(XL(4,J,3)-XL(5,J,3)).LT.1.D-5)GO TO 100
C**ROTATION ANGLE INFINITE
              IF(DABS(XL(4,J+1,3)-XL(5,J+1,3)).LT.1.D-5)GO TO 100
              TH=(XL(4,J,3)-XL(5,J,3))/
     1           (-XL(4,J+1,3)+XL(5,J+1,3))
              TH=DATAN(TH)
              CTH=DCOS(TH)
              STH=DSIN(TH)
              DO I=1,NATOM
                DO K=1,3
                  TEMP1=CTH*XL(I,J,K)+STH*XL(I,J+1,K)
                  TEMP2=-STH*XL(I,J,K)+CTH*XL(I,J+1,K)
                  XL(I,J,K)=TEMP1
                  XL(I,J+1,K)=TEMP2
                END DO
              END DO
100           CONTINUE
            END IF
          END IF
        END DO
C**FIND E SYMMETRIES
C**********************************************
C**FIND T SYMMETRIES
        DO J=1,NMODE
          IF(J.LT.NMODE-1)THEN
            IF(DABS(OMEGA(J)-OMEGA(J+1))*WAVENM.LT.1.D-3)THEN
            IF(DABS(OMEGA(J)-OMEGA(J+2))*WAVENM.LT.1.D-3)THEN
C**ROTATION ANGLE ZERO
              IF(DABS(XL(2,J,3)-XL(3,J,3)).LT.1.D-5)GO TO 200
C**ROTATION ANGLE INFINITE
              IF(DABS(XL(2,J+1,3)-XL(3,J+1,3)).LT.1.D-5)GO TO 200
              TH=(XL(2,J,3)-XL(3,J,3))/
     1           (-XL(2,J+1,3)+XL(3,J+1,3))
              TH=DATAN(TH)
              CTH=DCOS(TH)
              STH=DSIN(TH)
              DO I=1,NATOM
                DO K=1,3
                  TEMP1=CTH*XL(I,J,K)+STH*XL(I,J+1,K)
                  TEMP2=-STH*XL(I,J,K)+CTH*XL(I,J+1,K)
                  XL(I,J,K)=TEMP1
                  XL(I,J+1,K)=TEMP2
                END DO
              END DO
200           CONTINUE
            END IF
            END IF
          END IF
        END DO
C**FIND T SYMMETRIES
CCC   END IF
C**TEMPORARY
      RETURN
      END
C***********************************************************
C***********************************************************
      SUBROUTINE RTGEOM(NATOM,XX,E,WAVENM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XX(NATOM,3),E(3)
      COMMON/MOMI/XK(3,3),XMU(3,3)
C*************************************************
C**ROUTINE TO ROTATE VECTORS OF PRINCIPAL AXIS SYSTEM
C*************************************************
C**TEMPORARY
C**FIND E SYMMETRIES
C**TEMPORARY
CCC   IF(DABS(XX(4,1))-DABS(XX(4,2)).LT.1.D-6)THEN
CCC     DO J=1,3
CCC       IF(J.LT.3)THEN
CCC         IF(DABS(1/E(J)-1/E(J+1))*WAVENM.LT.1.D-3)THEN
CCC           IF(DABS(XK(1,J+1)-XK(2,J+1)).LT.1.D-8)GO TO 999
CCC           TH=(XK(1,J)-XK(2,J))/
CCC  1           (-XK(1,J+1)+XK(2,J+1))
CCC           TH=DATAN(TH)
CCC           CTH=DCOS(TH)
CCC           STH=DSIN(TH)
CCC           DO I=1,3
CCC             TEMP1=CTH*XK(I,J)+STH*XK(I,J+1)
CCC             TEMP2=-STH*XK(I,J)+CTH*XK(I,J+1)
CCC             XK(I,J)=TEMP1
CCC             XK(I,J+1)=TEMP2
CCC           END DO
999           CONTINUE
CCC         END IF
CCC       END IF
CCC     END DO
CCC   END IF
C**TEMPORARY
C**FIND E SYMMETRIES
      RETURN
      END
