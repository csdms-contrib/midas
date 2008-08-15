C**************************SUBROUTINE LOGDIST**************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES THE LOGRITHMIC VELOCITY DISTRIBUTION
C                   CALLLED FROM TRCALC
C----------------------------------------------------------------------
      SUBROUTINE LOGDIST
      IMPLICIT NONE
      REAL*8 YMAX(100),DELY(100),VEL(30,100),EPSY(30,100),PSI,BETA 
      REAL*8 AWPM(20,3,100),D50M(100)                             
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),
     #SVELUM(100)
      REAL*8   STRDATA(10,100)                                  
      INTEGER*4 NSIZE(3),NSIGMA,M,N                            
C---------------------------------------------------------------------
C     COMMON DECLARATION
C---------------------------------------------------------------------
      COMMON /COM19/YMAX,DELY,VEL,EPSY,PSI,BETA
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM26/STRDATA
      COMMON /COM28/NSIZE,NSIGMA,M,N
C---------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C---------------------------------------------------------------------
      REAL*8 Y,V
      INTEGER*4 L,JN
C---------------------------------------------------------------------
      DO 2 JN = 2, M
         DELY(JN) = (DCHM(JN)-2.0*D50M(JN))/(N-1)
         YMAX(JN) = 2.0*D50M(JN)
         DO 1 L = 1, N
            V = STRDATA(5,JN)
            Y = DCHM(JN) - (L-1)*DELY(JN)
            IF (L .EQ. N) Y = YMAX(JN)
            VEL(L,JN) = (SVELUM(JN)/0.4)*(DLOG(Y/DCHM(JN))+1.) + V
C----------------------------------------------------------------------
C CHECK VALIDITY OF THIS ASSUMPTION
C----------------------------------------------------------------------
            IF (VEL(L,JN) .LT. SVELUM(JN)) VEL(N,JN) = SVELUM(JN)
    1    CONTINUE
    2 CONTINUE
      RETURN 
      END
