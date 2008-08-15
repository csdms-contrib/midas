C***********************SUBROUTINE GSDIST*************************
C------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES THE PHINORMAL GRAIN SIZE
C                   DISTRIBUTION FOR EACH GRAIN DENSITY
C------------------------------------------------------------------
      SUBROUTINE GSDIST(J)
      IMPLICIT NONE
      REAL*8 DIMID(20,3),SIGMA(3)
      REAL*8 DIFEND(20,3),DIMIN(3)
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100),
     #AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 AWPM(20,3,100),D50M(100)
      REAL*8 MEAN(3),STDEV(3),PROPTJ(3)                              
      INTEGER*4 NSIZE(3),NSIGMA,M,N
C----------------------------------------------------------------------
C        COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM22/AWPM,D50M
      COMMON /COM24/MEAN,STDEV,PROPTJ
      COMMON /COM28/NSIZE,NSIGMA,M,N
C---------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 FO,FI,H,T,P,A1,A2,A3,A4,A5,H1,PISQR
      REAL*8 ALIMIT,BLIMIT,AWP(20,3),CONV,CLLENGTH,DISTMIN1,DIST
      INTEGER*4 I,J,JN,FLAG
C----------------------------------------------------------------------
C     THIS SUBROUTINE CALCULATES THE ERROR FUNCTION USING THE
C     FOLLOWING ALGORITHM:
C     ERF(X)= 1-(A1*T+A2*T^2+A3*T^3+A4*T^4+A5*T^5)*EXP(-X^2)
C     T=1/(1+P*X)
C     PROBABILITIES ( I.E. AWP'S ) OF A CERTAIN GRAIN OCCURRING
C     IN A SIZE INTERVAL ARE CALCULATED BY SUBTRACTING ERF'S
C
C     CONV=CONVERSION FACTOR BETWEEN LOG BASE E AND LOG BASE 2
C     DETERMINE NECESSARY COEFFICIENTS                                -
C----------------------------------------------------------------------
      P = .2316419
      A1 = .31938153
      A2 = -0.356563782
      A3 = 1.781477937
      A4 = -1.821255978
      A5 = 1.330274427
      CONV = 0.693147181
      PISQR = SQRT(2.0*3.1415926536)
C----------------------------------------------------------------------
C     MEAN AND STANDARD DEVIATION SHOULD BE IN PHI UNITS
C     CLASSES RANGE FROM MEAN-3*STDEV TO MEAN+3*STDEV
C     ALIMIT=SMALLEST GRAIN SIZE (=LARGEST PHI SIZE )
C     BLIMIT=LARGEST GRAIN SIZE
C----------------------------------------------------------------------
      ALIMIT = MEAN(J) + 3.0*STDEV(J)
      DIMIN(J) = EXP((-CONV*ALIMIT))/1000.
      BLIMIT = MEAN(J) - 3.0*STDEV(J)
      CLLENGTH = (BLIMIT-ALIMIT)/FLOAT(NSIZE(J))
      DISTMIN1 = (ALIMIT-MEAN(J))/STDEV(J)
      H1 = EXP((-DISTMIN1*DISTMIN1/2.))/PISQR
      T = 1/(1+P*DISTMIN1)
      H = ((((A5*T+A4)*T+A3)*T+A2)*T+A1)*T
      FO = H*H1
      DO 2 I = 1, NSIZE(J)
         FLAG = 0
         DIST = ((ALIMIT+I*CLLENGTH)-MEAN(J))/STDEV(J)
         DIMID(I,J) = EXP((-CONV*(ALIMIT+(I-0.5)*CLLENGTH)))/1000.
         DIFEND(I,J) = EXP((-CONV*(ALIMIT+I*CLLENGTH)))/1000.
         DIMID(I,J) = INT(DIMID(I,J)*1000000.+0.5)/1000000.
         DIFEND(I,J) = INT(DIFEND(I,J)*1000000.+0.5)/1000000.
         H1 = EXP((-DIST*DIST/2.))/PISQR
         IF (DIST .LT. 0.) THEN
            DIST = -DIST
            FLAG = 1
         ENDIF
         T = 1/(1+P*DIST)
         H = ((((A5*T+A4)*T+A3)*T+A2)*T+A1)*T
         FI = H*H1
         IF (FLAG .EQ. 1) FI = 1 - FI
         AWP(I,J) = FI - FO
         FO = FI
         DO 1 JN = 1, M
            AWPOM(I,J,JN) = AWP(I,J)*PROPTJ(J)
            AWPM(I,J,JN) = AWP(I,J)*PROPTJ(J)
    1    CONTINUE
    2 CONTINUE
      RETURN 
      END
