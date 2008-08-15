C*****************************SUBROUTINE TURB*************************
C----------------------------------------------------------------------
C THIS SUBROUTINE: CALCULATES THE GAUSSIAN OR LOG-GAUSSIAN DISTRIBUTION
C                  OF INSTANTANEOUS SHEAR STRESSES ON THE BED, GIVEN A
C                  MEAN AND COEFFICIENT OF VARIATION
C----------------------------------------------------------------------
      SUBROUTINE TURB
      IMPLICIT NONE
      REAL*8 RHO,G                                
      REAL*8 PROPT(50)  
      REAL*8 SHVEL(50,100),TO(50,100)            
      REAL*8 COVAR,A,B,CON1,CON2                
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),
     #SVELUM(100)
      INTEGER*4 NK
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*4 NSIZE(3),NSIGMA,M,N                             
C----------------------------------------------------------------------
C            COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM6/PROPT,NK
      COMMON /COM7/SHVEL,TO
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #OPTBINP,OPTSINP,OPTERO
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C            DECLARATION STATEMENTS-EXTERNAL VARIABLES NOT IN COMMON
C----------------------------------------------------------------------
        REAL*8 TOMEAN
        INTEGER*4 JN
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 CONV,DISTMIN1,DIST,DELTO,DELTOH,TOLO
      REAL*8 FO,FI,H,T,P,A1,A2,A3,A4,A5,STDEV,H1,PISQR
      INTEGER*4 K,FLAG
C----------------------------------------------------------------------
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
C    INITIALIZE INCREMENT IN INST BED SHEAR STRESS AND MIDPOINT
C    VALUES OF BED SHEAR STRESS RANGES
C----------------------------------------------------------------------
      DO 4 JN = 2, M
C----------------------------------------------------------------------
C       IF COVAR=0.333 THEN TO(K,JN) VARIES FROM 0 TO TOXM(JN)+3SIGMA
C       IF COVAR>0.333 THEN THE CALCULATED TOXM(JN) IS NO LONGER
C       EQUAL TO THE INPUT TOXM(JN) DUE TO TRUNCATION
C----------------------------------------------------------------------
         IF (OPTURB .EQ. 1) THEN
            TOMEAN = TOXM(JN)
         ELSE
            TOMEAN = LOG(TOXM(JN))
         ENDIF
         STDEV = COVAR*TOMEAN
         DELTO = (STDEV*6.0D0)/FLOAT(NK)
         DELTOH = DELTO*0.5D0
         IF (COVAR .LE. 0.333D0) THEN
            TO(1,JN) = TOMEAN - STDEV*3.0D0 + DELTOH
         ELSE
            TO(1,JN) = DELTOH
         ENDIF
         TOLO = TO(1,JN) - DELTOH + 0.00001
         DO 1 K = 2, NK
            TO(K,JN) = TO(K-1,JN) + DELTO
    1    CONTINUE
         IF (OPTURB .EQ. 2) THEN
            DO 2 K = 1, NK
               TO(K,JN) = EXP(TO(K,JN))
    2       CONTINUE
         ENDIF
C----------------------------------------------------------------------
C    CORRESPONDING VALUE OF SHEAR VELOCITY
C----------------------------------------------------------------------
         DO 3 K = 1, NK
            SHVEL(K,JN) = DSQRT(TO(K,JN)/RHO)
    3    CONTINUE
    4 CONTINUE
C----------------------------------------------------------------------
C     GAUSSIAN DISTRIBUTION OF BED SHEAR STRESS
C----------------------------------------------------------------------
      DISTMIN1 = (TOLO-TOMEAN)/STDEV
      H1 = EXP((-DISTMIN1*DISTMIN1/2.))/PISQR
      DISTMIN1 = -DISTMIN1
      T = 1/(1+P*DISTMIN1)
      H = ((((A5*T+A4)*T+A3)*T+A2)*T+A1)*T
      FO = H*H1
      DO 5 K = 1, NK
         FLAG = 0
         DIST = ((TOLO+K*DELTO)-TOMEAN)/STDEV
         H1 = EXP((-DIST*DIST/2.))/PISQR
         IF (DIST .LT. 0.) THEN
            DIST = -DIST
            FLAG = 1
         ENDIF
         T = 1/(1+P*DIST)
         H = ((((A5*T+A4)*T+A3)*T+A2)*T+A1)*T
         FI = 1 - H*H1
         IF (FLAG .EQ. 1) FI = 1 - FI
         PROPT(K) = FI - FO
         FO = FI
    5 CONTINUE
      RETURN 
      END
