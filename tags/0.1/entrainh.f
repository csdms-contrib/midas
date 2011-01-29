C************************SUBROUTINE ENTRAINH***************************
C----------------------------------------------------------------------
C THIS SUBROUTINE: CALCULATES THE CRITICAL SHIELDS THETA FOR THE MEDIAN
C                  SIZE OF THE DISTRIBUTION AND THEN CALCULATES THE
C                  CRITICAL SHEAR STRESS FOR THE I,J FRACTION USING
C                  A HIDING FUNCTION
C----------------------------------------------------------------------
      SUBROUTINE ENTRAINH
      IMPLICIT NONE
      REAL*8 RHO,G                                            
      REAL*8 VISKIN,THETAC                                   
      REAL*8 DIMID(20,3),SIGMA(3)                           
      REAL*8 DIFEND(20,3),DIMIN(3)                         
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)     
      REAL*8 THET50(100),VISC                            
      REAL*8 COVAR,A,B,CON1,CON2                        
      REAL*8 AWPM(20,3,100),D50M(100)                  
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*4 NSIZE(3),NSIGMA,M,N                   
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM2/VISKIN,THETAC
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM12/THET50,VISC
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #              OPTBINP,OPTSINP,OPTERO
      COMMON /COM22/AWPM,D50M
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 TCI,REG,THETLN,RELN,D50
      INTEGER*4 I,J,JN
C----------------------------------------------------------------------
C     INITIALIZE APPROXIMATION TO CRITICAL BED SHEAR STRESS
C----------------------------------------------------------------------
      DO 6 JN = 2, M
         D50 = D50M(JN)
         IF (D50 .LE. 0.0008D0) THEN
            TCI = 0.1D0*17.7D0*(D50*100.0D0)**0.59D0
         ELSE
            TCI = 0.1D0*65.62D0*(D50*100.0D0)**1.11D0
         ENDIF
         REG = DSQRT(TCI/RHO)*D50/VISKIN
         IF (REG .LE. 1.0D0) THEN
C----------------------------------------------------------------------
C     VALUE FOR HYDRAULICALLY SMOOTH BED
C----------------------------------------------------------------------
            WRITE (6, 100) 
  100            FORMAT(' '//'WARNING:SHIELDS THETA IN THE HIDING'/
     #'FUNCTION COMES FROM THE HYDRAULICALLY SMOOTH REGIME')
            THETLN = DLOG(THETAC)
         ELSE IF (REG .LE. 60.0) THEN
C----------------------------------------------------------------------
C     VALUE FOR HYDRAULICALLY TRANSITIONAL BED
C----------------------------------------------------------------------
            RELN = DLOG(REG)
            THETLN = (-2.26D0) - 0.905D0*RELN + 0.168*RELN*RELN
         ELSE
C----------------------------------------------------------------------
C     VALUE FOR HYDRAULICALLY ROUGH BED
C----------------------------------------------------------------------
            THETLN = DLOG(THETAC)
         ENDIF
         THET50(JN) = DEXP(THETLN)
C----------------------------------------------------------------------
C     CALCULATE THE TC FOR MIXTURES WITH HIDING,USING KOMAR EQN 6
C     OR EGIAZAROFF, DEPENDING UPON OPTFUNC
C----------------------------------------------------------------------
         DO 5 J = 1, NSIGMA
            DO 4 I = 1, NSIZE(J)
               GO TO (1,2) OPTFUNC
    1          CONTINUE
               TC(I,J,JN) = THET50(JN)*G*(SIGMA(J)-RHO)*DIMID(I,J)**(1.+
     #            CON1)*D50**(-CON1)
 
               GO TO 3
    2          CONTINUE
               IF (DIMID(I,J)/D50 .LE. 1./CON2) THEN
                  TC(I,J,JN) = 100.0
               ELSE
                  TC(I,J,JN) = (SIGMA(J)-RHO)*G*DIMID(I,J)*THET50(JN)*
     #               DLOG10(CON2)**2./DLOG10(CON2*DIMID(I,J)/D50)**2.
               ENDIF
    3          CONTINUE
               SHVELC(I,J,JN) = SQRT(TC(I,J,JN)/RHO)
    4       CONTINUE
    5    CONTINUE
    6 CONTINUE
      RETURN 
      END
