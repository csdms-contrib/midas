C*****************************SVELA************************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE:CALCULATES SHEAR VELOCITY ASSOCIATED WITH GRAIN
C                  ROUGHNESS
C----------------------------------------------------------------------
      SUBROUTINE SVELA
      IMPLICIT NONE
      REAL*8 RHO,G                                                 
      REAL*8 AWPM(20,3,100),D50M(100)                            
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100), 
     #SVELUM(100)
      REAL*8   STRDATA(10,100)                                  
      INTEGER*4 NSIZE(3),NSIGMA,M,N                            
C----------------------------------------------------------------------
C            DECLARATION STATEMENTS--COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM26/STRDATA
      COMMON /COM28/NSIZE,NSIGMA,M,N
C---------------------------------------------------------------------- -
C            DECLARATION STATEMENT--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 AK,DR1,DR2,DR3,DR,F3,RCH,SF,V
      INTEGER*4 JN
C----------------------------------------------------------------------
C  CALCULATE INITIAL DATA (RECTANGULAR CHANNEL SECTION IS ASSUMED)
C----------------------------------------------------------------------
      DO 5 JN = 2, M
         RCH = STRDATA(1,JN)
         V = STRDATA(5,JN)
         SF = STRDATA(3,JN)
C----------------------------------------------------------------------
C  CALCULATE REDUCED HYDRAULIC RADIUS
C----------------------------------------------------------------------
         AK = 2.5*D50M(JN)
         DR1 = 2.0*RCH
         DR2 = 0.1*RCH
         DR3 = 0.5*(DR1+DR2)
    1    CONTINUE
         F3 = 6.25 + 2.5*DLOG(DR3/AK) - V/SQRT(G*DR3*SF)
         IF (ABS(F3) .LE. 0.001) GO TO 4
         IF (DR1 .EQ. DR2) GO TO 4
         IF (F3 .GT. 0) GO TO 2
         IF (F3 .EQ. 0) GO TO 4
         DR2 = DR3
         GO TO 3
    2    CONTINUE
         DR1 = DR3
    3    CONTINUE
         DR3 = 0.5*(DR1+DR2)
         GO TO 1
    4    CONTINUE
         DR = DR3
         DR = AMIN1(RCH,DR)
C----------------------------------------------------------------------
C  CALCULATE SHEAR VELOCITY
C----------------------------------------------------------------------
         SVELM(JN) = SQRT(G*DR*SF)
         SVELUM(JN) = SQRT(G*RCH*SF)
C----------------------------------------------------------------------
C     CALCULATE THE BED SHEAR STRESS
C----------------------------------------------------------------------
         TOXM(JN) = RHO*G*SF*DR
    5 CONTINUE
      RETURN 
      END
