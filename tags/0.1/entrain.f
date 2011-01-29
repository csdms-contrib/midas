C*******************************SUBROUTINE ENTRAIN*********************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES THE CRITICAL SHEAR STRESS OF THE
C                   MEDIAN SIZE OF EACH SIZE-DENSITY FRACTION USING
C                   YALIN AND KARAHAN
C----------------------------------------------------------------------
      SUBROUTINE ENTRAIN
      IMPLICIT NONE
      REAL*8 RHO,G                                      
      REAL*8 VISKIN,THETAC                             
      REAL*8 DIMID(20,3),SIGMA(3)                     
      REAL*8 DIFEND(20,3),DIMIN(3)                   
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)
      INTEGER*4 NSIZE(3),NSIGMA,M,N                 
C---------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C---------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM2/VISKIN,THETAC
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C---------------------------------------------------------------------
      REAL*8 TCI,FTC,DFTC,TCIN,REG,THETLN,RELN,TC1(20,3)
      INTEGER*4 I,J,JN
C---------------------------------------------------------------------
      DO 5 J = 1, NSIGMA
         DO 4 I = 1, NSIZE(J)
C----------------------------------------------------------------------
C     INITIALIZE APPROXIMATION TO CRITICAL BED SHEAR STRESS
C     USING KOMAR APPROXIMATION  EQUATIONS
C----------------------------------------------------------------------
            IF (DIMID(I,J) .LE. 0.0008D0) THEN
               TCI = 0.1D0*17.7D0*(DIMID(I,J)*100.0D0)**0.59D0
            ELSE
               TCI = 0.1D0*65.62D0*(DIMID(I,J)*100.0D0)**1.11D0
            ENDIF
            REG = SQRT(TCI/RHO)*DIMID(I,J)/VISKIN
            IF (REG .LE. 1.0D0) THEN
C----------------------------------------------------------------------
C     VALUE FOR HYDRAULICALLY SMOOTH BED
C----------------------------------------------------------------------
               TC1(I,J) = 0.1D0*(SIGMA(J)-RHO)*G*RHO**0.15D0*VISKIN**
     #            0.3D0
               TC1(I,J) = (TC1(I,J)*DIMID(I,J)**0.7D0)**0.869565
            ELSE IF (REG .LE. 60.0) THEN
C----------------------------------------------------------------------
C     VALUE FOR HYDRAULICALLY TRANSITIONAL BED
C----------------------------------------------------------------------
    1          CONTINUE
               THETLN = DLOG(TCI/((SIGMA(J)-RHO)*G*DIMID(I,J)))
               RELN = DLOG(SQRT(TCI/RHO)*DIMID(I,J)/VISKIN)
               FTC=(-2.26D0)-0.905D0*RELN+0.168D0*RELN*RELN-THETLN
               DFTC = (-1.45D0+0.168D0*RELN)/TCI
               TCIN = TCI - FTC/DFTC
               IF (ABS(TCIN-TCI) - 0.001D0 .LE. 0) GO TO 2
               TCI = TCIN
               GO TO 1
    2          CONTINUE
               TC1(I,J) = TCIN
            ELSE
C----------------------------------------------------------------------
C     VALUE FOR HYDRAULICALLY ROUGH BED
C----------------------------------------------------------------------
               TC1(I,J) = THETAC*G*(SIGMA(J)-RHO)*DIMID(I,J)
            ENDIF
C----------------------------------------------------------------------
C     VALUE OF CRITICAL SHEAR VELOCITY
C----------------------------------------------------------------------
            DO 3 JN = 2, M
               TC(I,J,JN) = TC1(I,J)
               SHVELC(I,J,JN) = SQRT(TC1(I,J)/RHO)
    3       CONTINUE
    4    CONTINUE
    5 CONTINUE
      RETURN 
      END
