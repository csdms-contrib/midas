C*********************SUBROUTINE SETTLE********************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES THE CONSTANT TERMINAL SETTLING VELOCITY
C                   OF EACH SIZE-DENSITY FRACTION'S MEDIAN SIZE FROM
C                   DIETRICH'S EQUATION
C----------------------------------------------------------------------
      SUBROUTINE SETTLE
      IMPLICIT NONE
      REAL*8 RHO,G                                       
      REAL*8 VISKIN,THETAC                              
      REAL*8 DIMID(20,3),SIGMA(3)                      
      REAL*8 DIFEND(20,3),DIMIN(3)                    
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)
      INTEGER*4 NSIZE(3),NSIGMA,M,N                 
C----------------------------------------------------------------------
C      COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM2/VISKIN,THETAC
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C      DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 DIMDL,RHS
      INTEGER*4 I,J
C----------------------------------------------------------------------
      DO 2 J = 1, NSIGMA
         DO 1 I = 1, NSIZE(J)
            DIMDL = (SIGMA(J)-RHO)*G*DIMID(I,J)**3/(RHO*VISKIN*VISKIN)
            DIMDL = DLOG10(DIMDL)
            RHS = (-3.76715D0) + 1.92944D0*DIMDL - 0.09815D0*DIMDL**2
            RHS = RHS - 0.00575D0*DIMDL**3 + 0.00056D0*DIMDL**4
            SETV(I,J) = (((SIGMA(J)-RHO)*G*VISKIN*10.0D0**RHS)/RHO)**
     #         0.333
    1    CONTINUE
    2 CONTINUE
      RETURN 
      END
