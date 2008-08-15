C*************************SUBROUTINE DRANGEH**************************
C---------------------------------------------------------------------
C THIS SUBROUTINE: CALCULATES THE GRAIN SIZE FRACTIONS THAT CAN BE
C                  TRANSPORTED AS BEDLOAD FOR EACH SHEAR STRESS RANGE
C---------------------------------------------------------------------
      SUBROUTINE DRANGEH
      IMPLICIT NONE
      REAL*8 RHO,G
      REAL*8 VISKIN,THETAC                                              
      REAL*8 DIMID(20,3),SIGMA(3)                                      
      REAL*8 DIFEND(20,3),DIMIN(3)                                    
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)                
      REAL*8 PROPT(50)                                              
      INTEGER*4 NK
      REAL*8 SHVEL(50,100),TO(50,100)                              
      REAL*8 DMAX(50,3,100),DMIN(50,3,100)                        
      INTEGER*4 IMIN(50,3,100),IMAX(50,3,100)
      REAL*8 THET50(100),VISC                                    
      REAL*8 COVAR,A,B,CON1,CON2                                
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      REAL*8 AWPM(20,3,100),D50M(100)                          
      INTEGER*4 NSIZE(3),NSIGMA,M,N                           
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM2/VISKIN,THETAC
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM6/PROPT,NK
      COMMON /COM7/SHVEL,TO
      COMMON /COM8/DMAX,DMIN,IMIN,IMAX
      COMMON /COM12/THET50,VISC
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #OPTBINP,OPTSINP,OPTERO
      COMMON /COM22/AWPM,D50M
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C                  DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 DI,SUR,X1,X2,FMID,F,RTBIS,DX,SETVEL,DIMSV,DIMDL,FD
      REAL*8 DFD,DIN,XACC,FUN1,D50
      INTEGER*4 J,K,I,JN
      FUN1(DI,D50,SUR,CON2)=DLOG10(CON2*DI/D50)-DSQRT(DI/SUR)
C----------------------------------------------------------------------
      DO 19 JN = 2, M
         D50 = D50M(JN)
C----------------------------------------------------------------------
C     DETERMINE THE HIDING FUNCTION TO BE USED
C----------------------------------------------------------------------
         GO TO (1,4) OPTFUNC
C----------------------------------------------------------------------
C     KOMAR EQUATION 6 HIDING FUNCTION
C----------------------------------------------------------------------
    1    CONTINUE
         DO 3 J = 1, NSIGMA
            DO 2 K = 1, NK
               DMAX(K,J,JN) = (TO(K,JN)*D50**CON1/(THET50(JN)*(SIGMA(J)-
     #            RHO)*G))**(1/(1.+CON1))
    2       CONTINUE
    3    CONTINUE
         GO TO 10
C----------------------------------------------------------------------
C     EGIAZAROFF HIDING FUNCTION
C----------------------------------------------------------------------
    4    CONTINUE
         XACC = 0.0000001D0
         X1 = 2.7183**2*D50/CON2
         X2 = 80.0D0*D50
         DO 9 J = 1, NSIGMA
            DO 8 K = 1, NK
               SUR = TO(K,JN)/((SIGMA(J)-RHO)*G*DLOG10(CON2)**2.*THET50(
     #            JN))
C----------------------------------------------------------------------
C     USE BISECTION TO FIND THE ROOT OF THE EGIAZAROFF FUNCTION
C----------------------------------------------------------------------
               FMID = REAL(DLOG10(CON2*X2/D50)-DSQRT(X2/SUR))
               F = REAL(DLOG10(CON2*X1/D50)-DSQRT(X1/SUR))
               IF (F*FMID .GE. 0.0D0) THEN
C----------------------------------------------------------------------
C     IF NO ROOT EXISTS FOR THIS TO, THEN SET DMAX TO A VALUE THAT
C     WILL MAKE IMAX = NSIZE(J)+100
C----------------------------------------------------------------------
                  DMAX(K,J,JN) = DIFEND(1,J) - 1.0D0
                  GO TO 7
               ENDIF
               IF (F .LT. 0.0D0) THEN
                  RTBIS = X1
                  DX = X2 - X1
               ELSE
                  RTBIS = X2
                  DX = X1 - X2
               ENDIF
    5          CONTINUE
               DX = DX*0.5D0
               DI = RTBIS + DX
               FMID = REAL(DLOG10(CON2*DI/D50)-DSQRT(DI/SUR))
               IF (FMID .LT. 0.0D0) RTBIS = DI
               IF (DABS(DX).LT.XACC .OR. FMID.EQ.0.0D0) GO TO 6
               GO TO 5
    6          CONTINUE
               DMAX(K,J,JN) = DI
    7          CONTINUE
    8       CONTINUE
    9    CONTINUE
   10    CONTINUE
C----------------------------------------------------------------------
C     CALCULATE THE GRAIN DIAMETER AT THE SUSPENSION THRESHOLD
C----------------------------------------------------------------------
         DO 15 K = 1, NK
            SETVEL = B*SHVEL(K,JN)
            IF (SHVEL(K,JN) .LT. 1.0E-6) THEN
               DO 11 J = 1, NSIGMA
                  DMIN(K,J,JN) = 0.0D0
                  DMAX(K,J,JN) = 0.0D0
                  IMIN(K,J,JN) = 1
                  IMAX(K,J,JN) = NSIZE(J) + 100
   11          CONTINUE
               GO TO 14
            ENDIF
            DO 13 J = 1, NSIGMA
               DIMSV=RHO*SETVEL*SETVEL*SETVEL/((SIGMA(J)-RHO)*G*VISKIN)
C----------------------------------------------------------------------
C     FIRST APPROXIMATION
C----------------------------------------------------------------------
               IF (DIMSV .LE. 1.0D0) THEN
                  DI = DSQRT(18.0D0*SETVEL*VISC/(SIGMA(J)-RHO)*G)
               ELSE
                  DI = 0.3D0*SETVEL*SETVEL*RHO/((SIGMA(J)-RHO)*G)
               ENDIF
C----------------------------------------------------------------------
C     SOLVE DIETRICH EQUATION FOR GRAIN SIZE USING NEWTON-RAPHSON
C----------------------------------------------------------------------
               DIMDL = DLOG10(DI*DI*DI*G*(SIGMA(J)-RHO)/(RHO*VISKIN*
     #            VISKIN))
   12          CONTINUE
               FD = (-3.767D0) + 1.929D0*DIMDL - 0.098D0*DIMDL**2 - 
     #            0.00575D0*DIMDL**3 + 0.00056D0*DIMDL**4 - DLOG10(DIMSV
     #            )
               DFD = 1.929D0 - 0.196D0*DIMDL - 0.01725D0*DIMDL**2 + 
     #            0.00224D0*DIMDL**3
               DIN = DIMDL - FD/DFD
               IF (DABS(DIN-DIMDL) .GE. 0.000001D0) THEN
                  DIMDL = DIN
                  GO TO 12
               ELSE
                  DIMDL = 10.0D0**DIMDL
                  DMIN(K,J,JN) = (DIMDL*RHO*VISKIN**2/((SIGMA(J)-RHO)*G)
     #               )**0.333D0
               ENDIF
   13       CONTINUE
   14       CONTINUE
   15    CONTINUE
C----------------------------------------------------------------------
C     DEFINE GRAIN SIZE ARRAY INDICES FOR DMIN AND DMAX
C----------------------------------------------------------------------
         DO 18 K = 1, NK
            DO 17 J = 1, NSIGMA
               DO 16 I = 2, NSIZE(J)
                  IF (DMIN(K,J,JN) .GE. DIFEND(I-1,J)) THEN
                     IF (DMIN(K,J,JN) .LT. DIMID(I,J)) THEN
                        IMIN(K,J,JN) = I
                     ELSE
                        IMIN(K,J,JN) = I + 1
                     ENDIF
                  ENDIF
                  IF (DMAX(K,J,JN) .GE. DIFEND(I,J)) IMAX(K,J,JN) = I
   16          CONTINUE
               IF(DMAX(K,J,JN).LT.DIFEND(1,J))IMAX(K,J,JN)=NSIZE(J)+100
               IF (DMIN(K,J,JN) .LT. DIFEND(1,J)) IMIN(K,J,JN) = 1
   17       CONTINUE
   18    CONTINUE
   19 CONTINUE
      RETURN 
      END
