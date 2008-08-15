C**********************SUBROUTINE DRANGE*******************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES THE GRAIN SIZE FRACTIONS THAT CAN BE
C                   TRANSPORTED AS BEDLOAD FOR EACH SHEAR STRESS RANGE
C  CREATORS: R. SLINGERLAND & K. VOGEL
C----------------------------------------------------------------------
      SUBROUTINE DRANGE
      IMPLICIT NONE
      REAL*8 RHO,G                                
      REAL*8 VISKIN,THETAC                       
      REAL*8 DIMID(20,3),SIGMA(3)               
      REAL*8 DIFEND(20,3),DIMIN(3)             
      REAL*8 PROPT(50)                        
      REAL*8 SHVEL(50,100),TO(50,100)        
      REAL*8 DMAX(50,3,100),DMIN(50,3,100)  
      REAL*8 THET50(100),VISC              
      REAL*8 COVAR,A,B,CON1,CON2          
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100),
     #AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 AWPM(20,3,100),D50M(100)                                 
      INTEGER*4 NK
      INTEGER*4 IMIN(50,3,100),IMAX(50,3,100)
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*4 NSIZE(3),NSIGMA,M,N                                  
C----------------------------------------------------------------------
C     COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM2/VISKIN,THETAC
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM6/PROPT,NK
      COMMON /COM7/SHVEL,TO
      COMMON /COM8/DMAX,DMIN,IMIN,IMAX
      COMMON /COM12/THET50,VISC
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #OPTBINP,OPTSINP,OPTERO
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM22/AWPM,D50M
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C            DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 REG,AA,BB,CC,DIMSV,DLOG10,DIMDL,FD,DFD,DIN,DI,SETVEL
      INTEGER*4 I,K,J,JN
C----------------------------------------------------------------------
C     CALCULATE THE GRAIN DIAMETER AT ENTRAINMENT THRESHOLD (DMAX)
C     FOR A GIVEN GRAIN DENSITY USING BRIDGE(1981) REGRESSION EQNS
C----------------------------------------------------------------------
      DO 9 JN = 2, M
         DO 2 K = 1, NK
            REG = SHVEL(K,JN)*D50M(JN)/VISKIN
            DO 1 J = 1, NSIGMA
               IF (REG .LE. 1.0D0) THEN
C----------------------------------------------------------------------
C     HYDRAULICALLY SMOOTH BED
C----------------------------------------------------------------------
                  DMAX(K,J,JN) = (TO(K,JN)**1.15D0/(0.1D0*(SIGMA(J)-RHO)
     #               *G*VISKIN**0.3D0*RHO**0.15D0))**1.43D0
               ELSE IF (REG .LE. 60.0D0) THEN
C---------------------------------------------------------------------
C     HYDRAULICALLY TRANSITIONAL BED: SOLVE EQUATION FOR DMAX USING
C     QUADRATIC EQUATION
C---------------------------------------------------------------------
                  AA = 0.16841D0
                  BB = AA*2.0D0*DLOG(SHVEL(K,JN)/VISKIN) + 0.09512D0
                  CC = (-DLOG(TO(K,JN)/(G*(SIGMA(J)-RHO)))) - 2.26D0 - 
     #               0.90488D0*DLOG(SHVEL(K,JN)/VISKIN) + AA*DLOG(SHVEL(
     #               K,JN)/VISKIN)**2
                  IF (BB**2 - 4.0D0*AA*CC .LE. 0.0D0) THEN
                     DMAX(K,J,JN) = (TO(K,JN)**1.15D0/(0.1D0*(SIGMA(J)-
     #                  RHO)*G*VISKIN**0.3D0*RHO**0.15D0))**1.43D0
                  ELSE
                     DMAX(K,J,JN) = DEXP(((-BB)+DSQRT(BB**2-4.0D0*AA*CC)
     #                  )/(2.0D0*AA))
                  ENDIF
               ELSE
C---------------------------------------------------------------------
C     HYDRAULICALLY ROUGH BED
C---------------------------------------------------------------------
                  DMAX(K,J,JN) = TO(K,JN)/(THETAC*(SIGMA(J)-RHO)*G)
               ENDIF
    1       CONTINUE
    2    CONTINUE
C----------------------------------------------------------------------
C     CALCULATE DMIN, THE GRAIN DIAMETER AT THE SUSPENSION THRESHOLD
C----------------------------------------------------------------------
         DO 5 K = 1, NK
            SETVEL = B*SHVEL(K,JN)
            DO 4 J = 1, NSIGMA
               DIMSV=RHO*SETVEL*SETVEL*SETVEL/((SIGMA(J)-RHO)*G*VISKIN)
C----------------------------------------------------------------------
C     FIRST APPROXIMATION TO GRAIN SIZE, DI
C----------------------------------------------------------------------
               IF (DIMSV .LE. 1.0D0) THEN
                  DI = DSQRT(18.0D0*SETVEL*VISC/((SIGMA(J)-RHO)*G))
               ELSE
                  DI = 0.3D0*SETVEL*SETVEL*RHO/((SIGMA(J)-RHO)*G)
               ENDIF
C----------------------------------------------------------------------
C     SOLVE DIETRICH EQUATION FOR GRAIN SIZE USING NEWTON-RAPHSON
C----------------------------------------------------------------------
               DIMDL = DLOG10(DI*DI*DI*G*(SIGMA(J)-RHO)/(RHO*VISKIN*
     #            VISKIN))
    3          CONTINUE
               FD = (-3.767D0) + 1.929D0*DIMDL - 0.0981D0*DIMDL**2 - 
     #            0.00575D0*DIMDL**3 + 0.00056D0*DIMDL**4 - DLOG10(DIMSV
     #            )
               DFD = 1.929D0 - 0.196D0*DIMDL - 0.01725D0*DIMDL**2 + 
     #            0.00224D0*DIMDL**3
               DIN = DIMDL - FD/DFD
               IF (DABS(DIN-DIMDL) .GT. 0.000001D0) THEN
                  DIMDL = DIN
                  GO TO 3
               ELSE
                  DIMDL = 10.0D0**DIMDL
                  DMIN(K,J,JN) = (DIMDL*RHO*VISKIN**2/((SIGMA(J)-RHO)*G)
     #               )**0.333D0
               ENDIF
    4       CONTINUE
    5    CONTINUE
C----------------------------------------------------------------------
C     DEFINE GRAIN SIZE ARRAY INDICES FOR DMIN AND DMAX
C----------------------------------------------------------------------
         DO 8 K = 1, NK
            DO 7 J = 1, NSIGMA
               DO 6 I = 2, NSIZE(J)
                  IF (DMIN(K,J,JN) .GE. DIFEND(I-1,J)) THEN
                     IF (DMIN(K,J,JN) .LT. DIMID(I,J)) THEN
                        IMIN(K,J,JN) = I
                     ELSE
                        IMIN(K,J,JN) = I + 1
                     ENDIF
                  ENDIF
                  IF (DMAX(K,J,JN) .GE. DIFEND(I,J)) IMAX(K,J,JN) = I
    6          CONTINUE
               IF(DMAX(K,J,JN).LT.DIFEND(1,J))IMAX(K,J,JN)=NSIZE(J)+100
               IF (DMIN(K,J,JN) .LT. DIFEND(1,J)) IMIN(K,J,JN) = 1
    7       CONTINUE
    8    CONTINUE
    9 CONTINUE
      RETURN 
      END
