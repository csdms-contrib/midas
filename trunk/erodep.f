C**********************  SUBROUTINE ERODEP  *************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES THE TOTAL EROSION OR DEPOSITION
C                   OCCURRING IN A TIME STEP AND CALCULATES THE NEW
C                   PARTICLE SIZE DISTRIBUTIONS IN ACTIVE LAYER
C----------------------------------------------------------------------
      SUBROUTINE ERODEP
      IMPLICIT NONE
      REAL*8 DIMID(20,3),SIGMA(3)    
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100), 
     #       AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 AWPM(20,3,100),D50M(100)                                  
      REAL*8 WIDM(100), ELEVM(100), DCHM(100), TOXM(100), SVELM(100), 
     #       SVELUM(100)
      INTEGER*4 NSIZE(3),NSIGMA,M,N                                  
      REAL*8 ELEVOM(100)                                            
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #              AVDENS
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM32/ELEVOM
C----------------------------------------------------------------------
C            DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 MSERO,MSDEP,MSACTIVE,MSTHK,MSEROACT
      REAL*8 TOTNEWMS,MSTOT,MSDEPACT
      REAL*8 NEWMSACT(20,3),NEWACTIVE
      INTEGER*4 I,J,JN
C----------------------------------------------------------------------
      DO 19 JN = 2, M
         MSACTIVE = ACTIVE*VOLCON*NWAVDENS(JN)
         IF (ELEVOM(JN) .GE. ELEVM(JN)) THEN
            MSEROACT = ACTIVE*VOLCON*AVDENS(JN)
         ELSE
            MSEROACT = MSACTIVE
         ENDIF
         MSDEPACT = MSACTIVE
C----------------------------------------------------------------------
C       MSPRTHK=MASS ADDED OR REMOVED OF EACH SIZE AND DENSITY FRACTION
C        UNITS= KG/M2
C----------------------------------------------------------------------
         MSERO = 0.0D0
         MSDEP = 0.0D0
C----------------------------------------------------------------------
C     CALCULATE EROSION/DEPOSITION AND NEW BED COMPOSITION;
C     MSERO=SUM OF MASS/AREA ERODED; MSDEP=SUM OF MASS/AREA DEPOSITED;
C     MSPRTHK=MASS/AREA ERODED FOR MINERAL FRACTION IJ( >0 FOR EROSION)
C     MSACTIVE=TOTAL MASS/AREA IN THE ACTIVE LAYER
C----------------------------------------------------------------------
         DO 2 J = 1, NSIGMA
            DO 1 I = 1, NSIZE(J)
               PRTHICK(I,J,JN) = MSPRTHK(I,J,JN)/(SIGMA(J)*VOLCON)
               IF (MSPRTHK(I,J,JN) .GT. 0.0D0) THEN
                  MSERO = MSERO + MSPRTHK(I,J,JN)
               ELSE
                  MSDEP = MSDEP + MSPRTHK(I,J,JN)
               ENDIF
    1       CONTINUE
    2    CONTINUE
         MSTOT = MSDEP + MSERO
         THICKM(JN) = 0.0D0
         DO 4 J = 1, NSIGMA
            DO 3 I = 1, NSIZE(J)
               THICKM(JN) = THICKM(JN) + PRTHICK(I,J,JN)
    3       CONTINUE
    4    CONTINUE
C----------------------------------------------------------------------
C      REPLACE ERODED MATERIAL AT BASE OF ACTIVE LAYER
C      MSTHK=TOTAL MASS/AREA ADDED OR SUBTRACTED FROM THE ACTIVE LAYER
C----------------------------------------------------------------------
       MSTHK=MSACTIVE-TOTNEWMS
C----------------------------------------------------------------------
C      IF THICKM(JN).LT.0 THEN THERE IS DEPOSITION,
C----------------------------------------------------------------------
         IF (THICKM(JN) .LT. 0.0D0) THEN
C----------------------------------------------------------------------
C  IF A THICKNESS GREATER THAN THE ACTIVE LAYER HAS BEEN DEPOSITED THEN
C  THE NEW ACTIVE LAYER TAKES ON THE SIZE DISTRIBUTION OF THE DEPOSITED
C  MATERIAL.
C----------------------------------------------------------------------
            NEWACTIVE = ACTIVE
            TOTNEWMS = 0.0D0
C----------------------------------------------------------------------
C  CALCULATE THE NEW LAYER THICKNESS AND MASS AFTER THE ERODIBLE
C  FRACTIONS HAVE BEEN REMOVED.
C----------------------------------------------------------------------
            IF ((-THICKM(JN)) .GE. ACTIVE) THEN
               AWPM(I,J,JN) = MSPRTHK(I,J,JN)/MSDEP
               IF (AWPM(I,J,JN) .LE. 0.0D0) AWPM(I,J,JN) = 0.0D0
            ELSE
C----------------------------------------------------------------------
C  CALCULATE THE NEW LAYER THICKNESS AND MASS AFTER THE ERODIBLE
C  FRACTIONS HAVE BEEN REMOVED.
C----------------------------------------------------------------------
               DO 6 J = 1, NSIGMA
                  DO 5 I = 1, NSIZE(J)
C----------------------------------------------------------------------
C            CALCULATE THE NEW MASS/AREA IN THE NEW ACTIVE LAYER
C            NEWMSACT(I,J)=KG/M2
C            TOTNEWMS=TOTAL NEW MASS/AREA ADDED OR REMOVED AT A NODE
C            TOTNEWMS=KG/M2
C            THICKM(JN)=TOTAL THICKNESS ERODED OR DEP. AT NODE (M)
C----------------------------------------------------------------------
                     IF (MSPRTHK(I,J,JN) .GT. 0.0D0) NEWACTIVE = 
     #                  NEWACTIVE - PRTHICK(I,J,JN)
    5             CONTINUE
    6          CONTINUE
               DO 8 J = 1, NSIGMA
                  DO 7 I = 1, NSIZE(J)
                     NEWMSACT(I,J) = AWPM(I,J,JN)*MSACTIVE - MSPRTHK(I,J
     #                  ,JN)
                     IF (NEWMSACT(I,J) .LT. 0.0) NEWMSACT(I,J) = 0.0
                     TOTNEWMS = TOTNEWMS + NEWMSACT(I,J)
    7             CONTINUE
    8          CONTINUE
               DO 10 J = 1, NSIGMA
                  DO 9 I = 1, NSIZE(J)
                     AWPM(I,J,JN) = NEWMSACT(I,J)/TOTNEWMS
    9             CONTINUE
   10          CONTINUE
            ENDIF
         ELSE
C----------------------------------------------------------------------
C         IF THICKM(JN).GT.0) THEN THERE IS EROSION
C----------------------------------------------------------------------
            TOTNEWMS = MSACTIVE
C----------------------------------------------------------------------
C      CALCULATE MASS PER UNIT AREA OF EACH GRAIN FRACTION REMAINING IN
C      ACTIVE LAYER FOLLOWING EROSION IN THIS TIME INCREMENT
C      ACTERO=REVISED MASS/AREA LEFT IN THE ACTIVE LAYER
C      THICK=THICKNESS ERODED FROM THE ACTIVE LAYER( >0 FOR EROSION )
C      NEWMSACT(I,J)=REVISED MASS/AREA OF IJ IN THE ACTIVE LAYER
C----------------------------------------------------------------------
            DO 12 J = 1, NSIGMA
               DO 11 I = 1, NSIZE(J)
C----------------------------------------------------------------------
C            CALCULATE THE NEW MASS/AREA IN THE NEW ACTIVE LAYER
C            NEWMSACT(I,J)=KG/M2
C            TOTNEWMS=TOTAL NEW MASS/AREA ADDED OR REMOVED AT A NODE
C            TOTNEWMS=KG/M2
C            THICKM(JN)=TOTAL THICKNESS ERODED OR DEP. AT NODE (M)
C----------------------------------------------------------------------
                  NEWMSACT(I,J)=AWPM(I,J,JN)*MSACTIVE-MSPRTHK(I,J,JN)
                  IF (NEWMSACT(I,J) .LT. 0.0) NEWMSACT(I,J) = 0.0
                  TOTNEWMS = TOTNEWMS - MSPRTHK(I,J,JN)
   11          CONTINUE
   12       CONTINUE
            MSACTIVE = 0.0D0
            IF (ELEVOM(JN) .GE. ELEVM(JN)) THEN
               DO 14 J = 1, NSIGMA
                  DO 13 I = 1, NSIZE(J)
                     NEWMSACT(I,J) = NEWMSACT(I,J) + AWPOM(I,J,JN)*
     #                  MSEROACT*THICKM(JN)/ACTIVE
                     MSACTIVE = MSACTIVE + NEWMSACT(I,J)
   13             CONTINUE
   14          CONTINUE
C----------------------------------------------------------------------
C RESET THE ELEVOM VALUE TO THE NEW BED ELEVATION
C----------------------------------------------------------------------
               ELEVOM(JN) = ELEVM(JN) - THICKM(JN)
            ELSE
               DO 16 J = 1, NSIGMA
                  DO 15 I = 1, NSIZE(J)
                     NEWMSACT(I,J) = NEWMSACT(I,J) + AWPM(I,J,JN)*
     #                  MSEROACT*THICKM(JN)/ACTIVE
                     MSACTIVE = MSACTIVE + NEWMSACT(I,J)
   15             CONTINUE
   16          CONTINUE
            ENDIF
            DO 18 J = 1, NSIGMA
               DO 17 I = 1, NSIZE(J)
                  AWPM(I,J,JN) = NEWMSACT(I,J)/MSACTIVE
   17          CONTINUE
   18       CONTINUE
         ENDIF
C----------------------------------------------------------------------
C       UPDATE THE ELEVATION OF THE STREAM PROFILE
C----------------------------------------------------------------------
         ELEVM(JN) = ELEVM(JN) - THICKM(JN)
   19 CONTINUE
      ELEVM(M+1) = 2.*ELEVM(M) - ELEVM(M-1)
      RETURN 
      END
C***************************SUBROUTINE EXVERT**************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: COMPUTES VERTICAL DIFFUSIVITIES AT THE MIDPOINTS
C                   OF VERTICAL INTERVALS.
C----------------------------------------------------------------------
      SUBROUTINE EXVERT
      IMPLICIT NONE
      REAL*8 YMAX(100),DELY(100),VEL(30,100),EPSY(30,100),PSI,BETA 
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),
     #SVELUM(100)                                                 
      INTEGER*4 NSIZE(3),NSIGMA,M,N                              
C----------------------------------------------------------------------
C                COMMON DECLARATION
C----------------------------------------------------------------------
      COMMON /COM19/YMAX,DELY,VEL,EPSY,PSI,BETA
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C               DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 Y
      INTEGER*4 L,JN
C----------------------------------------------------------------------
      DO 2 JN = 2, M
         DO 1 L = 1, N - 1
            Y = DCHM(JN) - DELY(JN)/2 - (L-1)*DELY(JN)
            EPSY(L,JN) = 0.4*SVELUM(JN)*Y*(1.0-Y/DCHM(JN))
    1    CONTINUE
    2 CONTINUE
      RETURN 
      END
