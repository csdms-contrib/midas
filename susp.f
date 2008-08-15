C**************************SUBROUTINE SUSP*****************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES SUSPENDED LOAD TRANSPORT RATES
C----------------------------------------------------------------------
      SUBROUTINE SUSP()
      IMPLICIT NONE
      REAL*8 DIMID(20,3),SIGMA(3)                               
      REAL*8 DIFEND(20,3),DIMIN(3)                             
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)       
      REAL*8 SLWTM(20,3,100),SLWTAM(20,3,100)              
      REAL*8 BEDCON(20,3,100)                             
      REAL*8 COVAR,A,B,CON1,CON2                         
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100),
     #AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 CONIM(20,3,100),CVELM(20,3,100),PHI                      
      REAL*8 DBDX(20,3,100),XFLUX(20,3,100),TXFLUX(100),TDBDX(100)   
      REAL*8 CONCM(20,3,30),CONC(30)                                
      REAL*8 YMAX(100),DELY(100),VEL(30,100),EPSY(30,100),PSI,BETA 
      REAL*8 AWPM(20,3,100),D50M(100)                             
      REAL*8 BLWTM(20,3,100),BLWTAM(20,3,100),SUM1(100),SUM2(100)
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100), 
     #SVELUM(100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                             
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*2 OPTSUSP,OPTLNODE
      INTEGER*4 NSIZE(3),NSIGMA,M,N                             
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM9/SLWTM,SLWTAM
      COMMON /COM10/BEDCON
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #OPTBINP,OPTSINP,OPTERO
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM15/CONIM,CVELM,PHI,OPTSUSP,OPTLNODE
      COMMON /COM16/DBDX,XFLUX,TDBDX,TXFLUX
      COMMON /COM18/CONCM,CONC
      COMMON /COM19/YMAX,DELY,VEL,EPSY,PSI,BETA
      COMMON /COM22/AWPM,D50M
      COMMON /COM23/BLWTAM,BLWTM,SUM1,SUM2
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--EXTERNAL VARIABLES NOT IN COMMON BLOCKS
C----------------------------------------------------------------------
C     INTEGER*4 NT
C----------------------------------------------------------------------
      INTEGER*4 I,IN,JN,J
      REAL*8 B1,BCON,BWT,BWTA,SWT,SWTA,W
      REAL*8 MSACT(100)
      REAL*8 AV,Y,CA,Z,DY1,PHI1
      REAL*8 LOWLIM,HILIM,FACTOR
      INTEGER*4 IV,L,SIGNMSP,NEWSIGN,REDUC
C----------------------------------------------------------------------
      DO 14 JN = 2, M + 1
         NWAVDENS(JN) = 0.0D0
         MSACT(JN) = 0.0D0
         DO 2 J = 1, NSIGMA
            DO 1 I = 1, NSIZE(J)
               NWAVDENS(JN) = NWAVDENS(JN) + AWPM(I,J,JN)*SIGMA(J)
    1       CONTINUE
    2    CONTINUE
         MSACT(JN) = ACTIVE*VOLCON*NWAVDENS(JN)
         DO 13 J = 1, NSIGMA
            DO 12 I = 1, NSIZE(J)
C----------------------------------------------------------------------
C  CALCULATE PROFILES FOR SUSPENDED PARTICLES
C      IN JN=M+1 WE ARE AT THE LAST NODE AND CVELM & BLWTM ARE
C      SET AS BOUNDARY CONDITIONS
C-------------------------------------------------------------------
               REDUC = 0
    3          CONTINUE
               CA = BEDCON(I,J,JN)
               IF (JN .NE. M+1) THEN
                  W = SETV(I,J)
C-------------------------------------------------------------------
C        IF NO SUSPENDIBLE MATERIAL AVAILABLE SET SUSP LOAD=0.0
C-------------------------------------------------------------------
                  IF (CA.EQ.0.0D0 .AND. CONIM(I,J,JN-1).EQ.0.0D0) THEN
                     DO 4 IV = 1, N
                        CONC(IV) = 0.0D0
    4                CONTINUE
                     CONIM(I,J,JN) = 0.0D0
                     CVELM(I,J,JN) = 0.0D0
                     GO TO 8
                  ELSE
C---------------------------------------------------------------------
C      CALCULATE THE CONCENTRATIONS AT THE GRID POINTS
C      IF OPTSUSP=1 THEN ANDRE VAN NIEKERK'S ALGORITHM
C      INCORPORATING LONGITUDINALLY ADVECTIVE SUSPENDED
C      LOAD FLUX IS USED
C---------------------------------------------------------------------
                     IF (OPTSUSP .EQ. 1) THEN
                        DY1 = 1.0
C   5                   CONTINUE
                     ELSE
C----------------------------------------------------------------------
C     CALCULATE EQUILIBRIUM CONCENTRATIONS USING THE ROUSE EQUATIONS
C----------------------------------------------------------------------
                        Z = W/(0.4*SVELUM(JN))
                        AV = YMAX(JN)
                        DO 6 L = 1, N - 1
                           Y = DCHM(JN) - (L-1)*DELY(JN)
                           IF (L .EQ. 1) Y = DCHM(JN) - DELY(JN)/4
                           CONC(L) = CA*(((DCHM(JN)-Y)*AV)/(Y*(DCHM(JN)-
     #                        AV)))**Z
    6                   CONTINUE
                        CONC(N) = CA
                     ENDIF
                  ENDIF
C----------------------------------------------------------------------
C  NOW THAT CONCENTRATIONS OF SUSPENDED MATERIAL AT VARIOUS
C  STREAM LEVELS ARE KNOWN (CONC(IV)), INTEGRATE THESE
C  TO GET THE SUSPENDED LOAD FLUX.
C  CALCULATE THE CONCENTRATION PROFILE INTEGRALS
C    BY CALCULATING THE CONCENTRATION INTEGRALS AT EACH GRID
C----------------------------------------------------------------------
                  CONIM(I,J,JN) = 0.
                  CVELM(I,J,JN) = 0.
                  DO 7 IV = 1, N - 1
                     CONIM(I,J,JN) = CONIM(I,J,JN) + (CONC(IV)+CONC(IV+1
     #                  ))*(DELY(JN)/2.)
                     CVELM(I,J,JN) = CVELM(I,J,JN) + (CONC(IV)*VEL(IV,JN
     #                  )+CONC(IV+1)*VEL(IV+1,JN))*(DELY(JN)/2.)
    7             CONTINUE
                  CONIM(I,J,JN) = CONIM(I,J,JN) + CONC(N)*YMAX(JN)/2.
                  IF (REDUC .EQ. 0) BLWTM(I,J,JN) = BLWTM(I,J,JN) + (
     #               CONC(N)*VEL(N,JN))*YMAX(JN)/2.
               ENDIF
    8          CONTINUE
               IF (JN .EQ. 2) THEN
                  IF (OPTBINP .EQ. 2) BLWTM(I,J,1) = BLWTM(I,J,2)
                  IF (OPTSINP .EQ. 2) CVELM(I,J,1) = CVELM(I,J,2)
               ENDIF
C----------------------------------------------------------------------
C  CONTINUITY PORTION OF MASSFLUX CALCS
C  CALL CONTINUITY EQUATION CALCULATION TO DETERMINE
C  NET MASS ADDED OR REMOVED AT EACH NODE
C  DETERMINE THE COEFFICIENTS IN THE FINITE DIFFERENCE FORMULATION OF
C  THE SEDIMENT CONTINUITY EQUATION. NOTE: A WEIGHTING FACTOR OF PHI
C  IS USED.
C            [XFLUX]=KG/M, [CVEL*]=KG/MS
C            START BY CHECKING IF AT FIRST NODE
C----------------------------------------------------------------------
               IN = JN - 1
C--------------------------------------------------------------------
C     IF AT LAST NODE: SELECT BOUNDARY CONDITION USING OPTLNODE
C--------------------------------------------------------------------
               IF (IN .EQ. M) THEN
                  PHI1 = PHI
                  FACTOR = 1
C--------------------------------------------------------------------
C     IF OPTELNODE=0: CONDITION OF RIGID, UNERODABLE BED AT LAST
C     NODE: FLUX IN AT M = FLUX OUT AT M
C--------------------------------------------------------------------
                  IF (REDUC .EQ. 0) THEN
                     IF (OPTLNODE .EQ. 0) THEN
                        DBDX(I,J,IN) = 0.0
                        XFLUX(I,J,IN) = 0.0
                        GO TO 9
C--------------------------------------------------------------------
C      IF OPTLNODE = 1 THEN USE A BACKWARD DIFFERENCE USING CALCULATED
C      TRNASORT RATES, ON HALF THE DELX INTERVAL AT M
C--------------------------------------------------------------------
                     ELSE IF (OPTLNODE .EQ. 1) THEN
                        PHI1 = 1.
                        FACTOR = 2.
C--------------------------------------------------------------------
C      IF OPTLNODE = 2 DO A CENTRAL DIFFERENCE UNDER THE ASSUMPTION
C      THAT AT NODE M+1 THE TRANSPORT RATES ARE = TO THOSE AT M
C--------------------------------------------------------------------
                     ELSE IF (OPTLNODE .EQ. 2) THEN
                        BLWTM(I,J,M+1) = BLWTM(I,J,M)
                        CVELM(I,J,M+1) = CVELM(I,J,M)
                        PHI1 = 1.0
C--------------------------------------------------------------------
C      IF OPTLNODE=3 ASSUME TRANSPORT RATES AT M+1 ARE = 0
C--------------------------------------------------------------------
                     ELSE IF (OPTLNODE .EQ. 3) THEN
                        BLWTM(I,J,M+1)=2.*BLWTM(I,J,M)-BLWTM(I,J,M-1)
                        CVELM(I,J,M+1)=2.*CVELM(I,J,M)-CVELM(I,J,M-1)
                     ELSE IF (OPTLNODE .EQ. 4) THEN
                        DBDX(I,J,IN)=2.0*DBDX(I,J,IN-1)-DBDX(I,J,IN-2)
                        XFLUX(I,J,IN) = 2.0*XFLUX(I,J,IN-1) - XFLUX(I,J,
     #                     IN-2)
                        GO TO 9
                     ELSE IF (OPTLNODE .EQ. 5) THEN
                        DBDX(I,J,IN) = DBDX(I,J,IN-1)
                        XFLUX(I,J,IN) = XFLUX(I,J,IN-1)
                        GO TO 9
                     ENDIF
                     DBDX(I,J,IN) = ((1.-PHI1)*(BLWTM(I,J,IN+1)*WIDM(IN)
     #                  -BLWTM(I,J,IN)*WIDM(IN))+PHI1*FACTOR*(BLWTM(I,J,
     #                  IN)*WIDM(IN)-BLWTM(I,J,IN-1)*WIDM(IN-1)))*DELT/
     #                  DELX
                     XFLUX(I,J,IN) = ((1.-PHI1)*(CVELM(I,J,IN+1)*WIDM(IN
     #                  )-CVELM(I,J,IN)*WIDM(IN))+PHI1*FACTOR*(CVELM(I,J
     #                  ,IN)*WIDM(IN)-CVELM(I,J,IN-1)*WIDM(IN-1)))*DELT/
     #                  DELX
                  ELSE
                     MSPRTHK(I,J,IN) = AWPM(I,J,IN)*MSACT(IN)
                     PRTHICK(I,J,IN) = MSPRTHK(I,J,IN)/(SIGMA(J)*VOLCON)
                     GO TO 10
                  ENDIF
               ENDIF
C--------------------------------------------------------------------
C            IF NOT AT FIRST OR LAST NODE THEN
C--------------------------------------------------------------------
               IF (IN .EQ. 1) THEN
                  AWPM(I,J,1) = AWPM(I,J,2)
                  MSACT(1) = MSACT(2)
                  XFLUX(I,J,1) = (CVELM(I,J,IN+1)*WIDM(IN+1)-CVELM(I,J,
     #               IN)*WIDM(IN))*DELT/DELX
                  DBDX(I,J,1) = (BLWTM(I,J,2)*WIDM(2)-BLWTM(I,J,1)*WIDM(
     #               1))*DELT/DELX
               ELSE
                  DBDX(I,J,IN) = ((1.-PHI)*(BLWTM(I,J,IN+1)*WIDM(IN+1)-
     #               BLWTM(I,J,IN)*WIDM(IN))+PHI*(BLWTM(I,J,IN)*WIDM(IN)
     #               -BLWTM(I,J,IN-1)*WIDM(IN-1)))*DELT/DELX
                  XFLUX(I,J,IN) = ((1.-PHI)*(CVELM(I,J,IN+1)*WIDM(IN+1)-
     #               CVELM(I,J,IN)*WIDM(IN))+PHI*(CVELM(I,J,IN)*WIDM(IN)
     #               -CVELM(I,J,IN-1)*WIDM(IN-1)))*DELT/DELX
               ENDIF
    9          CONTINUE
               TXFLUX(IN) = TXFLUX(IN) + XFLUX(I,J,IN)
               TDBDX(IN) = TDBDX(IN) + DBDX(I,J,IN)
C-----------------------------------------------------------------------
C  CALCULATE AMOUNT AND SIZE DENSITY DISTRIBUTION OF MATERIAL TO BE
C  ERODED OR DEPOSITED
C-----------------------------------------------------------------------
               PRTHICK(I,J,IN) = (DBDX(I,J,IN)+XFLUX(I,J,IN))/(VOLCON*
     #            SIGMA(J)*WIDM(IN))
               MSPRTHK(I,J,IN) = (DBDX(I,J,IN)+XFLUX(I,J,IN))/WIDM(IN)
               IF (ABS(MSPRTHK(I,J,IN)) .LT. 1E-7) GO TO 10
               IF (MSPRTHK(I,J,IN).GT.AWPM(I,J,IN)*MSACT(IN) .OR. REDUC
     #            .EQ.1) THEN
                  IF (REDUC .EQ. 0) THEN
                     IF (MSPRTHK(I,J,IN) .NE. 0.0) THEN
                        SIGNMSP = MSPRTHK(I,J,IN)/ABS(MSPRTHK(I,J,IN))
                        NEWSIGN = SIGNMSP
                     ELSE
                        SIGNMSP = 1
                        NEWSIGN = 1
                     ENDIF
                  ELSE
                     IF (MSPRTHK(I,J,IN) .NE. 0.0) THEN
                        NEWSIGN = MSPRTHK(I,J,IN)/ABS(MSPRTHK(I,J,IN))
                     ELSE
                        NEWSIGN = 1
                     ENDIF
                  ENDIF
C-----------------------------------------------------------------------
C  IF AVAILABLE BED MATERIAL IS EXCEEDED THEN DO INTERVAL HALVING
C  TO DETERMINE THE ALLOWABLE TRANSPORT RATES AND MASS FLUXES
C  START BY HALFING THE FLUXES. REMEMBER BEDCON IS THE ACTIVE LAYER
C  THICKNESS.
C-----------------------------------------------------------------------
                  IF (REDUC .EQ. 0) THEN
                     B1 = BLWTM(I,J,IN+1) + SLWTM(I,J,IN+1)
                     BCON = BEDCON(I,J,IN+1)
                     BWT = BLWTM(I,J,IN+1)
                     BWTA = BLWTAM(I,J,IN+1)
                     SWT = SLWTM(I,J,IN+1)
                     SWTA = SLWTAM(I,J,IN+1)
                     LOWLIM = 0.
                     HILIM = 1.0
                     FACTOR = (HILIM-LOWLIM)/2.
                     BEDCON(I,J,IN+1) = FACTOR*BCON
                     BLWTM(I,J,IN+1) = FACTOR*BWT
                     BLWTAM(I,J,IN+1) = FACTOR*BWTA
                     SLWTM(I,J,IN+1) = FACTOR*SWT
                     SLWTAM(I,J,IN+1) = FACTOR*SWTA
C----------------------------------------------------------------------
C     RECALCULATE THE SUSPENDED LOAD TRANSPORT RATES BASED
C     ON THE NEW BEDLOAD TRANSPORT RATE
C     SET REDUC=1 TO LET SUSP KNOW WE ARE IN THE MASS REDUCTION
C     LOOP
C----------------------------------------------------------------------
                     REDUC = 1
                     GO TO 3
                  ENDIF
C----------------------------------------------------------------------
C         IF MSPRTHK HAS CHANGED SIGN DUE TO HALVING OR IF YOU ARE
C         APPROACHING FROM BELOW,BUT HAVE NOT APPROACHED THE TOLERANCE
C         THEN START HALVING THE DISTANCE UP TO PRTHICK.
C         SET THE LOW LIMIT TO THE PREVIOUS GUESS
C----------------------------------------------------------------------
                  IF (MSPRTHK(I,J,IN).NE.0.0 .AND. NEWSIGN.NE.SIGNMSP
     #                .OR. AWPM(I,J,IN)*MSACT(IN)-MSPRTHK(I,J,IN).GT.
     #               1.0E-9) THEN
                     LOWLIM = FACTOR
                     FACTOR = FACTOR + (HILIM-LOWLIM)/2.
                     BEDCON(I,J,IN+1) = FACTOR*BCON
                     BLWTM(I,J,IN+1) = FACTOR*BWT
                     BLWTAM(I,J,IN+1) = FACTOR*BWTA
                     SLWTM(I,J,IN+1) = FACTOR*SWT
                     SLWTAM(I,J,IN+1) = FACTOR*SWTA
                     GO TO 3
                  ELSE IF (ABS(HILIM) .LT. 1.0D-9) THEN
                     MSPRTHK(I,J,IN) = AWPM(I,J,IN)*MSACT(IN)
C----------------------------------------------------------------------
C            THE RIGHT SIDE OF THE IF CLAUSE IS:
C                  AWP*AVERAGE DENSITY OF THE ACTIVE LAYER/DENSITY OF
C                                                SIZE FRACTION
C                                                CONSIDERED
C----------------------------------------------------------------------
                  ELSE IF (MSPRTHK(I,J,IN) - AWPM(I,J,IN)*MSACT(IN)
     #                   .GT. 1.0E-9) THEN
C----------------------------------------------------------------------
C                  IF STILL TOO MUCH FLUX THEN START HALVING AGAIN
C                  UNTIL YOU BEGIN TO APPROACH PRTHICK FROM BELOW
C                  SET THE UPPER LIMIT TO THE PREVIOUS GUESS
C----------------------------------------------------------------------
                     HILIM = FACTOR
                     FACTOR = FACTOR - (HILIM-LOWLIM)/2.
                     BEDCON(I,J,IN+1) = FACTOR*BCON
                     BLWTM(I,J,IN+1) = FACTOR*BWT
                     BLWTAM(I,J,IN+1) = FACTOR*BWTA
                     SLWTM(I,J,IN+1) = FACTOR*SWT
                     SLWTAM(I,J,IN+1) = FACTOR*SWTA
                     GO TO 3
                  ELSE
                     MSPRTHK(I,J,JN) = MSACT(JN)*AWPM(I,J,JN)
                  ENDIF
               ENDIF
   10          CONTINUE
               DO 11 IV = 1, N
                  CONCM(I,J,IV) = CONC(IV)
   11          CONTINUE
   12       CONTINUE
   13    CONTINUE
   14 CONTINUE
      RETURN 
      END
