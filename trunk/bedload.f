C***********************SUBROUTINE BEDLOAD*****************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: 1) CALCULATES THE BEDLOAD TRANSPORT RATES AND
C                      WEIGHTS PER UNIT AREA FOR EACH SIZE-DENSITY
C                      FRACTION
C
C                   2) NB: BEDLOAD TRANSPORT OF DIFFERENT SIZE-
C                          DENSITIES IS PROPORTIONED ACCORDING
C                          TO THE VOLUMES IN THE BED.
C  CREATORS: R. SLINGERLAND $ K. VOGEL
C----------------------------------------------------------------------
      SUBROUTINE BEDLOAD(NT)
      IMPLICIT NONE
      REAL*8 RHO,G                                         
      REAL*8 DIMID(20,3),SIGMA(3)                         
      REAL*8 DIFEND(20,3),DIMIN(3)                       
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)   
      REAL*8 PROPT(50)                                 
      REAL*8 SHVEL(50,100),TO(50,100)                 
      REAL*8 DMAX(50,3,100),DMIN(50,3,100)           
      REAL*8 SLWTM(20,3,100),SLWTAM(20,3,100)       
      REAL*8 BEDCON(20,3,100)                      
      REAL*8 COVAR,A,B,CON1,CON2                  
      REAL*8 YMAX(100),DELY(100),VEL(30,100),EPSY(30,100),PSI,BETA   
      REAL*8 AWPM(20,3,100),D50M(100)                               
      REAL*8 BLWTM(20,3,100),BLWTAM(20,3,100),SUM1(100),SUM2(100)  
      REAL*8 DELX,DELT,XPRINT,TPRINT                             
      REAL*8 BLWTO(20,3)                                        
      INTEGER*4 NK
      INTEGER*4 IMIN(50,3,100),IMAX(50,3,100)
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*4 NSIZE(3),NSIGMA,M,N                               
      INTEGER*4 QSMOD,BCOUNT,THALFB,RENEW
C----------------------------------------------------------------------
C     COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM6/PROPT,NK
      COMMON /COM7/SHVEL,TO
      COMMON /COM8/DMAX,DMIN,IMIN,IMAX
      COMMON /COM9/SLWTM,SLWTAM
      COMMON /COM10/BEDCON
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #        OPTBINP,OPTSINP,OPTERO
      COMMON /COM19/YMAX,DELY,VEL,EPSY,PSI,BETA
      COMMON /COM22/AWPM,D50M
      COMMON /COM23/BLWTAM,BLWTM,SUM1,SUM2
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
      COMMON /COM36/BLWTO,QSMOD,BCOUNT,THALFB,RENEW
C----------------------------------------------------------------------
C     INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 VOL,PROPVOL,TMP,MINISIF,CONV
      INTEGER*4 J,I,K,IS,IF,JN,NT
C----------------------------------------------------------------------
C     CALCULATE SEDIMENT INFLUX AT THE HEAD (BOUNDARY CONDITION)
C----------------------------------------------------------------------
      IF (NT .EQ. 0) THEN
         IF (OPTBINP .EQ. 1) THEN
            DO 2 J = 1, NSIGMA
               DO 1 I = 1, NSIZE(J)
                  BLWTM(I,J,1) = 0.0D0
                  BLWTAM(I,J,1) = 0.0D0
    1          CONTINUE
    2       CONTINUE
         ELSE IF (OPTBINP .EQ. 3) THEN
            DO 4 J = 1, NSIGMA
               READ (4, 100) (BLWTO(I,J), I = 1, NSIZE(J))
               DO 3 I = 1, NSIZE(J)
                  BLWTM(I,J,1) = BLWTO(I,J)
    3          CONTINUE
    4       CONTINUE
            BCOUNT = 0
         ENDIF
      ENDIF
      IF (QSMOD .EQ. 1) THEN
         IF (BCOUNT .GE. INT(RENEW*3600.0/DELT)) BCOUNT = 0
         DO 6 J = 1, NSIGMA
            DO 5 I = 1, NSIZE(J)
               CONV = LOG(0.5)
               BLWTM(I,J,1) = BLWTO(I,J)*EXP(CONV*(REAL(BCOUNT*DELT)/
     #            REAL(THALFB*3600.0/DELT)))
    5       CONTINUE
    6    CONTINUE
         BCOUNT = BCOUNT + 1
      ENDIF
C----------------------------------------------------------------------
C     START BEDLOAD CALC'S
C----------------------------------------------------------------------
      DO 13 JN = 2, M
         VOL = 0.0D0
C----------------------------------------------------------------------
C     VOL IS A PROPORTIONALITY FACTOR NECESSARY
C     TO CALCULATE THE PROPORTION OF VOLUME A CERTAIN AWPM REPRESENTS
C----------------------------------------------------------------------
         DO 8 J = 1, NSIGMA
            DO 7 I = 1, NSIZE(J)
               VOL = VOL + AWPM(I,J,JN)/SIGMA(J)
C----------------------------------------------------------------------
C     SET ALL VARIABLES TO ZERO
C----------------------------------------------------------------------
               BLWTM(I,J,JN) = 0.0D0
               BLWTAM(I,J,JN) = 0.0D0
               SLWTM(I,J,JN) = 0.0D0
               SLWTAM(I,J,JN) = 0.0D0
               BEDCON(I,J,JN) = 0.0D0
    7       CONTINUE
    8    CONTINUE
C----------------------------------------------------------------------
C     BEDLOAD TRANSPORT RATE CALCULATIONS
C----------------------------------------------------------------------
         DO 12 K = 1, NK
            DO 11 J = 1, NSIGMA
               IS = IMIN(K,J,JN)
               IF = IMAX(K,J,JN)
               IF (IF .LE. NSIZE(J)) THEN
                  IF (IS .LE. NSIZE(J)) THEN
C----------------------------------------------------------------------
C     IS=MINIMUM SIZE THAT WILL STAY IN THE MOVING BED LAYER
C     IF=MAXIMUM ENTRAINABLE SIZE
C----------------------------------------------------------------------
                     DO 9 I = IS, IF
                        IF (TO(K,JN) .GT. TC(I,J,JN)) THEN
                           PROPVOL = AWPM(I,J,JN)/(SIGMA(J)*VOL)
                           TMP = PROPT(K)*PROPVOL*10.0D0*(TO(K,JN)-TC(I,
     #                        J,JN))*(SHVEL(K,JN)-SHVELC(I,J,JN))
C----------------------------------------------------------------------
C     CONVERSION OF "WET" TRANSPORT RATE TO "DRY" TRANSPORT RATE
C     BLWTM(I,J,JN)= THE DRY WEIGHT TRANSPORTED THROUGH A METER
C            WIDTH OF STREAM IN ONE SECOND AS BEDLOAD TRANSPORT
C----------------------------------------------------------------------
                           TMP = TMP*SIGMA(J)/((SIGMA(J)-RHO)*G)
                           BLWTM(I,J,JN) = BLWTM(I,J,JN) + TMP
C----------------------------------------------------------------------
C     BEDLOAD WEIGHT TRANSPORT PER UNIT AREA CALCULATIONS
C----------------------------------------------------------------------
                           BLWTAM(I,J,JN) = BLWTAM(I,J,JN) + TMP/(A*(
     #                        SHVEL(K,JN)-SHVELC(I,J,JN)))
                        ENDIF
    9                CONTINUE
                  ENDIF
                  MINISIF = MIN(IS-1,IF)
                  DO 10 I = 1, INT(MINISIF)
                     IF (TO(K,JN) .GT. TC(I,J,JN)) THEN
                        PROPVOL = AWPM(I,J,JN)/(SIGMA(J)*VOL)
C----------------------------------------------------------------------
C     CALCULATE THE IMMERSED WT FLUX.
C     [TMP]=KG/S3=N/MS
C----------------------------------------------------------------------
                        TMP = PROPT(K)*PROPVOL*10.00D0*(TO(K,JN)-TC(I,J,
     #                     JN))*(SHVEL(K,JN)-SHVELC(I,J,JN))
C----------------------------------------------------------------------
C     CALCULATE THE MASS FLUX PER UNIT WIDTH
C     [TMP]=KG/MS
C     [SLWTM]=KG/MS
C----------------------------------------------------------------------
                        TMP = TMP*SIGMA(J)/((SIGMA(J)-RHO)*G)
                        SLWTM(I,J,JN) = SLWTM(I,J,JN) + TMP
C----------------------------------------------------------------------
C     CALCULATE THE DRY MASS OF GRAINS MOVING OVER UNIT BED AREA
C     [TMP]=KG/M2 (THE TMP HERE REFERS TO WHAT IS BELOW THIS COMMENT
C---------------------------------------------------------------------
                        TMP = TMP/(A*(SHVEL(K,JN)-SHVELC(I,J,JN)))
                        SLWTAM(I,J,JN) = SLWTAM(I,J,JN) + TMP
C----------------------------------------------------------------------
C     IF OPTSALT.EQ. 1 THEN USE BRIDGE FUNCTION TO CALCULATE
C     BEDCON=(KG/M3)=MASS PER UNIT VOLUME.
C     CALCULATE THE SUSPENDABLE AMOUNT PRESENT IN THE MOVING BED
C     LAYER:SLWTM [KG M-1 S-1] AND SLWTAM [KG M-2]
C     (SEE BRIDGE AND DOMINIC (1984) FOR YMAX EQUATION)
C----------------------------------------------------------------------
                        IF (OPTSALT .EQ. 1) THEN
                           YMAX(JN) = DIMID(I,J)*(2.53*DSQRT((TO(K,JN)-
     #                        TC(I,J,JN))/((SIGMA(J)-RHO)*G*DIMID(I,J)))
     #                        +0.5)
                        ELSE
C----------------------------------------------------------------------
C     ELSE OPTSALT.EQ.2 AND YMAX IS CALCULATED USING EINSTEIN(2*D50)
C----------------------------------------------------------------------
                           YMAX(JN) = D50M(JN)*2.0
                        ENDIF
                        BEDCON(I,J,JN) = BEDCON(I,J,JN) + TMP/YMAX(JN)
                     ENDIF
                     YMAX(1) = YMAX(2)
   10             CONTINUE
               ENDIF
   11       CONTINUE
   12    CONTINUE
   13 CONTINUE
      RETURN 
  100 FORMAT(9F8.0/9F8.0/2F8.0)
      END
