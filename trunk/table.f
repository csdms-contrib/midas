C******************************SUBROUTINE TABLE************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: PRINTS A TABLE OF SEDIMENT TRANSPORT VALUES
C----------------------------------------------------------------------
      SUBROUTINE TABLE(NT)
      IMPLICIT NONE
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100), 
     #AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 CONIM(20,3,100),CVELM(20,3,100),PHI                       
      REAL*8 DBDX(20,3,100),XFLUX(20,3,100),TXFLUX(100),TDBDX(100)    
      REAL*8 SDBDX(20,3,100),SDBDT(20,3,100),SXFLUX(20,3,100),       
     #STFLUX(20,3,100),STHICK(100)
      REAL*8 BLWTM(20,3,100),BLWTAM(20,3,100),SUM1(100),SUM2(100)   
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),  
     #SVELUM(100)
      REAL*8 SBLWTA(20,3,100),SBLWT(20,3,100),SCVELIA(20,3,100),  
     1SCVELI(20,3,100),SPRTHICK(20,3,100),TOTTHICK(20,3,100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                             
      INTEGER*2 OPTSUSP,OPTLNODE
      INTEGER*4      COUNTER
      INTEGER*4 NSIZE(3),NSIGMA,M,N                             
C----------------------------------------------------------------------
C            COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM15/CONIM,CVELM,PHI,OPTSUSP,OPTLNODE
      COMMON /COM16/DBDX,XFLUX,TDBDX,TXFLUX
      COMMON /COM17/SDBDX,SDBDT,SXFLUX,STFLUX,STHICK
      COMMON /COM23/BLWTAM,BLWTM,SUM1,SUM2
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM27/SBLWTA,SBLWT,SCVELIA,SCVELI,SPRTHICK,
     #TOTTHICK,COUNTER
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
C----------------------------------------------------------------------
C            DECLARATION STATEMENTS--EXTERNAL VARIABLES NO IN COMMON
C----------------------------------------------------------------------
      INTEGER*4 NT
      INTEGER*4 JN
C----------------------------------------------------------------------
C            DECLARATION STATEMENTS--      INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 TIME,ERROR,TIMEPRES,TIMEMIN1
      INTEGER*4 J,I
C----------------------------------------------------------------------
C          SDBDX=TIME AVERAGED MASS EXCHANGE BETWEEN BEDLOAD LAYER
C                  AND ACTIVE BED (KG/M2S)
C          SDBDT=INSTANTANEOUS MASS EXCHANGE RATE BETWEEN BEDLOAD LAYER
C                  AND ACTIVE BED (KG/M2S)
C          SXFLUX=TIME AVERAGED MASS EXCHANGE BETWEEN SUSPENDED LOAD
C                  AND BED  (KG/M2S)
C          STFLUX=INSTANTANEOUS MASS EXCHANGE RATE BETWEEN BEDLOAD
C                 AND ACTIVE BED
C          SBLWT=INSTANTANEOUS BEDLOAD FLUX (KG/MS)
C          SBLWTA=AVERAGE BEDLOAD FLUX
C          SCVEL=INSTANTANEOUS SUSPENDED LOAD FLUX
C          SCVELIA=AVERAGE SUSPENDED LOAD FLUX
C          SPRTHICK=THICKNESS ERODED OR DEPOSITED OVER THE TIMESTEP
C          TOTTHICK=TOTAL THICKNESSES ERODED OR DEPOSITED SINCE THE
C                START OF SIMULATION.
C----------------------------------------------------------------------
      DO 7 JN = 2, M
         IF (NT .EQ. 0) THEN
            COUNTER = 1
            DO 2 J = 1, NSIGMA
               DO 1 I = 1, NSIZE(J)
                  SDBDX(I,J,JN) = 0.0D0
                  SDBDT(I,J,JN) = 0.0D0
                  SXFLUX(I,J,JN) = 0.0D0
                  STFLUX(I,J,JN) = 0.0D0
                  SBLWT(I,J,JN) = 0.0D0
                  SCVELI(I,J,JN) = 0.0D0
                  SBLWTA(I,J,JN) = 0.0D0
                  SCVELIA(I,J,JN) = 0.0D0
                  SPRTHICK(I,J,JN) = 0.0D0
                  TOTTHICK(I,J,JN) = 0.0D0
    1          CONTINUE
    2       CONTINUE
         ENDIF
C----------------------------------------------------------------------
C                  IF PREVIOUS TIMESTEP IS A TPRINT THEN CLEAR THE
C                  AVERAGING ARRAYS SO THAT A NEW TIME AVERAGE CAN
C                  BE DEVELOPED
C----------------------------------------------------------------------
         TIME = FLOAT(NT-1)*DELT/TPRINT + 0.000001
         ERROR = TIME - FLOAT(INT(TIME))
         IF (ERROR .LT. .001) THEN
            COUNTER = 1
            DO 4 J = 1, NSIGMA
               DO 3 I = 1, NSIZE(J)
                  SDBDX(I,J,JN) = 0.0D0
                  SXFLUX(I,J,JN) = 0.0D0
                  SBLWTA(I,J,JN) = 0.0D0
                  SCVELIA(I,J,JN) = 0.0D0
                  SPRTHICK(I,J,JN) = 0.0D0
    3          CONTINUE
    4       CONTINUE
         ENDIF
         TIMEPRES = FLOAT(COUNTER)
         TIMEMIN1 = FLOAT(COUNTER-1)
         DO 6 J = 1, NSIGMA
            DO 5 I = 1, NSIZE(J)
               SDBDX(I,J,JN) = ((SDBDX(I,J,JN)*TIMEMIN1*WIDM(JN)*DELT+
     #            DBDX(I,J,JN))/(TIMEPRES*DELT))/WIDM(JN)
               SDBDT(I,J,JN) = (DBDX(I,J,JN)/DELT)/WIDM(JN)
               SXFLUX(I,J,JN) = ((SXFLUX(I,J,JN)*TIMEMIN1*WIDM(JN)*DELT+
     #            XFLUX(I,J,JN))/(TIMEPRES*DELT))/WIDM(JN)
               STFLUX(I,J,JN) = (XFLUX(I,J,JN)/DELT)/WIDM(JN)
               SBLWT(I,J,JN) = BLWTM(I,J,JN)
               SBLWTA(I,J,JN) = (TIMEMIN1*SBLWTA(I,J,JN)+BLWTM(I,J,JN))/
     #            TIMEPRES
               SCVELI(I,J,JN) = CVELM(I,J,JN)
               SCVELIA(I,J,JN) = (TIMEMIN1*SCVELIA(I,J,JN)+CVELM(I,J,JN)
     #            )/TIMEPRES
               SPRTHICK(I,J,JN) = SPRTHICK(I,J,JN) + PRTHICK(I,J,JN)
               TOTTHICK(I,J,JN) = TOTTHICK(I,J,JN) + PRTHICK(I,J,JN)
    5       CONTINUE
    6    CONTINUE
         COUNTER = COUNTER + 1
    7 CONTINUE
      RETURN 
      END
