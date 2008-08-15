C*********************** SUBROUTINE RESULT.F **************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: WRITES TO FILE RESULTA THE BEDLOAD TRANSPORT
C                   INFORMATION AFTER THE ACTIVE BED LAYER IS TAKEN
C                   INTO ACCONT
C----------------------------------------------------------------------
      SUBROUTINE RESULT(NT)
      IMPLICIT NONE
      REAL*8 DIMID(20,3),SIGMA(3)                                   
      REAL*8 Q,MANN(100),DINIT,FLOWEL                              
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100),  
     #AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 SDBDX(20,3,100),SDBDT(20,3,100),SXFLUX(20,3,100),          
     #STFLUX(20,3,100),STHICK(100)
      CHARACTER*72 RUN                                                 
      REAL*8 AWPM(20,3,100),D50M(100)                                 
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),    
     #SVELUM(100)
      REAL*8 SBLWTA(20,3,100),SBLWT(20,3,100),SCVELIA(20,3,100),    
     #SCVELI(20,3,100),SPRTHICK(20,3,100),TOTTHICK(20,3,100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                               
      REAL*8  PLOTPRINT                                           
      INTEGER*4 NEXT,ICONT
      INTEGER*2 PRINT(7)
      INTEGER*4      COUNTER
      INTEGER*4 NSIZE(3),NSIGMA,M,N                              
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM11/Q,MANN,DINIT,ICONT,NEXT,FLOWEL
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM17/SDBDX,SDBDT,SXFLUX,STFLUX,STHICK
      COMMON /COM21/RUN,PRINT
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM27/SBLWTA,SBLWT,SCVELIA,SCVELI,COUNTER,SPRTHICK,
     #TOTTHICK
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
      COMMON /COM30/PLOTPRINT
C----------------------------------------------------------------------
C      DECLARATION STATEMENTS--EXTERNAL VARIABLES--NOT IN COMMON BLOCKS
C----------------------------------------------------------------------
      INTEGER*4 NT
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 TOTSDBDX,TOTSBDT,TOTSXFLUX,TOTSTFLUX,TOTSBLWTA
      REAL*8 TOTBLWT,TOTSCVELIA,TOTSCVELI,TOTSPRTHICK
      REAL*8 TOTALTHICK,INSTTHIC,STEP,TIME
      INTEGER*4 I,J,JN
C*******************************************************************
C     IF AT FIRST TIMESTEP THEN OPEN THE FILES.
C---------------------------------------------------------------------
      IF (NT .EQ. 0) THEN
         IF (PRINT(3) .EQ. 1) THEN
C----------------------------------------------------------------------
C     WRITE HEADER TO SUMMARY DATA FILE (UNIT=12)
C----------------------------------------------------------------------
            WRITE (12, 100) RUN, ICONT, DELX, DELT, XPRINT, TPRINT
  100    FORMAT (1X, 'SEDIMENT TRANSPORT MODEL FOR HETEROGENEOUS SIZE',
     #                 'DENSITY MIXTURES',/
     #                 '                          VERSION 3.0',//
     #                 ' TITLE : ',A72,//
     #                 '------------------------------------------',/
     #                 ' THE NUMBER OF NODES (ICONT)  ',I4,/
     #                 ' THE DIST BETW. NODES(DELX)   ',F8.3, ' M',/
     #                 ' THE SIZE OF TIMESTEP (DELT)  ',F8.3, ' S',/
     #                 ' PRINT EVERY (XPRINT)         ',F8.3, ' M',/
     #                 ' PRINT EVERY (TPRINT)         ',F10.0, ' S',/
     #                 '------------------------------------------',//)
         ENDIF
         IF (PRINT(4) .EQ. 1) THEN
C----------------------------------------------------------------------
C                  WRITE HEADER TO DETAILED DATA FILE (UNIT=13)
C----------------------------------------------------------------------
            WRITE (13, 101) RUN, Q, ICONT, DELX, DELT, XPRINT, TPRINT
  101    FORMAT (1X, 'SEDIMENT TRANSPORT MODEL FOR HETEROGENEOUS SIZE',
     #                 'DENSITY MIXTURES',/
     #                 '                          VERSION 3.0',//
     #                 ' TITLE : ',A72,//
     #                 '------------------------------------------',/
     #                 ' THE DISCHARGE IS             ',F8.3,' M3/S',/
     #                 ' THE NUMBER OF NODES (ICONT)  ',I4,/
     #                 ' THE DIST BETW. NODES(DELX)   ',F8.3, ' M',/
     #                 ' THE SIZE OF TIMESTEP (DELT)  ',F8.3, ' S',/
     #                 ' PRINT EVERY (XPRINT)         ',F8.3, ' M',/
     #                 ' PRINT EVERY (TPRINT)         ',F10.0, ' S',/
     #                 '------------------------------------------',//)
         ENDIF
      ENDIF
C----------------------------------------------------------------------
C            YOU HAVE NOW FINISHED WRITING THE HEADERS TO FILE 13 AND
C            FILE 12 IF THEY WERE TO BE WRITTEN TO.
C----------------------------------------------------------------------
      IF (PRINT(3) .EQ. 1) THEN
C----------------------------------------------------------------------
C                  WRITE MORE HEADER TO SUMMARY FILE
C----------------------------------------------------------------------
         WRITE (12, 102) 
  102              FORMAT(/,'  TIMESTEP  POSITION     AVERAGE    ',
     #       ' AVERAGE     ',
     #       ' AVERAGE       AVERAGE  ',
     #       '    CUMULATIVE    CUMULATIVE    CUMULATIVE   DFIFTY  ',)
         WRITE (12, 103) 
  103           FORMAT('                         BEDLOAD',
     #            '     SUSP. LOAD ',
     #            '    BEDLOAD       SUSP. LOAD ',
     #            ' AMT.ERODED(-) AMT.ERODED(-) AMT.ERODED(-)       ',)
         WRITE (12, 104) 
  104           FORMAT('                         SPATIAL',
     #            '       SPATIAL   ',
     #            '    MASS TRANS   MASS TRANS',
     #            '  OVER DELT     OVER TPRINT     OVER N*DELT       ',)
         WRITE (12, 105) 
  105           FORMAT('                         GRADIENT',
     #            '      GRADIENT  ',
     #             '     RATE PER     RATE PER',
     #            '                                              ',)
         WRITE (12, 106) 
  106           FORMAT('                        PER UNIT A',
     #            '     PER UNIT A  ',
     #            'UNIT WIDTH     UNIT WIDTH',)
         WRITE (12, 107) 
  107            FORMAT('               (M)       (KG/M2S) ',
     #       '     (KG/M2S) ',
     #       '     (KG/MS)       (KG/MS)',
     #       '            (M)              (M)            (M)'
     #       '           (M)')
         WRITE (12, 117) 
      ENDIF
      IF (PRINT(4) .EQ. 1) THEN
         WRITE (13, 118) 
         WRITE (13, 108) NT
  108           FORMAT(1X,'AT TIMESTEP',I5,/)
      ENDIF
      DO 5 JN = 2, M
         STEP = (JN-2)*DELX/XPRINT + 0.00001
         IF (FLOAT(STEP) - FLOAT(INT(STEP)) .LT. 0.0001) THEN
            IF (PRINT(4) .EQ. 1) THEN
               TIME = DELT*NT/(PRINT(6)*TPRINT) + 0.00001
               IF (TIME - FLOAT(INT(TIME)) .LE. 0.0001) THEN
                  WRITE (13, 109) FLOAT(JN-2)*DELX
  109            FORMAT(//'RESULTS AT ',E10.4,'METERS DOWN REACH')
                  WRITE (13, 110) 
  110              FORMAT(/,' GRAIN SIZE AVERAGE    INSTANT',
     #  '   AVERAGE    INSTANT',
     #  '    AVERAGE    INSTANT    AVERAGE    INSTANT',
     #  '   CUMULATIVE  CUMULATIVE  CUMULATIVE',)
C-----------------------------------------------------------------
                  WRITE (13, 111) 
  111           FORMAT(' MIDPOINT   BEDLOAD    BEDLOAD',
     #            '   SUSP. LOAD SUSP. LOAD',
     #            ' BEDLOAD    BEDLOAD SUSP. LOAD SUSP. LOAD',
     #            ' AMT.ERODED(-) AMT.ERODED(-) AMT.ERODED(-)',)
                  WRITE (13, 112) 
  112           FORMAT('            SPATIAL    SPATIAL',
     #            '   SPATIAL    SPATIAL',
     #            '   MASS TRANS MASS TRANS  MASS TRANS MASS TRANS',
     #            ' OVER DELT   OVER TPRINT OVER N*DELT',)
                  WRITE (13, 113) 
  113           FORMAT('            GRADIENT   GRADIENT',
     #            '  GRADIENT   GRADIENT  ',
     #             ' RATE PER   RATE PER    RATE PER   RATE PER',
     #       '    ',)
                  WRITE (13, 114) 
  114           FORMAT('            PER UNIT A PER UNIT A',
     #            'PER UNIT A PER UNIT A  ',
     #            'UNIT WIDTH UNIT WIDTH  UNIT WIDTH UNIT WIDTH',)
                  WRITE (13, 115) 
  115            FORMAT('   (M)       (KG/M2S) (KG/M2S) ',
     #            '   (KG/M2S)    (KG/M2S) ',
     #            ' (KG/MS)     (KG/MS)   (KG/MS)     (KG/MS)',
     #            '      (M)        (M)        (M)    ')
                  WRITE (13, 117) 
C---------------------------------------------------------------------
C                  WRITE THE FLUX VALUES TO THE DETAILED FILE
C---------------------------------------------------------------------
                  DO 2 J = 1, NSIGMA
                     DO 1 I = 1, NSIZE(J)
                        WRITE (13, 116) DIMID(I,J), SDBDX(I,J,JN), SDBDT
     #                     (I,J,JN), SXFLUX(I,J,JN), STFLUX(I,J,JN), 
     #                     SBLWTA(I,J,JN), SBLWT(I,J,JN), SCVELIA(I,J,JN
     #                     ), SCVELI(I,J,JN),  - PRTHICK(I,J,JN),  - 
     #                     SPRTHICK(I,J,JN),  - TOTTHICK(I,J,JN)
    1                CONTINUE
    2             CONTINUE
               ENDIF
            ENDIF
  116 FORMAT(12(1X,E10.4))
C----------------------------------------------------------------------
C            PREPARE TO SUMMARIZE THE FLUX VALUES FROM EACH SIZE
C            FRACTION. FIRST SET VALUES TO ZERO
C----------------------------------------------------------------------
            TOTSDBDX = 0.0D0
            TOTSBDT = 0.0D0
            TOTSXFLUX = 0.0D0
            TOTSTFLUX = 0.0D0
            TOTSBLWTA = 0.0D0
            TOTBLWT = 0.0D0
            TOTSCVELI = 0.0D0
            TOTSCVELIA = 0.0D0
            TOTSPRTHICK = 0.0D0
            TOTALTHICK = 0.0D0
            INSTTHIC = 0.0D0
C----------------------------------------------------------------------
C            SUM THE FLUX VALUES OVER EACH SIZE FRACTION
C----------------------------------------------------------------------
            DO 4 J = 1, NSIGMA
               DO 3 I = 1, NSIZE(J)
                  TOTSDBDX = TOTSDBDX + SDBDX(I,J,JN)
                  TOTSBDT = TOTSBDT + SDBDT(I,J,JN)
                  TOTSXFLUX = TOTSXFLUX + SXFLUX(I,J,JN)
                  TOTSTFLUX = TOTSTFLUX + STFLUX(I,J,JN)
                  TOTSBLWTA = TOTSBLWTA + SBLWTA(I,J,JN)
                  TOTBLWT = TOTBLWT + SBLWT(I,J,JN)
                  TOTSCVELIA = TOTSCVELIA + SCVELIA(I,J,JN)
                  TOTSCVELI = TOTSCVELI + SCVELI(I,J,JN)
                  TOTSPRTHICK = TOTSPRTHICK + SPRTHICK(I,J,JN)
                  TOTALTHICK = TOTALTHICK + TOTTHICK(I,J,JN)
                  INSTTHIC = INSTTHIC + PRTHICK(I,J,JN)
    3          CONTINUE
    4       CONTINUE
  117 FORMAT('','-------------------------------------------------------
     #------------------------------------------------------------------
     #-----------')
  118           FORMAT(/,'************************************',
     #                  '*************************************',
     #                  '************************************',
     #                  '********************')
            IF (PRINT(4) .EQ. 1) THEN
               TIME = DELT*NT/(PRINT(6)*TPRINT) + 0.00001
               IF (TIME - FLOAT(INT(TIME)) .LE. 0.0001) THEN
                  WRITE (13, 117) 
                  WRITE (13, 120) TOTSDBDX, TOTSBDT, TOTSXFLUX, 
     #               TOTSTFLUX, TOTSBLWTA, TOTBLWT, TOTSCVELIA, 
     #               TOTSCVELI,  - INSTTHIC,  - TOTSPRTHICK,  - 
     #               TOTALTHICK
               ENDIF
            ENDIF
            IF (PRINT(3) .EQ. 1) WRITE (12, 119) NT, FLOAT(JN-2)*DELX, 
     #         TOTSDBDX, TOTSXFLUX, TOTSBLWTA, TOTSCVELIA,  - INSTTHIC, 
     #          - TOTSPRTHICK,  - TOTALTHICK, D50M(JN)
         ENDIF
    5 CONTINUE
  119     FORMAT(2X,I3,9(2X,E10.4))
  120     FORMAT(6X,'TOTALS:',11(1X,E10.4))
      RETURN 
      END
