C**************************SUBROUTINE WRFLDTA**************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: WRITES INFORMATION ABOUT THE CHANNEL FORM AND
C                   PROFILE TO OUTPUT
C----------------------------------------------------------------------
      SUBROUTINE WRFLDTA(NT)
      IMPLICIT NONE
      REAL*8 Q,MANN(100),DINIT,FLOWEL
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100), 
     #AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 AWPM(20,3,100),D50M(100)                                  
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),     
     #SVELUM(100)
      REAL*8   STRDATA(10,100)                                       
      REAL*8 DELX,DELT,XPRINT,TPRINT                                
      INTEGER*4 NEXT,ICONT
      INTEGER*2 PRINT(7)
      INTEGER*4 NSIZE(3),NSIGMA,M,N                                
      CHARACTER*72 RUN                                            
C----------------------------------------------------------------------
C            COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM11/Q,MANN,DINIT,ICONT,NEXT,FLOWEL
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM21/RUN,PRINT
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM26/STRDATA
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--EXTERNAL VARIABLES--NOT IN COMMON BLOCKS
C----------------------------------------------------------------------
      INTEGER*4 NT
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 POS
      REAL*8 VV
      REAL*8 SF
      REAL*8      RCH
      REAL*8 FR
      REAL*8      POSSTEP
      INTEGER*4 JN
      INTEGER*4 COUNTER
C----------------------------------------------------------------------
C      DEFINE PARAMETERS
C----------------------------------------------------------------------
      IF (NT .EQ. 0) THEN
         WRITE (8, 100) RUN, Q, ICONT, DELX, DELT, XPRINT, TPRINT
  100    FORMAT (1X, 'SEDIMENT TRANSPORT MODEL FOR HETEROGENEOUS SIZE',
     #                  'DENSITY MIXTURES',/
     #                  '                          VERSION 4.0',//
     #                  ' TITLE : ',A72,//
     #                  '----------------------------------------',/
     #                  ' THE DISCHARGE IS             ',F8.3,' M3/S',/
     #                  ' THE NUMBER OF NODES (ICONT)  ',I4,/
     #                  ' THE DIST BETW. NODES(DELX)   ',F8.3, ' M',/
     #                  ' THE SIZE OF TIMESTEP (DELT)  ',E10.4, ' S',/
     #                  ' PRINT EVERY (XPRINT)         ',E10.4, ' M',/
     #                  ' PRINT EVERY (TPRINT)         ',E10.4, ' S',/
     #                  '----------------------------------------',//)
      ENDIF
      WRITE (8, 101) NT
  101          FORMAT (1X,' AT TIMESTEP',I9//
     #  ' DISTANCE   BED ELEVATION  FLOW DEPTH ',
     #  ' CHANNEL WIDTH  HYD.RAD   MEAN VEL.    FROUDE#  ',
     #  '  FRICTION SLOPE   ',
     #  ' T eff.        D50',/
     #  '      (M)             (M)         (M)       ',
     #  '     (M)      (M)       (M/S)                  ',
     #  '                 (Pa)         (M)')
      COUNTER = 0
      DO 1 JN = 2, M
C---------------------------------------------------------------------
C            CHECK TO SEE IF ITS IN POSITION TO PRINT
C---------------------------------------------------------------------
         POSSTEP = (JN-2)*DELX/XPRINT + 0.00001
         IF (POSSTEP - FLOAT(INT(POSSTEP)) .LE. 0.0001) THEN
            POS = STRDATA(4,JN)
            RCH = STRDATA(1,JN)
            SF = STRDATA(3,JN)
            VV = STRDATA(5,JN)
            FR = STRDATA(2,JN)
            WRITE (8, 104) POS, ELEVM(JN), DCHM(JN), WIDM(JN), RCH, VV, 
     #         FR, SF, TOXM(JN), D50M(JN)
C----------------------------------------------------------------------
C               WRITE TO PLOT FILE
C----------------------------------------------------------------------
            COUNTER = COUNTER + 1
         ENDIF
    1 CONTINUE
C 102           FORMAT(2E15.8)
      WRITE (8, 103) 
  103                 FORMAT(' ','----------------------------------',
     #                        '-------------------------------------',
     #                        '------------------------------',
     #                          '------------------------------')
  104     FORMAT(E10.4,7X,F8.3,4X,F8.3,7X,F8.3,3X,F6.3,2X,E10.3,4X,
     #           F6.3,8X,E10.3,3X,E10.3,2X,E10.3)
      RETURN 
      END
