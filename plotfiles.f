C*********************SUBROUTINE PLOTFILES*****************************
C----------------------------------------------------------------------
C THIS SUBROUTINE: WRITES PLOT DATA TO FILES
C----------------------------------------------------------------------
      SUBROUTINE PLOTFILES
      IMPLICIT NONE
      REAL*8 DIFEND(20,3),DIMIN(3)                              
      REAL*8 AWPM(20,3,100),D50M(100)                          
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),   
     #SVELUM(100)
      REAL*8 SBLWTA(20,3,100),SBLWT(20,3,100),SCVELIA(20,3,100),   
     1SCVELI(20,3,100),SPRTHICK(20,3,100),TOTTHICK(20,3,100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                             
      REAL*8  PLOTPRINT                                         
      REAL*8 RLONG,TFILE                                       
      INTEGER*4      COUNTER
      INTEGER*4 NSIZE(3),NSIGMA,M,N                               
      INTEGER*4 NT,OPTDIST(3),IFILE,OPTCONT,NTIM
      CHARACTER*24  FOUT
      CHARACTER*4 RUNNAME
C---------------------------------------------------------------------
C     COMMON DECLARATION
C---------------------------------------------------------------------
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM27/SBLWTA,SBLWT,SCVELIA,SCVELI,COUNTER,SPRTHICK,
     #              TOTTHICK
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
      COMMON /COM30/PLOTPRINT
      COMMON /COM33/RLONG,TFILE,NTIM,NT,OPTDIST,
     #              IFILE,OPTCONT,FOUT,RUNNAME
C---------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C---------------------------------------------------------------------
      INTEGER*4 JN,I,J,NUMNODES
      REAL*8 POS,CUMWT,TEMP
      REAL*8 TOTSBLWTA,TOTSCVELIA
C---------------------------------------------------------------------
      IF (NT .EQ. 0) THEN
         NUMNODES = INT((M-2)*DELX/PLOTPRINT) + 1
         WRITE (18, *) NUMNODES
         WRITE (20, *) NUMNODES
         DO 2 J = 1, NSIGMA
            WRITE (12 + J, *) NUMNODES, NSIZE(J)
            DO 1 I = 1, NSIZE(J)
               WRITE (12 + J, *) DIFEND(I,J)*1000
    1       CONTINUE
    2    CONTINUE
      ENDIF
  100  FORMAT(1X,I6)
  101  FORMAT(1X,I6,2X,I3)
      WRITE (20, *) CHAR(39), NT, CHAR(39)
      WRITE (22, *) CHAR(39), NT, CHAR(39)
      DO 3 J = 1, NSIGMA
         WRITE (12 + J, *) CHAR(39), NT, CHAR(39)
    3 CONTINUE
  102 FORMAT(1X,A1,I6,A1)
      TOTSBLWTA = 0.0D0
      TOTSCVELIA = 0.0D0
      DO 5 J = 1, NSIGMA
         DO 4 I = 1, NSIZE(J)
            TOTSBLWTA = TOTSBLWTA + SBLWTA(I,J,M)
            TOTSCVELIA = TOTSCVELIA + SCVELIA(I,J,M)
            WRITE (22, *) SBLWTA(I,J,M)
    4    CONTINUE
    5 CONTINUE
      WRITE (21, *) NT, TOTSBLWTA
      DO 8 JN = 2, M
         POS = DELX*(JN-2)
         TEMP = POS/PLOTPRINT - INT(POS/PLOTPRINT+0.0001)
         IF (TEMP .LE. 0.01) THEN
            WRITE (18, *) POS, ELEVM(JN), D50M(JN)*1000
            WRITE (20, *) POS, DCHM(JN)
  103     FORMAT(1X,F10.3,2X,F8.3)
            CUMWT = 0.0
            DO 7 J = 1, NSIGMA
               WRITE (12 + J, *) POS
               DO 6 I = NSIZE(J), 1, -1
                  WRITE (12 + J, 105) AWPM(I,J,JN)*100
    6          CONTINUE
    7       CONTINUE
  104     FORMAT(/)
         ENDIF
  105  FORMAT(1X,F8.3,2X,F8.2)
  106  FORMAT(1X,F8.0)
  107  FORMAT(1X,3F10.3)
    8 CONTINUE
      RETURN 
      END
