C***************************SUBROUTINE FLDTA***************************
C  THIS SUBROUTINE: 1) READS INFORMATION ABOUT AN OPEN CHANNEL AND
C                   2) CALCULATES THE FLOW CHARACTERISTICS BASED ON
C                      THE GRADUALLY VARIED FLOW EQUATION. A COMPLETE
C                      DERIVATION OF THE METHOD CAN BE FOUND IN
C                      HENDERSON (1966)
C
C  CREATORS: P. B. FLEMINGS & K. VOGEL
C----------------------------------------------------------------------
      SUBROUTINE FLDTA
      IMPLICIT NONE
      REAL*8 RHO,G
      REAL*8 QO,MANN(100),DINIT,FLOWEL
      INTEGER*4 NEXT,ICONT
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),
     #SVELUM(100)
      REAL*8   STRDATA(10,100)                                  
      INTEGER*4 NSIZE(3),NSIGMA,M,N                            
      REAL*8 DELX,DELT,XPRINT,TPRINT                          
      REAL*8 WIDTOP,DELWID,ELEVTOP,SLOPE                     
      INTEGER*4 OPTEL,OPTWID
      REAL*8 ELEVOM(100)   
      REAL*8 RLONG,TFILE 
      INTEGER*4 NT,OPTDIST,IFILE,OPTCONT,NTIM
      CHARACTER*24  FOUT
      CHARACTER*4 RUNNAME
      REAL*8 QFAC       
      INTEGER*4 QTMOD,THALFQ,QCOUNT,QSPMOD
C----------------------------------------------------------------------
C        COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM11/QO,MANN,DINIT,ICONT,NEXT,FLOWEL
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM26/STRDATA
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
      COMMON /COM31/OPTEL,OPTWID,WIDTOP,DELWID,ELEVTOP,SLOPE
      COMMON /COM32/ELEVOM
      COMMON /COM33/ RLONG,TFILE,NTIM,NT,OPTDIST,
     #IFILE,OPTCONT,FOUT,RUNNAME
      COMMON /COM35/ QTMOD,THALFQ,QCOUNT,QSPMOD,QFAC
C----------------------------------------------------------------------
C            DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 DCONT,A,P,V,SF1,H1,H2,H21,RCH,D2,A2,P2,V2,CONST,DELTAD
      REAL*8 RCH2,SF2,POS,FR,SLOP,YCRIT,DUR,QN(100),CONV,CONST2
      REAL*8 DMIN1
      INTEGER*4 JN,K,END,STEP,NODE,SIGN,COUNT
      CHARACTER*1 DUMMY
C----------------------------------------------------------------------
C  DEFINE PARAMETERS
C----------------------------------------------------------------------
      IF (NT .EQ. 0) THEN
         READ (9, 102) DCONT
         ICONT = INT(DCONT/DELX+2.1)
         READ (9, 102) ELEVM(ICONT)
         READ (9, 101) DUMMY
         READ (9, 100) DINIT, QO, DUR, OPTCONT, QTMOD, THALFQ
         QCOUNT = 0
         DUR = DUR*3600.
         IF (OPTCONT .EQ. 2) FLOWEL = DINIT + ELEVM(ICONT)
         NEXT = INT(DUR/DELT+1.5)
      ENDIF
      IF (NT .GE. NEXT) THEN
         READ (9, 100) DINIT, QO, DUR, OPTCONT, QTMOD, THALFQ
         QCOUNT = 0
         DUR = DUR*3600.
         NEXT = NEXT + INT(DUR/DELT+0.5)
         IF (OPTCONT .EQ. 2) FLOWEL = DINIT + ELEVM(ICONT)
      ENDIF
      ELEVM(1) = 2.*ELEVM(2) - ELEVM(3)
C----------------------------------------------------------------------
C      MODIFY QO IF IT IS A FUNCTION OF TIME ( QTMOD = 1 )
C----------------------------------------------------------------------
      QN(2) = QO
      IF (QTMOD .EQ. 1) THEN
         CONV = LOG(0.5)
         QN(2) = QO*EXP(CONV*FLOAT(QCOUNT*DELT)/FLOAT(THALFQ))
         QCOUNT = QCOUNT + 1
      ENDIF
      IF (QSPMOD .EQ. 1) THEN
         QN(ICONT) = QN(2)*(ABS(ICONT+1)*DELX)**QFAC
      ELSE
         QN(ICONT) = QN(2)
      ENDIF
C----------------------------------------------------------------------
C    DO UNIFORM FLOW CALCULATIONS IF OPTCONT=1
C----------------------------------------------------------------------
      IF (OPTCONT .EQ. 1) THEN
         YCRIT = ((QN(ICONT)/WIDM(ICONT))**2./9.81)**(1./3.)
C----------------------------------------------------------------------
C CALCULATE DINIT USING A NEWTON-RAPHSON SOLUTION FOR MANNING'S EQ.
C----------------------------------------------------------------------
         SLOP = (ELEVM(ICONT-1)-ELEVM(ICONT))/DELX
         DINIT = ((QN(ICONT)*MANN(ICONT))/(WIDM(ICONT)*SQRT(SLOP)))**(3.
     #      /5.)
         IF (DINIT .LT. YCRIT) DINIT = YCRIT
         IF (2.0*DINIT/WIDM(ICONT) .GT. 0.00001) THEN
            K = 0
            CONST = QN(ICONT)*MANN(ICONT)/(WIDM(ICONT)**(5./3.)*SQRT(
     #         SLOP))
    1       CONTINUE
            CONST2 = DINIT**(5./3.)/(WIDM(ICONT)+2.*DINIT)**(2./3.)
            DELTAD = ABS(CONST-CONST2)
            DELTAD = (DELTAD*WIDM(ICONT)**(2./3.))**(3./5.)
            IF (DELTAD.LT.1E-9 .OR. K.EQ.20) GO TO 2
            SIGN = (CONST-CONST2)/ABS(CONST-CONST2)
            DINIT = (DINIT**(5./3.)+SIGN*DELTAD**(5./3.))**(3./5.)
            K = K + 1
            GO TO 1
         ENDIF
      ELSE IF (OPTCONT .EQ. 2) THEN
         DINIT = FLOWEL - ELEVM(ICONT)
      ENDIF
C----------------------------------------------------------------------
C       CALCULATE INITIAL CROSS-SECT AREA, THE WETTED
C       PERIMETER(P), THE MEAN VEL. (V), THE HYD RAD (RCH)
C       AND THE DEPTH OF THE WATER (DCHM).
C----------------------------------------------------------------------
    2 CONTINUE
      A = DINIT*WIDM(ICONT)
      P = WIDM(ICONT) + 2.*DINIT
      V = QN(ICONT)/A
      RCH = A/P
      DCHM(ICONT) = DINIT
C----------------------------------------------------------------------
C            CALCULATE THE ENERGY SLOPE (SFL) AND THE HYDRAULIC HEAD
C            (H1), OR IF NOT AT CONTROL POINT TAKE THE PREVIOUS NODES
C            VALUES FOR THESE VALUES.
C----------------------------------------------------------------------
      DO 5 K = 1, 2
         IF (K .EQ. 1) THEN
            END = M - 1
            STEP = 1
         ELSE
            END = 3
            STEP = -1
         ENDIF
         DO 4 JN = ICONT, END, STEP
            COUNT = 0
            IF (JN .EQ. ICONT) THEN
               SF1 = ((MANN(JN)*V)/RCH**0.666667)**2.
               H1 = ELEVM(JN) + DINIT + V**2./(2.*9.81)
               DMIN1 = DINIT
            ELSE
               SF1 = SF2
               H1 = H2
            ENDIF
C----------------------------------------------------------------------
C              GUESS A NEW STREAM DEPTH AT THE NEXT NODE STARTING
C            WITH THE DEPTH AT THE PRESENT NODE. THE ONLY NEW DATA
C            USED IS THE STREAM WIDTH AT THE NEW NODE.
C----------------------------------------------------------------------
            D2 = DMIN1
            NODE = JN + STEP
            IF (QSPMOD .EQ. 1) THEN
               QN(NODE) = QN(2)*EXP(QFAC*(ABS(JN-2)*DELX))
            ELSE
               QN(NODE) = QN(2)
            ENDIF
    3       CONTINUE
            IF (D2 .LE. 1.0D-9) D2 = 1.0D-3
            A2 = D2*WIDM(NODE)
            P2 = WIDM(NODE) + 2.*D2
            V2 = QN(NODE)/A2
            RCH2 = A2/P2
C----------------------------------------------------------------------
C            CALCULATE ENERGY BY ENERGY EQUATION AND BY ENERGY SLOPE
C            THEN RECALCULATE UNTIL THE DIFFERENCE IS SMALL
C----------------------------------------------------------------------
            SF2 = ((MANN(NODE)*V2)/RCH2**0.666667)**2.
            H2 = ELEVM(NODE) + D2 + V2**2./(2.*9.81)
            H21 = H1 - STEP*0.5*DELX*(SF1+SF2)
            IF (ABS(H2-H21).GE.1E-6 .AND. COUNT.LE.50) THEN
               COUNT = COUNT + 1
               D2 = D2 + (H21-H2)/(1-V2**2./(9.81*D2)+(3.*SF2*DELX)/(2.*
     #            RCH2))
               GO TO 3
            ELSE
               DCHM(NODE) = D2
            ENDIF
C         IF((STEP.EQ.-1).AND.((V2**2.).GT.(9.81*D2)))THEN
C            DCHM(NODE)=((QN(NODE)/WIDM(NODE))**2./9.81)**(1./3.)
C            A2=DCHM(NODE)*WIDM(NODE)
C            P2=WIDM(NODE)+2.*DCHM(NODE)
C            V2=QN(NODE)/A2
C            RCH2=A2/P2
C            SF2=((MANN(NODE)*V2)/(RCH2**0.666667))**2.
C            H2=ELEVM(NODE)+D2+V2**2./(2.*9.81)
C         ENDIF
    4    CONTINUE
    5 CONTINUE
C----------------------------------------------------------------------
C      WRITE TO THE ARRAY STRDATA THE VALUES OF THE ENERGY SLOPE(SF)
C            THE HYDRAULIC RADIUS (RCH), THE POSITION(POS).
C----------------------------------------------------------------------
      DCHM(1) = DCHM(2)
      DCHM(M+1) = DCHM(M)
      QN(1) = QN(2)
      DO 6 JN = 1, M
         DCHM(JN) = INT((DCHM(JN)+0.0005)*1000.)/1000.
         A = DCHM(JN)*WIDM(JN)
         P = WIDM(JN) + 2*DCHM(JN)
         V = QN(JN)/A
         RCH = A/P
         SF1 = ((MANN(JN)*V)/RCH**0.666667)**2.
         POS = FLOAT(JN-2)*DELX
         FR = V/(G*DCHM(JN))**.5
         STRDATA(1,JN) = RCH
         STRDATA(2,JN) = FR
         STRDATA(3,JN) = SF1
         STRDATA(4,JN) = POS
         STRDATA(5,JN) = V
         STRDATA(6,JN) = QN(JN)
    6 CONTINUE
      RETURN 
  100   FORMAT(3(4X,F8.0),2(4X,I2),4X,F8.0)
  101   FORMAT(A1,8(/))
  102   FORMAT(59X,F8.3)
      END
