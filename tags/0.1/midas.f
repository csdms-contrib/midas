C***************************PROGRAM MAIN***************************
C  THIS PROGRAM IS A NONUNIFORM,QUASI-UNSTEADY,MOVABLE BED,SINGLE
C  CHANNEL FLOW MODEL FOR HETEROGENEOUS SIZE-DENSITY MIXTURES
C
C  CREATORS: J. S. BRIDGE, R. SLINGERLAND, K. VOGEL, AND
C            A. VAN NIEKERK
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C     VARIABLES IN THE COMMON BLOCKS
C----------------------------------------------------------------------
      REAL*8 RHO,G
      REAL*8 VISKIN,THETAC
      REAL*8 DIMID(20,3),SIGMA(3) 
      REAL*8 DIFEND(20,3),DIMIN(3)
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)     
      REAL*8 PROPT(50)                                   
      REAL*8 SHVEL(50,100),TO(50,100)                    
      REAL*8 DMAX(50,3,100),DMIN(50,3,100)              
      REAL*8 SLWTM(20,3,100),SLWTAM(20,3,100)          
      REAL*8 BEDCON(20,3,100)                         
      REAL*8 Q,MANN(100),DINIT,FLOWEL                
      REAL*8 THET50(100),VISC                       
      REAL*8 COVAR,A,B,CON1,CON2                   
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100),
     #       AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 CONIM(20,3,100),CVELM(20,3,100),PHI
      REAL*8 DBDX(20,3,100),XFLUX(20,3,100),TXFLUX(100),TDBDX(100)
      REAL*8 SDBDX(20,3,100),SDBDT(20,3,100),SXFLUX(20,3,100),   
     #       STFLUX(20,3,100),STHICK(100)
      REAL*8 CONCM(20,3,30),CONC(30)                            
      REAL*8 YMAX(100),DELY(100),VEL(30,100),EPSY(30,100),PSI,BETA 
      REAL*8 DENDM(60)                                            
      REAL*8 AWPM(20,3,100),D50M(100)                           
      REAL*8 BLWTM(20,3,100),BLWTAM(20,3,100),SUM1(100),SUM2(100) 
      REAL*8 MEAN(3),STDEV(3),PROPTJ(3)                          
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),
     #       SVELUM(100)
      REAL*8 STRDATA(10,100)                                  
      REAL*8 SBLWTA(20,3,100),SBLWT(20,3,100),SCVELIA(20,3,100),
     #       SCVELI(20,3,100),SPRTHICK(20,3,100),TOTTHICK(20,3,100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                             
      REAL*8  PLOTPRINT                                         
      REAL*8 WIDTOP,DELWID,ELEVTOP,SLOPE                       
      REAL*8 ELEVOM(100)                                      
      REAL*8 RLONG,TFILE                                     
      REAL*8 QFAC                                          
      REAL*8 BLWTO(20,3)                                  
      INTEGER*4 NK
      INTEGER*4 IMIN(50,3,100),IMAX(50,3,100)
      INTEGER*4 NEXT,ICONT
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*2 OPTSUSP,OPTLNODE
      INTEGER*4 NTOT
      INTEGER*2 PRINT(7)
      INTEGER*4 COUNTER
      INTEGER*4 NSIZE(3),NSIGMA,M,N                               
      INTEGER*4 OPTEL,OPTWID
      INTEGER*4 NT,OPTDIST(3),IFILE,OPTCONT,NTIM
      INTEGER*4 QTMOD,THALFQ,QCOUNT,QSPMOD
      INTEGER*4 QSMOD,BCOUNT,THALFB,RENEW
      CHARACTER*72 RUN                                           
      CHARACTER*24  FOUT
      CHARACTER*4 RUNNAME
C----------------------------------------------------------------------
C     COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM2/VISKIN,THETAC
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM6/PROPT,NK
      COMMON /COM7/SHVEL,TO
      COMMON /COM8/DMAX,DMIN,IMIN,IMAX
      COMMON /COM9/SLWTM,SLWTAM
      COMMON /COM10/BEDCON
      COMMON /COM11/Q,MANN,DINIT,ICONT,NEXT,FLOWEL
      COMMON /COM12/THET50,VISC
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #        OPTBINP,OPTSINP,OPTERO
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #        AVDENS
      COMMON /COM15/CONIM,CVELM,PHI,OPTSUSP,OPTLNODE
      COMMON /COM16/DBDX,XFLUX,TDBDX,TXFLUX
      COMMON /COM17/SDBDX,SDBDT,SXFLUX,STFLUX,STHICK
      COMMON /COM18/CONCM,CONC
      COMMON /COM19/YMAX,DELY,VEL,EPSY,PSI,BETA
      COMMON /COM20/DENDM,NTOT
      COMMON /COM21/RUN,PRINT
      COMMON /COM22/AWPM,D50M
      COMMON /COM23/BLWTAM,BLWTM,SUM1,SUM2
      COMMON /COM24/MEAN,STDEV,PROPTJ
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM26/STRDATA
      COMMON /COM27/SBLWTA,SBLWT,SCVELIA,SCVELI,SPRTHICK,
     #        TOTTHICK,COUNTER
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
      COMMON /COM30/PLOTPRINT
      COMMON /COM31/OPTEL,OPTWID,WIDTOP,DELWID,ELEVTOP,SLOPE
      COMMON /COM32/ELEVOM
      COMMON /COM33/ RLONG,TFILE,NTIM,NT,OPTDIST,
     #        IFILE,OPTCONT,FOUT,RUNNAME
      COMMON /COM35/ QTMOD,THALFQ,QCOUNT,QSPMOD,QFAC
      COMMON /COM36/BLWTO,QSMOD,BCOUNT,THALFB,RENEW
C----------------------------------------------------------------------
C   INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 TIME,AWP(20,3),MANN1
      INTEGER*4 IRSTRT(20),NTSTAR
      INTEGER*4 NRUN,IRUN,K,I,J,NTFIN
      CHARACTER*24 F1(20),F2(20),F3(20)
      CHARACTER*4 MINSORT
      CHARACTER*8 FILENAME
C----------------------------------------------------------------------
      INTEGER*4 JN
C**********************************************************************
C----------------------------------------------------------------------
C     FILE MANAGEMENT PORTION
C----------------------------------------------------------------------
      NRUN = 1
      IRUN = 0
      OPEN(3, FILE='./files')
    1 CONTINUE
      READ (3, 100, END=2) F1(NRUN)
      READ (3, 100) F2(NRUN)
      IF (F1(NRUN) .EQ. F2(NRUN)) THEN
         IRSTRT(NRUN) = 1
         READ (3, 100) F2(NRUN)
      ENDIF
      READ (3, 100) F3(NRUN)
      NRUN = NRUN + 1
      GO TO 1
    2 CONTINUE
      NRUN = NRUN - 1
    3 CONTINUE
      IRUN = IRUN + 1
      OPEN(4, FILE='./'//F1(IRUN))
      OPEN(9, FILE='./'//F2(IRUN))
      DO 4 K = 1, 24
         IF (F3(IRUN)(K:K) .EQ. ' ') GO TO 5
    4 CONTINUE
    5 CONTINUE
      K = K - 1
      FOUT = F3(IRUN)(1:K)//'a'//CHAR(0)
      OPEN(8, FILE='./'//FOUT)
      FOUT = F3(IRUN)(1:K)//'b'//CHAR(0)
      OPEN(7, FILE='./'//FOUT)
      FOUT = F3(IRUN)(1:K)//'c'//CHAR(0)
      OPEN(12, FILE='./'//FOUT)
      FOUT = F3(IRUN)(1:K)//'d'//CHAR(0)
      OPEN(13, FILE='./'//FOUT)
      FOUT = F3(IRUN)(1:K)//'e'//CHAR(0)
C----------------------------------------------------------------------
C     READ INPUT PARAMETERS
C          >BED SHEAR STRESS DISTRIBUTION AND DERIVITIVES. IF OPTURB IS
C               1, THEN BED SHEAR STRESS IS GAUSSIAN DISTRIBUTED,
C               OTHERWISE IT IS LOG-GAUSSIAN DISTRIBUTED.
C          >SEDIMENT TRANSPORT PARAMETERS. IF OPTHIDE IS 1 THEN NO
C               GRAIN HIDING IS ASSUMED.
C          >HIDING FUNCTION OPTION. IF OPTFUNC IS 1 THEN THE KOMAR FUNC.
C               IS USED OTHERWISE THE EGAZIAROFF FUNCTION IS USED.
C          >SALTATION HEIGHT OPTION. IF OPTSALT IS 1 THEN THE BRIDGE
C               FUNCTION IS USED OTHERWISE THE EINSTEIN FUNCTION
C               IS USED.
C          >GRAIN SIZE AND DENSITY DISTRIBUTION CONTROL PARAMETERS
C          >GRAIN DENSITY VALUES
C          >WEIGHT PROPORTION OF AVAILABLE BED MATERIAL IN EACH
C               SIZE-DENSITY FRACTION. IF OPTDIST= 1, THEN WEIGHTS
C               ARE CALCULATED ASSUMING A LOG-GAUSSIAN DISTRIBUTION.
C          >INITIAL BOUNDARY CONDITION. IF OPT*INP IS 1 THEN CLEAR WATER
C               INFLOW IS ASSUMED. IF OPT*INP IS 2 THEN AN EQUILIBRIUM
C               CONDITION IS ASSUMED, IF OPT*INP IS 3 THEN A USER-
C               DEFINED CONCENTRATION PROFILE IS ASSUMED. IF OPTERO
C               IS 1 THEN EROSION OR DEPOSITION IS ALLOWED AT THE
C               HEAD OF THE REACH, IF OPTERO IS 2 THEN NO EROSION OR
C               DEPOSITION IS ALLOWED AT THE HEAD OF THE REACH.
C----------------------------------------------------------------------
      READ (4, 101) N, NK, OPTURB, COVAR
      READ (4, 102) RLONG, DELX, DELT
      READ (4, 103) RHO, G, VISC, MANN1
      READ (4, 104) A, B, VOLCON, THETAC, ACTIVE
      READ (4, 105) OPTHIDE, OPTFUNC, OPTSALT, OPTLNODE
      READ (4, 106) CON1, CON2, PHI, PSI, BETA
      READ (4, 107) OPTBINP, OPTSINP, OPTERO, OPTEL, OPTWID, OPTSUSP
      M = INT(RLONG/DELX+0.1) + 2
      IF (OPTEL .NE. 1) THEN
         READ (4, 108) ELEVTOP, SLOPE
      ELSE
         READ (4, 110) FILENAME
         OPEN(17, FILE='./'//FILENAME)
         READ (17, 111) (ELEVM(JN), JN = 2, M)
         READ (17, 112) ELEVM(M+1)
         CLOSE(17)
      ENDIF
      IF (OPTWID .NE. 1) THEN
         READ (4, 108) WIDTOP, DELWID
      ELSE
         READ (4, 110) FILENAME
         OPEN(17, FILE='./'//FILENAME)
         READ (17, 111) (WIDM(JN), JN = 2, M)
         CLOSE(17)
      ENDIF
      CALL ELWIDCALC
      READ (4, 113) NTIM
      READ (4, 114) QSPMOD, QFAC, QSMOD, THALFB, RENEW
      READ (4, 115) (PRINT(J), J = 1, 7)
      READ (4, 116) XPRINT, TPRINT, TFILE, PLOTPRINT
      READ (4, 117) NSIGMA
      DO 6 J = 1, NSIGMA
         READ (4, 118) OPTDIST(J), SIGMA(J), NSIZE(J)
    6 CONTINUE
       DELT=DELT*3600.
       TPRINT=TPRINT*3600.
      TFILE = TFILE*3600.
      VISKIN = VISC/RHO
C----------------------------------------------------------------------
C     IF RESTART FILE THEN READ IN VARIABLES & SKIP PRELIMINARIES
C----------------------------------------------------------------------
      IF (IRSTRT(IRUN) .EQ. 1) THEN
         CALL READRSTRT
         NTSTAR = NT + 1
         NTFIN = NT + NTIM
         GO TO 15
      ENDIF
C----------------------------------------------------------------------
C     IF OPTDIST = 1 THEN A PHI-NORMAL DISTRIBUTION IS GENERATED FOR
C     MINERAL J.  MEAN, DIMIN, AND STDEV MUST BE IN PHI UNITS
C----------------------------------------------------------------------
      READ (4, 119) FILENAME
      DO 10 J = 1, NSIGMA
         READ (4, 120) FILENAME
         IF (OPTDIST(J) .EQ. 1) THEN
            READ (4, 121) MEAN(J), STDEV(J), DIMIN(J), PROPTJ(J)
            CALL GSDIST (J)
         ELSE
C----------------------------------------------------------------------
C     IF OPTDIST= 2, THEN USE A USER-SPECIFIED SIZE DISTRIBUTION.
C     WITH ALL PARAMETERS TO BE READ IN.
C----------------------------------------------------------------------
            DO 7 I = 1, NSIZE(J)
               READ (4, 121) AWP(I,J), DIMID(I,J), DIFEND(I,J)
    7       CONTINUE
            READ (4, 121) DIMIN(J), PROPTJ(J)
            DO 9 JN = 1, M
               DO 8 I = 1, NSIZE(J)
                  AWPM(I,J,JN) = AWP(I,J)*PROPTJ(J)/100.
                  AWPOM(I,J,JN) = AWPM(I,J,JN)
    8          CONTINUE
    9       CONTINUE
         ENDIF
   10 CONTINUE
      READ (4, 122) RUN
      READ (4, 123) RUNNAME
      OPEN(18, FILE='./plotd.'//RUNNAME)
      OPEN(20, FILE='./dch.'//RUNNAME)
      OPEN(21, FILE='./trans.'//RUNNAME)
      OPEN(22, FILE='./gsdtr.'//RUNNAME)
      DO 11 J = 0, NSIGMA - 1
         MINSORT = 'min'//CHAR(48+J+1)
         OPEN(13 + J, FILE='./'//MINSORT//'.'//RUNNAME
     #      )
   11 CONTINUE
C----------------------------------------------------------------------
C     DEFINE D50 OF AVAILABLE SEDIMENT
C----------------------------------------------------------------------
      CALL DISCR
      CALL DFIFTY
C----------------------------------------------------------------------
C     CALCULATE VALUES OF CRITICAL BED SHEAR STRESS, CRITICAL
C     SHEAR VELOCITY, AND SETTLING VELOCITY FOR EACH GRAIN FRACTION
C----------------------------------------------------------------------
      CALL ENTRAIN
      CALL SETTLE
C----------------------------------------------------------------------
C            CALCULATE AVDENS-- THIS IS THE MASS/UNIT VOLUME IN THE
C            ORIGINAL ACTIVE LAYER
C---------------------------------------------------------------------
      DO 14 JN = 1, M
         AVDENS(JN) = 0.0
         DO 13 J = 1, NSIGMA
            DO 12 I = 1, NSIZE(J)
               AVDENS(JN) = AVDENS(JN) + AWPOM(I,J,JN)*SIGMA(J)
   12       CONTINUE
   13    CONTINUE
   14 CONTINUE
C----------------------------------------------------------------------
C     CALCULATE TIME-LOOP PARAMETERS
C----------------------------------------------------------------------
      NTSTAR = 0
      NTFIN = NTIM
C----------------------------------------------------------------------
C     BEGIN COMPUTATIONS OF BEDLOAD TRANSPORT AND BED EROSION AND
C     DEPOSITION, STARTING AT TIMESTEP NT = 0
C----------------------------------------------------------------------
   15 CONTINUE
      DO 20 NT = NTSTAR, NTFIN
         OPEN(19, FILE='./tst')
         WRITE (19, 124) NT
         CLOSE(19)
C----------------------------------------------------------------------
C     CLEAR ALL TRANSPORT ARRAYS
C----------------------------------------------------------------------
         DO 16 JN = 1, M
            MANN(JN) = MANN1
   16    CONTINUE
C----------------------------------------------------------------------
C  DETERMINE ALL THE NECESSARY CHANNEL DATA, OR READ VALUES. THIS
C  SUBROUTINE WILL BE CHANGED, DEPENDING ON THE SITUATION BEING MODELLED
C----------------------------------------------------------------------
         CALL FLDTA
         CALL SVELA
C----------------------------------------------------------------------
C  CALCULATE MASS FLUXES FOR SUSPENSION AND BEDLOAD. IF MORE OF A SIZE-
C  FRACTION IS ERODED THAN IS AVAILABLE IN THE BED, THEN
C----------------------------------------------------------------------
         TIME = DELT*NT/TPRINT + 0.00001
         IF (TIME - FLOAT(INT(TIME)) .LE. 0.0001) THEN
            IF (PRINT(1) .EQ. 1) CALL WRFLDTA (NT)
            IF (PRINT(7) .EQ. 1) CALL PLOTFILES (NT)
         ENDIF
C----------------------------------------------------------------------
C       BEGIN COMPUTATION OF BEDLOAD TRANSPORT
C       CALCULATE BED SHEAR STRESS AT EVERY NODE
C----------------------------------------------------------------------
         CALL TURB
C----------------------------------------------------------------------
C     CALCULATE BEDLOAD TRANSPORT RATES FOR EACH GRAIN FRACTION
C     AT EVERY NODE
C----------------------------------------------------------------------
         GO TO (17,18) OPTHIDE
   17    CONTINUE
         CALL DRANGE
         GO TO 19
   18    CONTINUE
         CALL ENTRAINH
         CALL DRANGEH
   19    CONTINUE
         CALL BEDLOAD (NT)
         TIME = TIME/PRINT(6)
         IF (TIME - FLOAT(INT(TIME)) .LE. 0.0001) THEN
            IF (PRINT(2) .EQ. 1) CALL TOPRINT (NT)
         ENDIF
C----------------------------------------------------------------------
C  BEGIN COMPUTATION OF SUSPENDED LOAD CONCENTRATION PROFILES
C  (CALCULATE OR GET THE BED LAYER THICKNESS, DETERMINE DELY)
C----------------------------------------------------------------------
C  GENERATE LOGARITHMIC VELOCITY DISTRIBUTION
C----------------------------------------------------------------------
         CALL LOGDIST
C----------------------------------------------------------------------
C  COMPUTE VERTICAL DIFFUSIVITIY PROFILE
C----------------------------------------------------------------------
         CALL EXVERT
         IF (NT .EQ. 0) CALL SUSP1
         CALL SUSP
C----------------------------------------------------------------------
C  CALCULATE EROSION/DEPOSITION AND NEW WEIGHT PROPORTIONS IN ACTIVE LAYER
C----------------------------------------------------------------------
         CALL ERODEP
C----------------------------------------------------------------------
C     DEFINE D50 OF AVAILABLE SEDIMENT
C----------------------------------------------------------------------
         CALL DFIFTY
         CALL TABLE (NT)
C----------------------------------------------------------------------
C  WRITE OUT RESULTS TO THE .c AND .d FILES
C----------------------------------------------------------------------
         IF (TIME - FLOAT(INT(TIME)) .LE. 0.0001) THEN
            IF (PRINT(3).EQ.1 .OR. PRINT(4).EQ.1) CALL RESULT (NT)
         ENDIF
C----------------------------------------------------------------------
C  PRINT RESULTS FOR RESTART FILE
C----------------------------------------------------------------------
         IF (PRINT(5) .EQ. 1) THEN
            IF (ABS(NT*DELT/TFILE-INT(NT*DELT/TFILE)) .LT. 0.000001) 
     #         THEN
               IF (NT .GT. 0) CALL WRITERSTRT
            ENDIF
         ENDIF
   20 CONTINUE
      IF (IRUN .LT. NRUN) GO TO 3
      STOP 
  100  FORMAT(A24)
  101  FORMAT(3(59X,I2,/),59X,F6.4)
  102  FORMAT(2(59X,F10.0,/),59X,F10.0)
  103  FORMAT(59X,F10.0,/,59X,F5.3,/,59X,F8.6/,59X,F8.6)
  104  FORMAT(59X,F6.3,4(/,59X,F8.6))
  105  FORMAT(3(59X,I1,/),59X,I1)
  106  FORMAT(59X,F6.3,4(/,59X,F7.4))
  107  FORMAT(5(59X,I1,/),59X,I1)
  108  FORMAT(/,/,/,59X,F8.3,/,59X,F8.6)
C 109  FORMAT(59X,F10.5,/,59X,F8.6)
  110  FORMAT(/,/,/,59X,A8,/)
  111  FORMAT(/,/,99F8.3)
  112  FORMAT(F8.0)
  113  FORMAT(59X,I10)
  114  FORMAT(2(59X,I1,/,59X,F8.4,/),59X,F8.4)
  115  FORMAT(/,6(59X,I1,/),59X,I1)
  116  FORMAT(3(59X,F10.2,/),59X,F10.2)
  117  FORMAT(59X,I1,4(/))
  118  FORMAT(15X,I1,14X,F10.2,5X,I2)
  119  FORMAT(A8,13(/))
  120  FORMAT(A8,/)
  121  FORMAT(3(6X,F10.0),6X,F10.0)
  122  FORMAT(/,1X,A72)
  123  FORMAT(59X,A4,/,/,/,/)
  124  FORMAT(1X,'TIMESTEP   ',I9)
      END
