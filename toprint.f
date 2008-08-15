C****************************SUBROUTINE TOPRINT************************
C  THIS SUBROUTINE: PRINTS ENTRAINMENT VALUES
C----------------------------------------------------------------------
      SUBROUTINE TOPRINT(NT)
      IMPLICIT NONE
      REAL*8 RHO,G                                                       
      REAL*8 VISKIN,THETAC                                              
      REAL*8 DIMID(20,3),SIGMA(3)                                      
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)                 
      REAL*8 PROPT(50)                                               
      REAL*8 SHVEL(50,100),TO(50,100)                               
      REAL*8 DMAX(50,3,100),DMIN(50,3,100)                         
      REAL*8 THET50(100),VISC                                     
      REAL*8 COVAR,A,B,CON1,CON2                                 
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100),
     #AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 AWPM(20,3,100),D50M(100)                                 
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100),    
     #SVELUM(100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                                
      INTEGER*4 NK
      INTEGER*4 IMIN(50,3,100),IMAX(50,3,100)
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*2 PRINT(7)
      INTEGER*4 NSIZE(3),NSIGMA,M,N                                
      CHARACTER*72 RUN                                            
C----------------------------------------------------------------------
C               COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM1/ RHO,G
      COMMON /COM2/VISKIN,THETAC
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM6/PROPT,NK
      COMMON /COM7/SHVEL,TO
      COMMON /COM8/DMAX,DMIN,IMIN,IMAX
      COMMON /COM12/THET50,VISC
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #OPTBINP,OPTSINP,OPTERO
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM21/RUN,PRINT
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--EXTERNAL VARIABLES NOT IN COMMON BLOCKS
C----------------------------------------------------------------------
      INTEGER*4 NT,JN
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 POS
      INTEGER*4 J
      INTEGER*4 I
      INTEGER*4 K,DUM
C----------------------------------------------------------------------
      IF (NT .EQ. 0) THEN
         WRITE (7, 100) 
  100     FORMAT(' ',7X, 'SEDIMENT TRANSPORT MODEL FOR HETEROGENEOUS SIZ
     #E-DENSITY MIXTURES'/,30X,'VERSION 4.0'/)
         WRITE (7, 101) RUN
  101     FORMAT(' ',A72//)
         WRITE (7, 111) 
         GO TO (1,2) OPTURB
    1    CONTINUE
         WRITE (7, 102) 
  102     FORMAT(' ','BED SHEAR STRESS IS GAUSSIAN DISTRIBUTED WITH')
         GO TO 3
    2    CONTINUE
         WRITE (7, 103) 
  103   FORMAT(' ','BED SHEAR STRESS IS LOG-GAUSSIAN DISTRIBUTED WITH')
    3    CONTINUE
         WRITE (7, 104) COVAR
  104    FORMAT(' ',' COEFFICIENT OF VARIATION',25X,F10.2)
         WRITE (7, 105) RHO, G, VISC, VISKIN, VOLCON, D50M(1), A, B, 
     #      THETAC, ACTIVE, DELT
  105     FORMAT(' ','FLUID DENSITY',37X,F10.1,'KG M-3'/1X,    'GRAVITAT
     #IONAL ACCELERATION',24X,F10.2,'M S-2'/1X,'MOLECULAR VISCOSITY',31
     #X,F10.5,'PA S'/1X,'KINEMATIC VISCOSITY',31X,F10.7,'M2 S-1'/1X, 'VO
     #LUME CONCENTRATION OF BED SEDIMENT',14X,F10.2/1X, 'MEDIAN GRAIN DI
     #AMETER OF AVAILABLE SEDIMENT',7X,F10.7,'M'/1X, 'A IN BEDLOAD EQUAT
     #ION',29X,F10.2/1X,'B IN SUSPENSION EQUATION',26X,F10.2/1X, 'DIMENS
     #IONLESS BED SHEAR STRESS FOR ROUGH BOUNDARY',1X,F10.3/1X, 'ACTIVE
     #LAYER THICKNESS',28X,F10.3,'M'/1X,'TIME INCREMENT STEP (DELT)',23X
     #,F10.1)
         IF (COVAR.GT.0.35 .AND. OPTURB.EQ.1) THEN
            WRITE (7, 106) 
  106    FORMAT(' ','WARNING!CHOSEN COEFFICIENT OF VARIATION CAUSES NEGA
     #TIVE SHEAR STRESSES'//' IN A NORMAL DISTRIBUTION.A LOGNORMAL DISTR
     #IBUTION SHOULD BE USED')
         ENDIF
         WRITE (7, 111) 
      ENDIF
C----------------------------------------------------------------------
C     WRITE OUT SHEAR STRESS VALUES AND THEIR DURATION, PLUS ARRAY
C     INDICES FOR EACH MINERAL SPECIES
C----------------------------------------------------------------------
      DO 8 JN = 2, M
         POS = (JN-2)*DELX/XPRINT + 0.000001
         IF (FLOAT(POS) - FLOAT(INT(POS)) .LE. 0.00001) THEN
            WRITE (7, 112) 
            WRITE (7, 107) (JN - 2)*DELX, NT, TOXM(JN), D50M(JN)
            WRITE (7, 112) 
  107 FORMAT(' '//10X,'METRES ALONG THE REACH = ',F10.4,' M',
     #/24X,'TIMESTEP = ',I8,/11X,'MEAN BED SHEAR STRESS = ',F10.2,' PA'
     #,/11X,'MEDIAN GRAIN DIAMETER = ',F10.8/)
            DO 7 J = 1, NSIGMA
               WRITE (7, 108) J, SIGMA(J)
  108   FORMAT(' ',//,27X,'MINERAL SPECIES #',I2,/20X,'MINERAL DENSITY'
     #,F10.1,'KG M-3',//1X,'BED SHEAR STRESS  MINIMUM  ARRAY   MAXIMUM A
     #RRAY',22X,'GRAIN SIZE AV-BLE   CRITICAL BED    CRITICAL    SETTLIN
     #G  ARRAY'/
     #3X,'MIDPOINT (PA)   SIZE(M)  INDEX   SIZE(M) INDEX  DURATION',10X,
     #'  MIDPOINT  WT. PROP. SHEAR STRESS SHEAR VELOCITY VELOCITY NUMBER
     #'/
     #70X,'    (M)                  (PA)        (M S-1)      (M S-1)'/)
               DUM = MIN(NK,NSIZE(J))
               DO 4 K = 1, DUM
                  I = K
                  WRITE (7, 109) TO(K,JN), DMIN(K,J,JN), IMIN(K,J,JN), 
     #               DMAX(K,J,JN), IMAX(K,J,JN), PROPT(K), DIMID(I,J), 
     #               AWPM(I,J,JN), TC(I,J,JN), SHVELC(I,J,JN), SETV(I,J)
     #               , I
    4          CONTINUE
               DO 5 K = DUM + 1, NK
                  WRITE (7, 109) TO(K,JN), DMIN(K,J,JN), IMIN(K,J,JN), 
     #               DMAX(K,J,JN), IMAX(K,J,JN), PROPT(K)
    5          CONTINUE
               DO 6 I = DUM + 1, NSIZE(J)
                  WRITE (7, 109) DIMID(I,J), AWPM(I,J,JN), TC(I,J,JN), 
     #               SHVELC(I,J,JN), SETV(I,J), I
    6          CONTINUE
  109    FORMAT(F10.2,7X,F10.7,I4,2X,F10.7,I4,2X,F10.7,10X,
     #2F10.6,F9.3,F13.3,6X,F9.3,2X,I3)
  110    FORMAT(85X,2F10.6,F9.3,F13.3,6X,F9.3,2X,I3)
    7       CONTINUE
            WRITE (7, 112) 
  111 FORMAT('','-------------------------------------------------------
     #---------')
  112 FORMAT('','-------------------------------------------------------
     #------------------------------------------------------------------
     #-----------')
         ENDIF
    8 CONTINUE
      RETURN 
      END
