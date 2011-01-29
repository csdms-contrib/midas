C**************************SUBROUTINE WRITERSTRT***********************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: WRITES THE .e FILE ALL THE PARAMETERS NECESSARY TO
C                   DO A RESTART
C----------------------------------------------------------------------
      SUBROUTINE WRITERSTRT
      IMPLICIT NONE
      REAL*8 RHO,G
      REAL*8 VISKIN,THETAC                                              
      REAL*8 DIMID(20,3),SIGMA(3)                                      
      REAL*8 DIFEND(20,3),DIMIN(3)                                    
      REAL*8 PROPT(50)                                               
      REAL*8 Q,MANN(100),DINIT,FLOWEL                               
      REAL*8 THET50(100),VISC                                      
      REAL*8 COVAR,A,B,CON1,CON2                                  
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100), 
     #       AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 CONIM(20,3,100),CVELM(20,3,100),PHI                       
      REAL*8 SDBDX(20,3,100),SDBDT(20,3,100),SXFLUX(20,3,100),        
     #       STFLUX(20,3,100),STHICK(100)
      REAL*8 YMAX(100),DELY(100),VEL(30,100),EPSY(30,100),PSI,BETA   
      REAL*8 DENDM(60)                                              
      REAL*8 AWPM(20,3,100),D50M(100)                              
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100), 
     #       SVELUM(100)
      REAL*8 SBLWTA(20,3,100),SBLWT(20,3,100),SCVELIA(20,3,100),        
     #       SCVELI(20,3,100),SPRTHICK(20,3,100),TOTTHICK(20,3,100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                                   
      REAL*8 PLOTPRINT                                                
      REAL*8 WIDTOP,DELWID,ELEVTOP,SLOPE                             
      REAL*8 RLONG,TFILE                                            
      INTEGER*4 NTOT
      INTEGER*4 NK
      INTEGER*4 NEXT,ICONT
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*2 OPTSUSP,OPTLNODE
      INTEGER*2 PRINT(7)
      INTEGER*4 COUNTER
      INTEGER*4 NSIZE(3),NSIGMA,M,N                                
      INTEGER*4 OPTEL,OPTWID
      INTEGER*4 NT,OPTDIST,IFILE,OPTCONT,NTIM
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
      COMMON /COM6/PROPT,NK
      COMMON /COM11/Q,MANN,DINIT,ICONT,NEXT,FLOWEL
      COMMON /COM12/THET50,VISC
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #OPTBINP,OPTSINP,OPTERO
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #AVDENS
      COMMON /COM15/CONIM,CVELM,PHI,OPTSUSP,OPTLNODE
      COMMON /COM17/SDBDX,SDBDT,SXFLUX,STFLUX,STHICK
      COMMON /COM19/YMAX,DELY,VEL,EPSY,PSI,BETA
      COMMON /COM20/DENDM,NTOT
      COMMON /COM21/RUN,PRINT
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM27/SBLWTA,SBLWT,SCVELIA,SCVELI,SPRTHICK,
     #TOTTHICK,COUNTER
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
      COMMON /COM30/PLOTPRINT
      COMMON /COM31/OPTEL,OPTWID,WIDTOP,DELWID,ELEVTOP,SLOPE
      COMMON /COM33/ RLONG,TFILE,NTIM,NT,OPTDIST,
     #IFILE,OPTCONT,FOUT,RUNNAME
C-------------------------------------------------------------------------------
      INTEGER*4 REPR,I,J,K,JN
      CHARACTER*24  F,FF
C-------------------------------------------------------------------------------
      REPR = REPR + 1
      DO 1 K = 1, 24
         IF (FOUT(K:K) .EQ. ' ') GO TO 2
    1 CONTINUE
    2 CONTINUE
      F = CHAR(48+INT(REPR/10))//CHAR(48+REPR-INT(REPR/10)*10)
      FF = FOUT(1:K-1)//F
      F = FF
      OPEN(2, FILE='/users/vogel/data/'//F)
      WRITE (2, 100) N, NK, OPTURB, COVAR, RLONG, DELX, DELT/3600., MANN
      WRITE (2, 101) RHO, G, VISC, A, B, VOLCON, THETAC, ACTIVE
      WRITE (2, 102) NSIGMA, OPTDIST, OPTHIDE, OPTFUNC, OPTSALT, 
     #   OPTLNODE, CON1, CON2, PHI, PSI, BETA
      WRITE (2, 103) (SIGMA(J), J = 1, NSIGMA)
      WRITE (2, 104) (NSIZE(J), J = 1, NSIGMA)
      WRITE (2, 104) (PRINT(J), J = 1, 7)
      WRITE (2, 105) XPRINT, TPRINT/3600., TFILE/3600., PLOTPRINT
      WRITE (2, 106) NTIM
      WRITE (2, 107) OPTBINP, OPTSINP, OPTERO, OPTEL, OPTWID, OPTSUSP
      WRITE (2, 108) RUN
      WRITE (2, 109) RUNNAME
      WRITE (2, 110) (WIDM(JN), JN = 1, M + 1)
      WRITE (2, 110) (ELEVM(JN), JN = 1, M + 1)
      DO 5 JN = 1, M
         DO 4 J = 1, NSIGMA
            DO 3 I = 1, NSIZE(J)
               WRITE (2, 114) AWPM(I,J,JN), AWPOM(I,J,JN)
    3       CONTINUE
    4    CONTINUE
    5 CONTINUE
      WRITE (2, 106) NTOT
      WRITE (2, 111) (DENDM(K), K = 1, NTOT)
      WRITE (2, 112) DINIT, Q, OPTCONT, ICONT, NEXT
      IF (OPTCONT .EQ. 2) WRITE (2, 116) FLOWEL
      WRITE (2, 113) NT, NTOT, COUNTER
      WRITE (2, 103) (DIMIN(J), J = 1, NSIGMA)
      DO 8 J = 1, NSIGMA
         DO 7 I = 1, NSIZE(J)
            WRITE (2, 114) DIMID(I,J), DIFEND(I,J)
            DO 6 IFILE = 2, M
               WRITE (2, 115) SDBDX(I,J,IFILE), SBLWTA(I,J,IFILE), 
     #            SXFLUX(I,J,IFILE), SCVELIA(I,J,IFILE), SPRTHICK(I,J,
     #            IFILE), TOTTHICK(I,J,IFILE)
    6       CONTINUE
    7    CONTINUE
    8 CONTINUE
      WRITE (2, 111) (PROPT(K), K = 1, NK)
  100    FORMAT(1X,2I2,I1,F10.8,2F10.3,F10.0,F10.5)
  101    FORMAT(1X,F8.2,7F8.5)
  102    FORMAT(1X,6I2,5F8.3)
  103    FORMAT(1X,9F8.2)
  104    FORMAT(1X,9I2)
  105    FORMAT(1X,4F8.2)
  106    FORMAT(1X,I9)
  107    FORMAT(1X,5I1)
  108    FORMAT(1X,A72)
  109    FORMAT(1X,A4)
  110    FORMAT(1X,9F8.3/1X,9F8.3/1X,2F8.3)
  111    FORMAT(6(1X,9F8.6/),6F8.6)
  112   FORMAT(1X,2F8.3,I2,I5,I12)
  113    FORMAT(1X,4(I9))
  114    FORMAT(1X,2F8.6)
  115    FORMAT(1X,6E10.4)
  116    FORMAT(1X,F8.3)
C 117    FORMAT(1X,3F8.2)
      RETURN 
      END
