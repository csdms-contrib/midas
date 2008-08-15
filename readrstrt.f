C*********************SUBROUTINE READRSTRT****************************
C---------------------------------------------------------------------
C  THIS SUBROUTINE: READS THE .e FILE TO RESTART A RUN
C---------------------------------------------------------------------
      SUBROUTINE READRSTRT
      IMPLICIT NONE
      REAL*8 DIMID(20,3),SIGMA(3)                                    
      REAL*8 DIFEND(20,3),DIMIN(3)                
      REAL*8 PROPT(50)                           
      REAL*8 Q,MANN(100),DINIT,FLOWEL           
      REAL*8 COVAR,A,B,CON1,CON2               
      REAL*8 VOLCON,ACTIVE,PRTHICK(20,3,100),NWAVDENS(100),THICKM(100), 
     #       AWPOM(20,3,100),MSPRTHK(20,3,100),AVDENS(100)
      REAL*8 CONIM(20,3,100),CVELM(20,3,100),PHI                       
      REAL*8 SDBDX(20,3,100),SDBDT(20,3,100),SXFLUX(20,3,100),        
     #       STFLUX(20,3,100),STHICK(100)
      REAL*8 DENDM(60)                                               
      REAL*8 AWPM(20,3,100),D50M(100)                              
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100), 
     #       SVELUM(100)
      REAL*8 SBLWTA(20,3,100),SBLWT(20,3,100),SCVELIA(20,3,100), 
     #       SCVELI(20,3,100),SPRTHICK(20,3,100),TOTTHICK(20,3,100)
      REAL*8  PLOTPRINT                                        
      REAL*8 WIDTOP,DELWID,ELEVTOP,SLOPE                      
      REAL*8 RLONG,TFILE                                     
      INTEGER*4 NK
      INTEGER*4 NEXT,ICONT
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*2 OPTSUSP,OPTLNODE
      INTEGER*4 NTOT
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
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM6/PROPT,NK
      COMMON /COM11/Q,MANN,DINIT,ICONT,NEXT,FLOWEL
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #              OPTBINP,OPTSINP,OPTERO
      COMMON /COM14/VOLCON,ACTIVE,PRTHICK,NWAVDENS,THICKM,AWPOM,MSPRTHK,
     #              AVDENS
      COMMON /COM15/CONIM,CVELM,PHI,OPTSUSP,OPTLNODE
      COMMON /COM17/SDBDX,SDBDT,SXFLUX,STFLUX,STHICK
      COMMON /COM20/ NTOT,DENDM
      COMMON /COM21/RUN,PRINT
      COMMON /COM22/AWPM,D50M
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM27/SBLWTA,SBLWT,SCVELIA,SCVELI,COUNTER,SPRTHICK,
     #              TOTTHICK
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM30/PLOTPRINT
      COMMON /COM31/OPTEL,OPTWID,WIDTOP,DELWID,ELEVTOP,SLOPE
      COMMON /COM33/RLONG,TFILE,NTIM,NT,OPTDIST,
     #              IFILE,OPTCONT,FOUT,RUNNAME
C-------------------------------------------------------------------------------
      INTEGER*4 I,J,K,JN
C-------------------------------------------------------------------------------
      READ (4, 102) NTIM
      READ (4, 103) OPTBINP, OPTSINP, OPTERO, OPTEL, OPTWID, OPTSUSP
      READ (4, 104) RUN
      READ (4, 105) RUNNAME
      READ (4, 106) (WIDM(JN), JN = 1, M + 1)
      READ (4, 106) (ELEVM(JN), JN = 1, M + 1)
      DO 3 JN = 1, M
         DO 2 J = 1, NSIGMA
            DO 1 I = 1, NSIZE(J)
               READ (4, 110) AWPM(I,J,JN), AWPOM(I,J,JN)
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
      READ (4, 102) NTOT
      READ (4, 107) (DENDM(K), K = 1, NTOT)
      READ (4, 108) DINIT, Q, OPTCONT, ICONT, NEXT
      IF (OPTCONT .EQ. 2) READ (4, 112) FLOWEL
      READ (4, 109) NT, NTOT, COUNTER
      READ (4, 100) (DIMIN(J), J = 1, NSIGMA)
      DO 6 J = 1, NSIGMA
         DO 5 I = 1, NSIZE(J)
            READ (4, 110) DIMID(I,J), DIFEND(I,J)
            DO 4 IFILE = 2, M
               READ (4, 111) SDBDX(I,J,IFILE), SBLWTA(I,J,IFILE), SXFLUX
     #            (I,J,IFILE), SCVELIA(I,J,IFILE), SPRTHICK(I,J,IFILE), 
     #            TOTTHICK(I,J,IFILE)
    4       CONTINUE
    5    CONTINUE
    6 CONTINUE
      READ (4, 107) (PROPT(K), K = 1, NK)
  100    FORMAT(9F8.5)
  101    FORMAT(4F10.2)
  102    FORMAT(I9)
  103    FORMAT(5I1)
  104    FORMAT(A72)
  105    FORMAT(A4)
  106    FORMAT(9F8.3/1X,9F8.3/1X,2F8.3)
  107    FORMAT(6(1X,9F8.6/),6F8.6)
  108   FORMAT(2F8.3,I2,I5,I12)
  109    FORMAT(4(I9))
  110    FORMAT(2F8.6)
  111    FORMAT(6E10.4)
  112    FORMAT(F8.3)
  113    FORMAT(3F8.2)
      RETURN 
      END
