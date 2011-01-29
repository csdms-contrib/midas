C***************************SUBROUTINE DFIFTY**************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: 1) CALCULATES THE D50 OF THE BED MATERIAL
C
C  CREATORS: A. VAN NIEKERK & K. VOGEL
C----------------------------------------------------------------------
      SUBROUTINE DFIFTY
      IMPLICIT NONE
      REAL*8 DIMID(20,3),SIGMA(3)           
      REAL*8 DIFEND(20,3),DIMIN(3)         
      REAL*8 DENDM(60)                    
      REAL*8 AWPM(20,3,100),D50M(100)    
      INTEGER*4 NTOT
      INTEGER*4 NSIZE(3),NSIGMA,M,N     
C----------------------------------------------------------------------
C               COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM3/DIMID,SIGMA
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM20/DENDM,NTOT
      COMMON /COM22/AWPM,D50M
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C               INTERNAL VARIABLES
C----------------------------------------------------------------------
      REAL*8 VOLTOT(60),VOL,SUM,LOWLIM,LOWDEN
      INTEGER*4 J,I,K,JN
C----------------------------------------------------------------------
C     CALCULATE THE TOTAL VOLUME PRESENT IN EACH SEQUENTIAL SIZE
C     INTERVAL
C----------------------------------------------------------------------
      DO 9 JN = 2, M
         DO 1 K = 1, NTOT
            VOLTOT(K) = 0.0D0
    1    CONTINUE
         DO 5 J = 1, NSIGMA
            LOWLIM = DIMIN(J)
            LOWDEN = DIMIN(J)
            K = 1
            DO WHILE(DENDM(K) .LE. DIMIN(J))
               K = K + 1
            END DO
C----------------------------------------------------------------------
C     CALCULATE THE TOTAL VOLUME BELONGING TO A CERTAIN DENDM(K)
C     SIZE INTERVAL
C----------------------------------------------------------------------
            DO 4 I = 1, NSIZE(J)
               DO WHILE(DENDM(K).LE.DIFEND(I,J) .AND. K.LE.NTOT)
                  VOLTOT(K) = VOLTOT(K) + AWPM(I,J,JN)/SIGMA(J)*(DENDM(K
     #               )-LOWDEN)/(DIFEND(I,J)-LOWLIM)
                  LOWDEN = DENDM(K)
                  K = K + 1
               END DO
               LOWLIM = DIFEND(I,J)
    4       CONTINUE
    5    CONTINUE
C----------------------------------------------------------------------
C  CALCULATE D50
C----------------------------------------------------------------------
         VOL = 0.0D0
         DO 6 K = 1, NTOT
            VOL = VOL + VOLTOT(K)
    6    CONTINUE
         SUM = 0.0D0
         DO 7 K = 1, NTOT
            SUM = SUM + VOLTOT(K)
            IF (SUM .GE. 0.5*VOL) THEN
C----------------------------------------------------------------------
C  LINEAR INTERPOLATION BETWEEN SIZES WHOSE PERCENTILES ARE
C  RESPECTIVELY JUST LESS THAN AND JUST GREATER THAN 50%
C----------------------------------------------------------------------
               IF (K .EQ. 1) THEN
                  D50M(JN) = DENDM(1)
               ELSE
                  D50M(JN) = DENDM(K-1) + (DENDM(K)-DENDM(K-1))/VOLTOT(K
     #               )*(0.5D0*VOL-(SUM-VOLTOT(K)))
               ENDIF
               GO TO 8
            ENDIF
    7    CONTINUE
    8    CONTINUE
    9 CONTINUE
      D50M(M+1) = D50M(M)
      D50M(1) = D50M(2)
      RETURN 
      END
