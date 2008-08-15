
C***********************SUBROUTINE DISCR*******************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: DISCRETISES THE PARTICLE SIZE RANGES THROUGH
C                   ALL SIZE RANGES OF ALL THE COMBINED DENSITY
C                   FRACTIONS FOR CALCULATION OF THE NEW D50'S
C
C  CREATORS: K. VOGEL & P. B. FLEMINGS
C
C----------------------------------------------------------------------
      SUBROUTINE DISCR
      IMPLICIT NONE
      REAL*8 DIFEND(20,3),DIMIN(3)
      REAL*8 DENDM(60)           
      INTEGER*4 NTOT
      INTEGER*4 NSIZE(3),NSIGMA,M,N  
C----------------------------------------------------------------------
C       COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM4/DIFEND,DIMIN
      COMMON /COM20/DENDM,NTOT
      COMMON /COM28/NSIZE,NSIGMA,M,N
C-------------------------------------------------------------------
C       INTERNAL VARIABLES
C-------------------------------------------------------------------
      REAL*8 TEMP
      INTEGER*4 J,I,K1,K2,NRP,K
C----------------------------------------------------------------------
      NTOT = 0
      DO 2 J = 1, NSIGMA
         DO 1 I = 1, NSIZE(J)
            NTOT = NTOT + 1
            DENDM(NTOT) = DIFEND(I,J)
    1    CONTINUE
    2 CONTINUE
C-------------------------------------------------------------------
C         DO A BUBBLE SORT TO HAVE DENDM IN ASCENDING ORDER.
C-------------------------------------------------------------------
      DO 4 K1 = 1, NTOT - 1
         DO 3 K2 = K1 + 1, NTOT
            IF (DENDM(K1) .GT. DENDM(K2)) THEN
               TEMP = DENDM(K1)
               DENDM(K1) = DENDM(K2)
               DENDM(K2) = TEMP
            ENDIF
    3    CONTINUE
    4 CONTINUE
C-------------------------------------------------------------------
C          REMOVE REPEAT SIZE FRACTIONS
C-------------------------------------------------------------------
      NRP = 0
      DO 6 K = 1, NTOT - 1
         IF (DENDM(K-NRP+1) - DENDM(K-NRP) .LT. 0.1E-7) THEN
            NTOT = NTOT - 1
            NRP = NRP + 1
            DO 5 K1 = K - NRP + 1, NTOT
               DENDM(K1) = DENDM(K1+1)
    5       CONTINUE
            DENDM(NTOT+1) = 0.0D0
         ENDIF
    6 CONTINUE
      RETURN 
      END
