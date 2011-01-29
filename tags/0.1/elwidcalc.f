C**************************ELWIDCALC.F**********************************
C----------------------------------------------------------------------
C THIS SUBROUTINE: CALCULATES THE ELEVATIONS AND WIDTHS OF THE NODE
C                  POINTS FOLLOWING A USER SPECIFIED ALGORITHM
C----------------------------------------------------------------------
      SUBROUTINE ELWIDCALC
      IMPLICIT NONE
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100), 
     #       SVELUM(100)
      REAL*8 DELX,DELT,XPRINT,TPRINT                            
      REAL*8 WIDTOP,DELWID,ELEVTOP,SLOPE                       
      REAL*8 ELEVOM(100)                                      
      INTEGER*4 NSIZE(3),NSIGMA,M,N                              
      INTEGER*4 OPTEL,OPTWID
C----------------------------------------------------------------------
C     COMMON BLOCKS
C----------------------------------------------------------------------
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM28/NSIZE,NSIGMA,M,N
      COMMON /COM29/DELX,DELT,XPRINT,TPRINT
      COMMON /COM31/OPTEL,OPTWID,WIDTOP,DELWID,ELEVTOP,SLOPE
      COMMON /COM32/ELEVOM
C----------------------------------------------------------------------
C     DECLARATION STATEMENTS--INTERNAL VARIABLES
C---------------------------------------------------------------------
      INTEGER*4 JN
C----------------------------------------------------------------------
      DO 3 JN = M, 1, -1
         IF (OPTWID .EQ. 1) THEN
            GO TO 1
         ELSE IF (OPTWID .EQ. 2) THEN
            WIDM(JN) = WIDTOP + DELWID*(JN-2)*DELX
         ELSE IF (OPTWID .EQ. 3) THEN
            WIDM(JN) = WIDTOP*EXP((JN-2)*DELX*DELWID)
         ELSE IF (OPTWID .EQ. 4) THEN
            WIDM(JN) = WIDTOP + LOG((JN-2)*DELX*DELWID)
         ELSE IF (OPTWID .EQ. 5) THEN
            WIDM(JN) = WIDTOP*(FLOAT(JN-2)/FLOAT(M-2))**DELWID
            IF (JN .EQ. 2) WIDM(2) = WIDM(3)
         ENDIF
    1    CONTINUE
         IF (OPTEL .EQ. 1) THEN
            GO TO 2
         ELSE IF (OPTEL .EQ. 2) THEN
            ELEVM(JN) = ELEVTOP - SLOPE*(JN-2)*DELX
         ELSE IF (OPTEL .EQ. 3) THEN
            ELEVM(JN) = ELEVTOP*EXP((JN-2)*DELX*SLOPE)
         ELSE IF (OPTEL .EQ. 4) THEN
            ELEVM(JN) = ELEVTOP - LOG((JN-2)*DELX*SLOPE)
         ENDIF
    2    CONTINUE
    3 CONTINUE
      WIDM(1) = WIDM(2) - (WIDM(3)-WIDM(2))
      WIDM(M+1) = WIDM(M)
      DO 4 JN = 2, M
         ELEVOM(JN) = ELEVM(JN)
    4 CONTINUE
      RETURN 
      END
