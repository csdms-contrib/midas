C**************************SUBROUTINE SUSP1****************************
C----------------------------------------------------------------------
C  THIS SUBROUTINE: CALCULATES SUSPENDED LOAD TRANSPORT RATES
C----------------------------------------------------------------------
      SUBROUTINE SUSP1()
      IMPLICIT NONE
      REAL*8 SETV(20,3),TC(20,3,100),SHVELC(20,3,100)                  
      REAL*8 BEDCON(20,3,100)                                         
      REAL*8 COVAR,A,B,CON1,CON2                                     
      REAL*8 CONIM(20,3,100),CVELM(20,3,100),PHI                    
      REAL*8 CONCM(20,3,30),CONC(30)                               
      REAL*8 YMAX(100),DELY(100),VEL(30,100),EPSY(30,100),PSI,BETA
      REAL*8 BLWTM(20,3,100),BLWTAM(20,3,100),SUM1(100),SUM2(100)  
      REAL*8 WIDM(100),ELEVM(100),DCHM(100),TOXM(100),SVELM(100), 
     #SVELUM(100)
      INTEGER*4 OPTURB,OPTHIDE,OPTSALT,OPTFUNC,OPTBINP,OPTSINP,OPTERO
      INTEGER*2 OPTSUSP,OPTLNODE
      INTEGER*4 NSIZE(3),NSIGMA,M,N                              
C----------------------------------------------------------------------
C     COMMON BLOCK DECLARATION
C----------------------------------------------------------------------
      COMMON /COM5/SETV,TC,SHVELC
      COMMON /COM10/BEDCON
      COMMON /COM13/COVAR,A,B,CON1,CON2,OPTURB,OPTHIDE,OPTSALT,OPTFUNC,
     #OPTBINP,OPTSINP,OPTERO
      COMMON /COM15/CONIM,CVELM,PHI,OPTSUSP,OPTLNODE
      COMMON /COM18/CONCM,CONC
      COMMON /COM19/YMAX,DELY,VEL,EPSY,PSI,BETA
      COMMON /COM23/BLWTAM,BLWTM,SUM1,SUM2
      COMMON /COM25/WIDM,ELEVM,DCHM,TOXM,SVELM,SVELUM
      COMMON /COM28/NSIZE,NSIGMA,M,N
C----------------------------------------------------------------------
C       DECLARATION STATEMENTS--EXTERNAL VARIABLES NOT IN COMMON BLOCKS
C----------------------------------------------------------------------
C     INTEGER*4 NT
C----------------------------------------------------------------------
      INTEGER*4 I,J,IV,L,NN
      REAL*8 CA,W,Z,AV,Y
      CHARACTER*1 DUMMY
C----------------------------------------------------------------------
C      THIS SUBROUTINE IS DESIGNED TO CALCULATE THE BED CONCENTRATION
C      FOR ALL SIZE-DENSITY FRACTIONS AT THE FIRST NODE (JN=1), THE
C      INPUT NODE. AS THIS NEVER CHANGES, IT ONLY NEEDS TO BE
C      CALCULATED AT TIMESTEP 0
C      IF OPTSINP=2 HOWEVER, THE EQUILIBRIUM CASE, INPUT AT NODE 1
C      IS DYNAMIC
C----------------------------------------------------------------------
      DO 2 J = 1, NSIGMA
         DO 1 I = 1, NSIZE(J)
            W = SETV(I,J)
            CONIM(I,J,1) = 0.0D0
            CVELM(I,J,1) = 0.0D0
            IF (OPTSINP .EQ. 1) BEDCON(I,J,1) = 0.0
C----------------------------------------------------------------------
C       OPTSINP=2 : SUSPENDED LOAD CONCENTRATIONS AT THE HEAD ARE
C       IN EQUILIBRIUM WITH THE FLOW CONDITIONS AT THE HEAD
C       THIS IS DEALT WITH IN SUSP
C----------------------------------------------------------------------
    1    CONTINUE
    2 CONTINUE
C----------------------------------------------------------------------
C      IF OPTSINP=3 A USER DEFINED BED- AND SUSPENDED LOAD IS INPUT
C      AT THE HEAD. VALUES OF BEDCON ARE READ IN AT THE HEAD AND
C      AN EQUILIBRIUM PROFILE IS THEN CALCULATED
C----------------------------------------------------------------------
      IF (OPTSINP .EQ. 3) THEN
         READ (4, 101) DUMMY
         DO 6 J = 1, NSIGMA
            READ (4, 100) (BEDCON(I,J,1), I = 1, NSIZE(J))
            DO 5 I = 1, NSIZE(J)
               CA = BEDCON(I,J,1)
               W = SETV(I,J)
               Z = W/(0.4*SVELUM(2))
               AV = YMAX(2)
               NN = N - 1
               DO 3 L = 1, NN
C----------------------------------------------------------------------
C     CALCULATE EQUILIBRIUM CONCENTRATIONS USING ROUSE EQUATIONS
C----------------------------------------------------------------------
                  Y = DCHM(2) - (L-1)*DELY(2)
                  IF (L .EQ. 1) Y = DCHM(2) - DELY(2)/4
                  CONC(L) = CA*(((DCHM(2)-Y)*AV)/(Y*(DCHM(2)-AV)))**Z
    3          CONTINUE
               CONC(N) = CA
               DO 4 IV = 1, N - 1
                  CONIM(I,J,1) = CONIM(I,J,1) + (CONC(IV)+CONC(IV+1))*
     #               DELY(2)/2.
                  CVELM(I,J,1) = CVELM(I,J,1) + (CONC(IV)*VEL(IV,2)+CONC
     #               (IV+1)*VEL(IV+1,2))*DELY(2)/2.
    4          CONTINUE
               CONIM(I,J,1) = CONIM(I,J,1) + (CONC(N)*YMAX(2))/2.
               BLWTM(I,J,1)=BLWTM(I,J,1)+(CONC(N)*VEL(N,2))*YMAX(2)/2.
    5       CONTINUE
    6    CONTINUE
      ENDIF
  100  FORMAT(9F8.0/9F8.0/2F8.0)
  101  FORMAT(A1,/)
      RETURN 
      END
