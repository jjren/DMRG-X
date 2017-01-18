C
C @(#)random.f	3.4 (BNP) 10/30/90; from random.f 2.2 10/13/87
C
      DOUBLE PRECISION FUNCTION RANDOM(IY)
      INTEGER IY
C
C        random is a uniform random number generator based  on  theory  and
C     suggestions  given  in  d.e.knuth (1969),  vol  2.  the integer  iy
C     should be initialized to an arbitrary integer prior to the first call
C     to random.the calling program should  not  alter  the  value  of  iy
C     between subsequent calls to random.values of random will be returned
C     in the interval (0,1).
C
      SAVE IA,IC,MIC,S,M2
C
      INTEGER IA,IC,MIC
      DOUBLE PRECISION S
C
      INTEGER M2
      DATA M2/0/
C
      IF (M2.EQ.0) THEN
C
C        if first entry, compute machine integer word length
C
         MIC = 1
10       M2 = MIC
         MIC = M2+M2
         IF (MIC.GT.M2) GOTO 10
         S = DBLE(M2)
C
C        compute multiplier and increment for linear congruential method
C
         IA = 8*INT(S*ATAN(1.D0)/8.D0)+5
         IC = 2*INT(S*(0.5D0-SQRT(3.D0)/6.D0))+1
         MIC = (M2-IC)+M2
C
C        s is the scale factor for converting to floating point
C
         S = 0.5D0/S
      ENDIF
C
C        compute next random number
C
      IY = IY*IA
C
C     the following statement is for computers which do not allow
C     integer overflow on addition
C
      IF (IY.GT.MIC) IY = (IY-M2)-M2
C
      IY = IY+IC
C
C     the following statement is for computers where the
C     word length for addition is greater than for multiplication
C
      IF (IY/2.GT.M2) IY = (IY-M2)-M2
C
C     the following statement is for computers where integer
C     overflow affects the sign bit
C
      IF (IY.LT.0) IY = (IY+M2)+M2
      RANDOM = DBLE(IY)*S
      RETURN
      END
