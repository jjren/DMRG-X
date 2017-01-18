C
C @(#)startv.f	1.10 (BNP) 10/30/90
C
      DOUBLE PRECISION FUNCTION STARTV(N,J,R,WRK,NQ,EPS,MSGLVL)
      INTEGER N,J,NQ(4),MSGLVL
      DOUBLE PRECISION R(5*N),WRK(N),EPS
C
C.... This routine delivers a starting vector in R and returns |R|;
C.... it returns ZERO if range is spanned or if no starting vector
C.... within range of operator can be found.
C
C.... N      dimension of the eigenproblem
C.... J      starting index for a Lanczos run
C.... R      an array containing [r(j),q(j),q(j-1),p(j),p(j-1)/Mr(j)]
C.... NQ(4)  location pointers for the array R
C
C.... BLAS routines:    DAXPY,DDOT
C.... subroutines:      RANDOM
C.... user-supplied:    OP,OPM,STORE
C
      INTEGER RETRQ
      PARAMETER (RETRQ = 2)
C
      INTEGER IRAND,I,ID,LOOP
      DOUBLE PRECISION RNM2,T,RANDOM,DDOT
C
      DOUBLE PRECISION ZERO
      DATA ZERO/0.0D0/
C
C.... get initial vector, default is random
C
      RNM2 = DDOT(N,R,1,R,1)
      IRAND = 918272+J
      DO 60 ID = 1,3
         IF (ID.GT.1.OR.J.GT.1.OR.RNM2.EQ.ZERO) THEN
            DO 20 I = 1,N
               R(I) = RANDOM(IRAND)-0.5D0
20          CONTINUE
         ENDIF
CKW     make sure R is in the range of OPM
         CALL LANCZOS_OPM(N,R,R(NQ(3)))
         CALL LANCZOS_OPM(N,R(NQ(3)),R)
         CALL LANCZOS_OPM(N,R,R(NQ(3)))
         RNM2 = DDOT(N,R,1,R(NQ(3)),1)
         IF (RNM2.EQ.ZERO) THEN
            GOTO 60
         ELSE IF (J.GT.1) THEN
            DO 50 LOOP = 1,2
               DO 40 I = 1,J-1
                  CALL STORE(N,RETRQ,I,WRK)
                  T = -DDOT(N,R(NQ(3)),1,WRK,1)
                  CALL DAXPY(N,T,WRK,1,R,1)
40             CONTINUE
               CALL LANCZOS_OPM(N,R,R(NQ(3)))
               T = DDOT(N,R(NQ(3)),1,R,1)
               IF (T.LE.EPS*RNM2) THEN
                  T = ZERO
                  GOTO 55
               ENDIF
50          CONTINUE
55          RNM2 = T
         ENDIF
         IF (RNM2.GT.ZERO) GOTO 80
60    CONTINUE
80    STARTV = SQRT(RNM2)
      RETURN
      END
