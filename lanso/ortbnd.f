C
C @(#)ortbnd.f	3.10 (BNP) 5/19/89; from ortbnd.f 2.7 10/17/87
C
      SUBROUTINE ORTBND(J,ALF,BET,ETA,OLDETA)
      INTEGER J
      DOUBLE PRECISION ALF(J),BET(J),ETA(J),OLDETA(J)
C
C.... Update the eta recurrence.
C
C.... inputs
C.... J      dimension of T
C.... ALF    diagonal elements of the tridiagonal T
C.... BET    off-diagonal elements of T
C.... ETA    orthogonality estimate of Lanczos vectors at step J
C.... OLDETA orthogonality estimate of Lanczos vectors at step J-1
C
C.... outputs
C.... ETA    orthogonality estimate of Lanczos vectors at step J+1
C.... OLDETA orthogonality estimate of Lanczos vectors at step J
C
C.... BLAS routines:    DSWAP
C
      DOUBLE PRECISION RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
C
      INTEGER I
      DOUBLE PRECISION T
C
      IF (J.LE.1) RETURN
      IF (RNM.EQ.0.0D0) THEN
         OLDETA(J-1) = 1.0D0
         GOTO 200
      ENDIF
      T = ANORM*RCEPS1
      IF (J.GT.2) THEN
         OLDETA(1) = (BET(2)*ETA(2)+(ALF(1)-ALF(J))*
     *      ETA(1)-BET(J)*OLDETA(1))
         OLDETA(1) = (OLDETA(1)+SIGN(T,OLDETA(1)))/RNM
      ENDIF
      DO 100 I = 2,J-2
         OLDETA(I) = BET(I+1)*ETA(I+1)+(ALF(I)-ALF(J))*ETA(I)+
     *      BET(I)*ETA(I-1)-BET(J)*OLDETA(I)
         OLDETA(I) = (OLDETA(I)+SIGN(T,OLDETA(I)))/RNM
100   CONTINUE
      OLDETA(J-1) = T/RNM
200   CALL DSWAP(J-1,OLDETA,1,ETA,1)
      ETA(J) = EPS1
      RETURN
      END
