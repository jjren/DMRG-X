C
C @(#)purge.f	3.16 (BNP) 5/11/89; from purge.f 2.13 6/25/88
C
      SUBROUTINE PURGE(N,LL,J,R,Q,RA,QA,WRK,ETA,OLDETA,MSGLVL)
      INTEGER N,LL,J,MSGLVL
      DOUBLE PRECISION R(N),Q(N),RA(N),QA(N),WRK(N),ETA(J),OLDETA(J)
C
C.... This routine examines ETA to decide whether
C.... re-orthogonalization should be performed.
C
C.... N      dimension of the eigenproblem
C.... LL     no. of initial Lanczos vectors in local orthog.
C.... J      current Lanczos step
C.... R      the residual vector to become the next Lanczos vector
C.... Q      the current Lanczos vector
C.... RA     the product of the mass matrix and r
C.... QA     the product of the mass matrix and q
C.... WRK    a temporary vector to hold the previous Lanczos vectors
C.... ETA    state of orthogonality between r and previous Lanczos vectors
C.... OLDETA state of orthogonality between q and previous Lanczos vectors
C
C.... BLAS routines:    DAXPY,DDOT,IDAMAX
C.... subroutines:      none
C.... user-supplied:    OPM,STORE
C
      INTEGER RETRQ
      PARAMETER (RETRQ = 2)
C
      DOUBLE PRECISION RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      INTEGER K,I,LOOP,IDAMAX
      DOUBLE PRECISION T,TQ,TR,REPS1,DDOT
C
      IF (J.LE.LL+1) RETURN
      K = IDAMAX(J-(LL+1),ETA(LL+1),1)+LL
C
      IF (ABS(ETA(K)).GT.REPS) THEN
         REPS1 = EPS1/REPS
         DO 55 LOOP = 1,2
            IF (RNM.GT.TOL) THEN
C
C....       Bring in a Lanczos vector t and orthogonalize both r and q
C....       against it
C
               TQ = 0.0D0
               TR = 0.0D0
               DO 50 I = 1,J-1
                  CALL STORE(N,RETRQ,I,WRK)
                  T = -DDOT(N,QA,1,WRK,1)
                  TQ = TQ+ABS(T)
                  CALL DAXPY(N,T,WRK,1,Q,1)
                  T = -DDOT(N,RA,1,WRK,1)
                  TR = TR+ABS(T)
                  CALL DAXPY(N,T,WRK,1,R,1)
50             CONTINUE
               CALL LANCZOS_OPM(N,Q,QA)
C
C....          restore local orthogonality
C
               T = -DDOT(N,R,1,QA,1)
               TR = TR+ABS(T)
               CALL DAXPY(N,T,Q,1,R,1)
C
               CALL LANCZOS_OPM(N,R,RA)
               RNM = SQRT(DDOT(N,RA,1,R,1))
               IF (TQ.LE.REPS1.AND.TR.LE.REPS1*RNM) GOTO 58
            ENDIF
55       CONTINUE
58       DO 60 I = LL+1,J
            ETA(I) = EPS1
            OLDETA(I) = EPS1
60       CONTINUE
      ENDIF
      RETURN
      END
