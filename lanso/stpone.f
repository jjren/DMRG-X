C
C @(#)stpone.f	3.13 (BNP) 6/3/89; from stpone.f 2.8 6/17/88
C
      SUBROUTINE STPONE(N,R,WRK,ALF,NQ,MSGLVL)
      INTEGER N,NQ(4),MSGLVL
      DOUBLE PRECISION R(5*N),WRK(N),ALF(1)
C
C.... This routine performs the first step of the Lanczos algorithm.
C.... It performs a step of extended local re-orthogonalization.
C
C.... N      dimension of the eigenproblem
C.... R      an array containing [r(j),q(j),q(j-1),p(j),p(j-1)/Mr(j)]
C.... ALF    diagonal elements of T
C.... NQ(4)  location pointers for the array R
C
C.... BLAS routines:    DATX,DAXPY,DDOT,DSCAL
C.... subroutines:      STARTV
C.... user-supplied:    OP,OPM,STORE
C
      DOUBLE PRECISION RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
C
      DOUBLE PRECISION T,STARTV,DDOT
C
      DOUBLE PRECISION ONE,ZERO
      DATA ONE,ZERO/1.0D0,0.0D0/
C
C.... get initial vector, default is random
C
      RNM = STARTV(N,1,R,WRK,NQ,EPS,MSGLVL)
      IF (RNM.EQ.ZERO) RETURN
C
C.... normalize starting vector
C
      T = ONE/RNM
      CALL DATX(N,T,R,1,R(NQ(1)),1)
      CALL DSCAL(N,T,R(NQ(3)),1)
C
C.... take the first step
C
      CALL LANCZOS_OP(N,R(NQ(3)),R(NQ(1)),R)
      ALF(1) = DDOT(N,R,1,R(NQ(3)),1)
      CALL DAXPY(N,-ALF(1),R(NQ(1)),1,R,1)
C
C.... restore local orthogonality
C
      T = DDOT(N,R,1,R(NQ(3)),1)
      CALL DAXPY(N,-T,R(NQ(1)),1,R,1)
      ALF(1) = ALF(1)+T
      CALL LANCZOS_OPM(N,R,R(NQ(4)))
      RNM = SQRT(DDOT(N,R,1,R(NQ(4)),1))
      ANORM = RNM+ABS(ALF(1))
      TOL = EPSN*ANORM
C
      RETURN
      END
