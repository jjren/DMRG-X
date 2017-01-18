C
C @(#)ritvec.f	1.7 (BNP) 5/18/89
C
      SUBROUTINE RITVEC(N,J,EV,KAPPA,RITZ,BND,ALF,BET,S,WRK1,WRK2,
     *   IERR,MSGLVL,EIGVEC,MAXPRS)
      INTEGER N,J,EV,IERR,MSGLVL,MAXPRS
      ! MAXPRS is the number of targeted states or highest level
      ! states
      DOUBLE PRECISION KAPPA,RITZ(J),BND(J),ALF(J),BET(J),
     *   S(J,J),WRK1(N),WRK2(N),EIGVEC(N*MAXPRS)
C
C.... subroutines:      TQL2,STORE
C.... BLAS routines:    DAXPY,DCOPY,DSCAL
C
      INTEGER RETRQ
      PARAMETER (RETRQ = 2)
C
      INTEGER I,K
      LOGICAL alive
C
      CALL DSCAL(J*J,0.0D0,S,1)
      DO 10 I = 1,J
         S(I,I) = 1.0D0
10    CONTINUE
      CALL DCOPY(J,ALF,1,WRK1,-1)
      IF (J.GT.1) CALL DCOPY(J-1,BET(2),1,WRK2(2),-1)
      !CALL TQL2(J,J,WRK1,WRK2,S,IERR)
      CALL NEWTQL2(J,WRK1,WRK2(2),S,IERR)
      ! jjren the newtql2 is more robust than the original version
      
      ! jjren change output the diagonal and off-diagonal lanczos element
      ! output the first element of each eigenvector and energy
      inquire(file="lanczos.out",exist=alive)
      if(alive) then
         open(unit=1001, file="lanczos.out", status="old", 
     *       position="append")
      else
         open(unit=1001,file="lanczos.out",status="replace")
      end if

      WRITE(1001,*) "========================="
      WRITE(1001,*) "eigenvalue" 
      WRITE(1001,*) WRK1(1:J)
      WRITE(1001,*) "eigenvector"
      WRITE(1001,*) S(J,1:J)
      WRITE(1001,*) "diagonal"
      WRITE(1001,*) ALF(1:J)
      WRITE(1001,*) "off-diagonal"
      WRITE(1001,*) BET(2:J)

      CLOSE(1001)

      
      IF (IERR.NE.0) RETURN
C
C.... on return WRK1 contains eigenvalues in ascending order
C....       and S contains the corresponding eigenvectors
C
      OPEN(EV,FORM='UNFORMATTED')
      REWIND(EV)
      WRITE(EV)N,J,KAPPA
      DO 50 K = 1,J
         IF ((BND(K).LE.KAPPA*ABS(RITZ(K))) .or. K <= MAXPRS) THEN
C.... Set WRK1 to zeros
            CALL DSCAL(N,0.0D0,WRK1,1)
            DO 20 I = 1,J
               CALL STORE(N,RETRQ,I,WRK2)
               CALL DAXPY(N,S(J-I+1,K),WRK2,1,WRK1,1)
20          CONTINUE
            WRITE(EV) RITZ(K),BND(K),(WRK1(I),I=1,N)
            if(K <= MAXPRS) EIGVEC((K-1)*N+1:K*N) = WRK1(1:N)
         ENDIF
50    CONTINUE
      CLOSE(EV)
      RETURN
      END
