C
C @(#)lanso.f	3.24 (BNP) 6/3/89; from lanso.f 2.18 7/7/88
C
      SUBROUTINE LANSO(N,LANMAX,MAXPRS,ENDL,ENDR,J,NEIG,RITZ,BND,
C     SUBROUTINE LANSO(...,R,WRK,ALF,BET,ETA,OLDETA,NQ,IERR,MSGLVL)
     *   R,WRK,ALF,BET,ETA,OLDETA,NQ,IERR,MSGLVL)
      INTEGER N,LANMAX,MAXPRS,J,NEIG,NQ(4),IERR,MSGLVL
      DOUBLE PRECISION ENDL,ENDR,R(5*N),WRK(N),
     *  ALF(LANMAX),BET(LANMAX+1),RITZ(LANMAX),BND(LANMAX),
     *  ETA(LANMAX),OLDETA(LANMAX)
C
C.... inputs
C.... N      dimension of the eigenproblem
C.... LANMAX upper limit to the number of Lanczos steps
C.... MAXPRS upper limit to the number of wanted eigenpairs
C.... ENDL   left end of the interval containing the wanted eigenvalues
C.... ENDR   right end of the interval containing the wanted eigenvalues
C
C.... work space
C.... R      holds 5 vectors of length n. see the text for details.
C.... NQ(4)  contains the pointers to the begining of each vector in R.
C.... ALF    array to hold diagonal of the tridiagonal T
C.... BET    array to hold off-diagonal of T
C.... ETA    orthogonality estimate of lanczos vectors at step j
C.... OLDETA orthogonality estimate of lanczos vectors at step j-1
C
C.... outputs
C.... J      actual number of Lanczos steps taken
C.... NEIG   number of computed eigenpairs
C.... RITZ   array to hold the converged Ritz values
C.... BND    array to hold the error bounds
C.... IERR   error flag
C
C.... BLAS routines:    DATX,DAXPY,DCOPY,DDOT,DSCAL,IDAMAX
C.... subroutines:      DSORT2,TQLB,ORTBND,PURGE,STARTV,STPONE
C.... user-supplied:    OP,OPM,STORE
C
      INTEGER STORQ,RETRQ,STORP,RETRP
      PARAMETER (STORQ = 1,RETRQ = 2,STORP = 3,RETRP = 4)
C
      INTEGER MAXLL
      DOUBLE PRECISION FOUR
      PARAMETER (MAXLL = 2,FOUR = 4.0D0)
C
      DOUBLE PRECISION RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
      COMMON/RDATA/RNM,ANORM,TOL,EPS,EPS1,RCEPS1,EPSN,REPS,EPS34
C
      LOGICAL ENOUGH
      INTEGER LL,I,L,FIRST,LAST,MID,ID1,ID2,ID3,IDAMAX
      DOUBLE PRECISION T,GAPL,GAP,DDOT,STARTV
C
      DOUBLE PRECISION ONE,ZERO
      DATA ONE,ZERO/1.0D0,0.0D0/
C
      CALL STPONE(N,R,WRK,ALF,NQ,MSGLVL)
      J = 1
      IF (N.EQ.1) THEN
         NEIG = 1
         RITZ(1) = ALF(1)
         BND(1) = ZERO
      ENDIF
C
      ETA(1) = EPS1
      OLDETA(1) = EPS1
C
      LL = 0
      FIRST = 2
      ! jjren
      ! default the lanso use a different LOOP mode
      ! in the dynamic calculation mode we just simplify it
      !LAST = MIN(MAXPRS+MAX(100,MAXPRS),LANMAX)
      LAST = LANMAX
      ENOUGH = .FALSE.
      DO 100 ID1 = 1,MAXPRS
         IF (ENOUGH) GOTO 200
         IF (RNM.LE.TOL) RNM = ZERO
C
C....    lanczos loop
         !WRITE(*,*) "FIRST, LAST", FIRST, LAST
         DO 10 J = FIRST,LAST
            MID = NQ(2)
            NQ(2) = NQ(1)
            NQ(1) = MID
            MID = NQ(3)
            NQ(3) = NQ(4)
            NQ(4) = MID
            CALL STORE(N,STORQ,J-1,R(NQ(2)))
            IF (J-1.LE.MAXLL) CALL STORE(N,STORP,J-1,R(NQ(4)))
            BET(J) = RNM
C
C...        Restart if invariant subspace is found
C
            IF (BET(J).EQ.ZERO) THEN
               RNM = STARTV(N,J,R,WRK,NQ,EPS,MSGLVL)
               ENOUGH = RNM.EQ.ZERO
               IF (ENOUGH) GOTO 15
            ENDIF
C
C....       take a Lanczos step
C
            T = ONE/RNM
            CALL DATX(N,T,R,1,R(NQ(1)),1)
            CALL DSCAL(N,T,R(NQ(3)),1)
C           
            
            CALL LANCZOS_OP(N,R(NQ(3)),R(NQ(1)),R)
            IF (BET(J).GT.ZERO) CALL DAXPY(N,-RNM,R(NQ(2)),1,R,1)
C
            ALF(J) = DDOT(N,R,1,R(NQ(3)),1)
            CALL DAXPY(N,-ALF(J),R(NQ(1)),1,R,1)
C
C....       orthogonalize against initial Lanczos vectors
C
            IF (J.LE.MAXLL+1.AND.ABS(ALF(J-1)).GT.FOUR*ABS(ALF(J))) THEN
               LL = J-1
            ENDIF
            DO 5 I = 1,MIN(LL,J-2)
               CALL STORE(N,RETRP,I,WRK)
               T = DDOT(N,WRK,1,R,1)
               CALL STORE(N,RETRQ,I,WRK)
               CALL DAXPY(N,-T,WRK,1,R,1)
               ETA(I) = EPS1
               OLDETA(I) = EPS1
5           CONTINUE
C
C....       extended local reorthogonalization
C
            T = DDOT(N,R,1,R(NQ(4)),1)
            CALL DAXPY(N,-T,R(NQ(2)),1,R,1)
            IF (BET(J).GT.ZERO) BET(J) = BET(J)+T
            T = DDOT(N,R,1,R(NQ(3)),1)
            CALL DAXPY(N,-T,R(NQ(1)),1,R,1)
            ALF(J) = ALF(J)+T
            CALL LANCZOS_OPM(N,R,R(NQ(4)))
            RNM = SQRT(DDOT(N,R,1,R(NQ(4)),1))
C$$$ Rasmus munk Larsen <rmunk@daimi.aau.dk> observed that ANORM is
C$$$ computed based on local estimate rather than global estimate
C$$$ The following is the suggested change -- Aug. 1998
            ANORM = MAX(BET(J)+ABS(ALF(J))+RNM, ANORM)
C$$$ Previous version
C$$$            ANORM = BET(J)+ABS(ALF(J))+RNM
            TOL = EPSN*ANORM
C
C....       update the orthogonality bounds
C
            CALL ORTBND(J,ALF,BET,ETA,OLDETA)
C
C....       restore the orthogonality state when needed
C
            CALL PURGE(N,LL,J,R,R(NQ(1)),R(NQ(4)),R(NQ(3)),WRK,
C           CALL PURGE(...,ETA,OLDETA,MSGLVL)
     *         ETA,OLDETA,MSGLVL)
            IF (RNM.LE.TOL) RNM = ZERO
10       CONTINUE
         J = LAST
15       IF (ENOUGH) J = J-1
         FIRST = J+1
         BET(J+1) = RNM
C
C....    Now analyze T
C
         L = 1
         DO 40 ID2 = 1,J
            IF (L.GT.J) GOTO 50
            DO 20 I = L,J
               IF (BET(I+1).EQ.ZERO) GOTO 30
20          CONTINUE
            I = J
C....       Now i is at the end of an unreduced submatrix
C
30          CALL DCOPY(I-L+1,ALF(L),1,RITZ(L),-1)
            IF (I.GT.L) CALL DCOPY(I-L,BET(L+1),1,WRK(L+1),-1)
            CALL TQLB(I-L+1,RITZ(L),WRK(L),BND(L),IERR)
            IF (IERR.NE.0) THEN
               PRINT *,' TQLB failed to converge (ierr =',IERR,')'
               PRINT *,' L =',L,' I =',I
               PRINT *,(ID3,RITZ(ID3),WRK(ID3),BND(ID3),ID3 = L,I)
            ENDIF
            DO 35 ID3 = L,I
               BND(ID3) = RNM*ABS(BND(ID3))
35          CONTINUE
            L = I+1
40       CONTINUE
C
C....    Sort eigenvalues into increasing order
C
50       CALL DSORT2(J,RITZ,BND)
C
C....    Massage error bounds for very close Ritz values
C
         MID = IDAMAX(J,BND,1)
         DO 70 L = -1,1,2
            DO 60 I = ((J+1)-L*(J-1))/2,MID-L,L
               IF (ABS(RITZ(I+L)-RITZ(I)).LT.EPS34*ABS(RITZ(I))) THEN
                  IF (BND(I).GT.TOL.AND.BND(I+L).GT.TOL) THEN
                     BND(I+L) = SQRT(BND(I)**2+BND(I+L)**2)
                     BND(I) = ZERO
                  ENDIF
               ENDIF
60          CONTINUE
70       CONTINUE
C
C....    Refine the error bounds
C
         NEIG = 0
         GAPL = RITZ(J)-RITZ(1)
         DO 80 I = 1,J
            GAP = GAPL
            IF (I.LT.J) GAPL = RITZ(I+1)-RITZ(I)
            GAP = MIN(GAP,GAPL)
            IF (GAP.GT.BND(I)) THEN
               BND(I) = BND(I)*(BND(I)/GAP)
            ENDIF
            IF (BND(I).LE.16.0D0*EPS*ABS(RITZ(I))) THEN
               NEIG = NEIG+1
               ENOUGH = ENOUGH.OR.ENDL.LT.RITZ(I).AND.RITZ(I).LT.ENDR
            ENDIF
80       CONTINUE
C
C....    Should we stop?
C
         IF (NEIG.LT.MAXPRS) THEN
            !WRITE(*,*) "NEIG", NEIG, MAXPRS, J
            IF (NEIG.EQ.0) THEN
               LAST = FIRST+100
            ELSE
               LAST = FIRST+MAX(2,((J-6)*(MAXPRS-NEIG))/NEIG)
            ENDIF
            LAST = MIN(LAST,LANMAX)
         ELSE
            ENOUGH = .TRUE.
         ENDIF
         ENOUGH = ENOUGH.OR.FIRST.GT.LANMAX
100   CONTINUE
200   CALL STORE(N,STORQ,J,R(NQ(1)))
      RETURN
      END
