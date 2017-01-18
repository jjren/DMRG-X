C
C @(#)store.f	3.3 (BNP) 3/16/89; from store.f 2.2 10/13/87
C
      SUBROUTINE STORE(N,ISW,J,S)
      INTEGER N,ISW,J
      DOUBLE PRECISION S(N)
C
      INTEGER MAXLL
      PARAMETER (MAXLL = 2)
C
      INTEGER STORQ,RETRQ,STORP,RETRP
      PARAMETER (STORQ = 1,RETRQ = 2,STORP = 3,RETRP = 4)
C
      DOUBLE PRECISION A
      ! COMMON/GETPUT/A(125000)
      ! jjren
      ! This memory is set to store the Lanczos vector
      ! can be stored in the disk if the memory is limited 
      ! if lanczos vector dimension is very large
      ! here we only increase the memory space
      COMMON/GETPUT/A(100000000)
C
      IF (ISW.EQ.STORQ) THEN
         CALL DCOPY(N,S,1,A((J+MAXLL-1)*N+1),1)
      ELSE IF (ISW.EQ.RETRQ) THEN
         CALL DCOPY(N,A((J+MAXLL-1)*N+1),1,S,1)
      ELSE IF (ISW.EQ.STORP) THEN
         IF (J.GT.MAXLL) STOP 'STORE: (STORP) J.GT.MAXLL'
         CALL DCOPY(N,S,1,A((J-1)*N+1),1)
      ELSE IF (ISW.EQ.RETRP) THEN
         IF (J.GT.MAXLL) STOP 'STORE: (RETRP) J.GT.MAXLL'
         CALL DCOPY(N,A((J-1)*N+1),1,S,1)
      ENDIF
      RETURN
      END
