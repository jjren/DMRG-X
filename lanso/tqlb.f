C
C @(#)tqlb.f	1.1 (BNP) 5/9/89
C
      SUBROUTINE TQLB(N,D,E,BND,IERR) 
C 
      INTEGER I,J,L,M,N,II,L1,L2,MML,IERR 
      DOUBLE PRECISION D(N),E(N),BND(N)
      DOUBLE PRECISION C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2,TST1,TST2,PYTHAG 
C 
C     THIS SUBROUTINE IS A MODIFICATION OF THE ALGOL PROCEDURE TQL1, 
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND 
C     WILKINSON. 
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971). 
C 
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC 
C     TRIDIAGONAL MATRIX BY THE QL METHOD. 
C 
C     ON INPUT 
C 
C        N IS THE ORDER OF THE MATRIX. 
C 
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX. 
C 
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY. 
C 
C      ON OUTPUT 
C 
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN 
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND 
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE 
C          THE SMALLEST EIGENVALUES. 
C 
C        E HAS BEEN DESTROYED. 
C 
C        BND WILL HOLD THE TOP ELEMENTS OF THE NORMALIZED EIGENVECTORS.
C 
C        IERR IS SET TO 
C          ZERO       FOR NORMAL RETURN, 
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN 
C                     DETERMINED AFTER 30 ITERATIONS. 
C 
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) . 
C 
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW, 
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY 
C 
C     THIS VERSION DATED AUGUST 1983 (MAY 1989).
C 
C     ------------------------------------------------------------------ 
C 
      IERR = 0 
      BND(1) = 1.0D0
      IF (N .EQ. 1) GO TO 1001 
C 
      DO 100 I = 2, N 
         BND(I) = 0.0D0
  100 E(I-1) = E(I) 
C 
      F = 0.0D0 
      TST1 = 0.0D0 
      E(N) = 0.0D0 
C 
      DO 290 L = 1, N 
         J = 0 
         H = ABS(D(L)) + ABS(E(L)) 
         IF (TST1 .LT. H) TST1 = H 
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT .......... 
         DO 110 M = L, N 
            TST2 = TST1 + ABS(E(M)) 
            IF (TST2 .EQ. TST1) GO TO 120 
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT 
C                THROUGH THE BOTTOM OF THE LOOP .......... 
  110    CONTINUE 
C 
  120    IF (M .EQ. L) GO TO 210 
  130    IF (J .EQ. 30) GO TO 1000 
         J = J + 1 
C     .......... FORM SHIFT .......... 
         L1 = L + 1 
         L2 = L1 + 1 
         G = D(L) 
         P = (D(L1) - G) / (2.0D0 * E(L)) 
         R = PYTHAG(P,1.0D0) 
         D(L) = E(L) / (P + SIGN(R,P)) 
         D(L1) = E(L) * (P + SIGN(R,P)) 
         DL1 = D(L1) 
         H = G - D(L) 
         IF (L2 .GT. N) GO TO 145 
C 
         DO 140 I = L2, N 
  140    D(I) = D(I) - H 
C 
  145    F = F + H 
C     .......... QL TRANSFORMATION .......... 
         P = D(M) 
         C = 1.0D0 
         C2 = C 
         EL1 = E(L1) 
         S = 0.0D0 
         MML = M - L 
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- .......... 
         DO 200 II = 1, MML 
            C3 = C2 
            C2 = C 
            S2 = S 
            I = M - II 
            G = C * E(I) 
            H = C * P 
            R = PYTHAG(P,E(I)) 
            E(I+1) = S * R 
            S = E(I) / R 
            C = P / R 
            P = C * D(I) - S * G 
            D(I+1) = H + S * (C * G + S * D(I)) 
            H = BND(I+1)
            BND(I+1) = S*BND(I)+C*H
            BND(I) = C*BND(I)-S*H
  200    CONTINUE 
C 
         P = -S * S2 * C3 * EL1 * E(L) / DL1 
         E(L) = S * P 
         D(L) = C * P 
         TST2 = TST1 + ABS(E(L)) 
         IF (TST2 .GT. TST1) GO TO 130 
  210    P = D(L) + F 
         H = BND(L)
C     .......... ORDER EIGENVALUES .......... 
         IF (L .EQ. 1) GO TO 250 
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- .......... 
         DO 230 II = 2, L 
            I = L + 2 - II 
            IF (P .GE. D(I-1)) GO TO 270 
            D(I) = D(I-1) 
            BND(I) = BND(I-1)
  230    CONTINUE 
C 
  250    I = 1 
  270    D(I) = P 
         BND(I) = H
  290 CONTINUE 
C 
      GO TO 1001 
C     .......... SET ERROR -- NO CONVERGENCE TO AN 
C                EIGENVALUE AFTER 30 ITERATIONS .......... 
 1000 IERR = L 
 1001 RETURN 
      END 
