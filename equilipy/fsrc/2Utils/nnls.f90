!
!
SUBROUTINE nnls (a, m, n, b, x, rnorm, w, indx, mode)
  !     SUBROUTINE nnls(a, m, n, b, x, rnorm, w, indx, mode)
  !
  !  Algorithm NNLS: NONNEGATIVE LEAST SQUARES
  !
  !  The original version of this code was developed by
  !  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
  !  1973 JUN 15, and published in the book
  !  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
  !  Revised FEB 1995 to accompany reprinting of the book by SIAM.
  !
  !  This translation into Fortran 90 by Alan Miller, February 1997
  !  Latest revision - 15 April 1997
!
  !  N.B. The following call arguments have been removed:
  !       mda, zz
  !
  !  GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
  !  N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM
  !
  !                   A * X = B  SUBJECT TO X  >=  0
  !  ------------------------------------------------------------------
  !                  Subroutine Arguments
  !
  !  A(), M, N   ON ENTRY, A() CONTAINS THE M BY N MATRIX, A.
  !              ON EXIT, A() CONTAINS THE PRODUCT MATRIX, Q*A , WHERE Q IS AN
  !              M x M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY THIS SUBROUTINE.
  !  B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CONTAINS Q*B.
  !  X()     ON ENTRY X() NEED NOT BE INITIALIZED.
  !          ON EXIT X() WILL CONTAIN THE SOLUTION VECTOR.
  !  RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE RESIDUAL VECTOR.
  !  W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN THE DUAL
  !          SOLUTION VECTOR.   W WILL SATISFY W(I) = 0. FOR ALL I IN SET P
  !          AND W(I) <= 0. FOR ALL I IN SET Z
  !  INDX()  AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
  !          ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS P AND Z
  !          AS FOLLOWS..
  !              INDX(1)   THRU INDX(NSETP) = SET P.
  !              INDX(IZ1) THRU INDX(IZ2)   = SET Z.
  !              IZ1 = NSETP + 1 = NPP1
  !              IZ2 = N
  !  MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS.
  !          1   THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
  !          2   THE DIMENSIONS OF THE PROBLEM ARE BAD.
  !              EITHER M <= 0 OR N <= 0.
  !          3   ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS.
  !
  !  ------------------------------------------------------------------
  !     ------------------------------------------------------------------
  USE precision
  IMPLICIT NONE
  INTEGER, INTENT(IN)       :: m, n
  INTEGER, INTENT(OUT)      :: mode
  INTEGER, INTENT(OUT), dimension(n)     :: indx
  REAL(8), INTENT(IN OUT), dimension(m,n) :: a
  REAL(8), INTENT(IN OUT), dimension(m) :: b
  REAL(8), INTENT(OUT), dimension(n)    :: x
  REAL(8), INTENT(OUT)    :: rnorm
  REAL(8), INTENT(OUT), dimension(n)    :: w
!
  INTERFACE
    SUBROUTINE g1(a, b, cterm, sterm, sig)
      USE precision
      IMPLICIT NONE
      REAL (dp), INTENT(IN)  :: a, b
      REAL (dp), INTENT(OUT) :: cterm, sterm, sig
    END SUBROUTINE g1
!
    SUBROUTINE h12(mode, lpivot, l1, m, u, up, c, ice, icv, ncv)
      USE precision
      IMPLICIT NONE
      INTEGER, INTENT(IN)                     :: mode, lpivot, l1, m, ice, icv, &
                                                ncv
      REAL (dp), DIMENSION(:), INTENT(IN OUT) :: u, c
      REAL (dp), INTENT(IN OUT)               :: up
    END SUBROUTINE h12
  END INTERFACE
!
  ! Local variables
!
  INTEGER                 :: i, ii, ip, iter, itmax, iz, iz1, iz2, izmax,   &
                            j, jj, jz, l, mda, npp1, nsetp
  REAL (dp), DIMENSION(m) :: zz
  REAL (dp), DIMENSION(1) :: dummy
  REAL (dp)               :: alpha, asave, cc, factor = 0.01_dp, sm, &
                            ss, t, temp, two = 2.0_dp, unorm, up, wmax,    &
                            zero = 0.0_dp, ztest
  !     ------------------------------------------------------------------
  mode = 1
  IF (m <= 0 .OR. n <= 0) THEN
    mode = 2
    RETURN
  END IF
  iter = 0
  itmax = 3*n
!
  !                    INITIALIZE THE ARRAYS indx() AND X().
!
  DO i = 1,n
    x(i) = zero
    indx(i) = i
  END DO
!
  iz2 = n
  iz1 = 1
  nsetp = 0
  npp1 = 1
  !                             ******  MAIN LOOP BEGINS HERE  ******
  !                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
  !                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.
!
  30 IF (iz1 > iz2 .OR. nsetp >= m) GO TO 350
!
  !         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
!
  DO iz = iz1,iz2
    j = indx(iz)
    w(j) = DOT_PRODUCT(a(npp1:m,j), b(npp1:m))
  END DO
!
  !                                   FIND LARGEST POSITIVE W(J).
  60 wmax = zero
  DO iz = iz1,iz2
    j = indx(iz)
    IF (w(j) > wmax) THEN
      wmax = w(j)
      izmax = iz
    END IF
  END DO
!
  !             IF WMAX  <=  0. GO TO TERMINATION.
  !             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
!
  IF (wmax <= zero) GO TO 350
  iz = izmax
  j = indx(iz)
!
  !     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
  !     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
  !     NEAR LINEAR DEPENDENCE.
!
  asave = a(npp1,j)
  CALL h12 (1, npp1, npp1+1, m, a(:,j), up, dummy, 1, 1, 0)
  unorm = zero
  IF (nsetp  /=  0) THEN
    unorm = SUM( a(1:nsetp,j)**2 )
  END IF
  unorm = SQRT(unorm)
  IF (unorm + ABS(a(npp1,j))*factor - unorm  >  zero) THEN
!
  !        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
  !        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).
!
    zz(1:m) = b(1:m)
    CALL h12 (2, npp1, npp1+1, m, a(:,j), up, zz, 1, 1, 1)
    ztest = zz(npp1)/a(npp1,j)
!
  !                                     SEE IF ZTEST IS POSITIVE
!
    IF (ztest > zero) GO TO 140
  END IF
!
  !     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
  !     RESTORE A(NPP1,J), SET W(J) = 0., AND LOOP BACK TO TEST DUAL
  !     COEFFS AGAIN.
!
  a(npp1,j) = asave
  w(j) = zero
  GO TO 60
!
  !     THE INDEX  J = indx(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
  !     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
  !     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
  !     COL J,  SET W(J) = 0.
!
  140 b(1:m) = zz(1:m)
!
  indx(iz) = indx(iz1)
  indx(iz1) = j
  iz1 = iz1+1
  nsetp = npp1
  npp1 = npp1+1
!
  mda = SIZE(a,1)
  IF (iz1  <=  iz2) THEN
    DO jz = iz1,iz2
      jj = indx(jz)
      CALL h12 (2, nsetp, npp1, m, a(:,j), up, a(:,jj), 1, mda, 1)
    END DO
  END IF
!
  IF (nsetp /= m) THEN
    a(npp1:m,j) = zero
  END IF
!
  w(j) = zero
  !                                SOLVE THE TRIANGULAR SYSTEM.
  !                                STORE THE SOLUTION TEMPORARILY IN ZZ().
  CALL solve_triangular(zz)
!
  !                       ******  SECONDARY LOOP BEGINS HERE ******
!
  !                          ITERATION COUNTER.
!
  210 iter = iter+1
  IF (iter > itmax) THEN
    mode = 3
    WRITE (*,'(/a)') ' NNLS quitting on iteration count.'
    GO TO 350
  END IF
!
  !                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.
  !                                  IF NOT COMPUTE ALPHA.
!
  alpha = two
  DO ip = 1,nsetp
    l = indx(ip)
    IF (zz(ip)  <=  zero) THEN
      t = -x(l)/(zz(ip)-x(l))
      IF (alpha > t) THEN
        alpha = t
        jj = ip
      END IF
    END IF
  END DO
!
  !          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL
  !          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.
!
  IF (alpha == two) GO TO 330
!
  !          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO
  !          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.
!
  DO ip = 1,nsetp
    l = indx(ip)
    x(l) = x(l) + alpha*(zz(ip)-x(l))
  END DO
!
  !        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I
  !        FROM SET P TO SET Z.
!
  i = indx(jj)
  260 x(i) = zero
!
  IF (jj /= nsetp) THEN
    jj = jj+1
    DO j = jj,nsetp
      ii = indx(j)
      indx(j-1) = ii
      CALL g1 (a(j-1,ii), a(j,ii), cc, ss, a(j-1,ii))
      a(j,ii) = zero
      DO l = 1,n
        IF (l /= ii) THEN
          !Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))
          temp = a(j-1,l)
          a(j-1,l) = cc*temp + ss*a(j,l)
          a(j,l)   = -ss*temp + cc*a(j,l)
        END IF
      END DO
!
  !                 Apply procedure G2 (CC,SS,B(J-1),B(J))
!
      temp = b(j-1)
      b(j-1) = cc*temp + ss*b(j)
      b(j)   = -ss*temp + cc*b(j)
    END DO
  END IF
!
  npp1 = nsetp
  nsetp = nsetp-1
  iz1 = iz1-1
  indx(iz1) = i
!
  !        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
  !        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
  !        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
  !        THAT ARE NONPOSITIVE WILL BE SET TO ZERO
  !        AND MOVED FROM SET P TO SET Z.
!
  DO jj = 1,nsetp
    i = indx(jj)
    IF (x(i) <= zero) GO TO 260
  END DO
!
  !         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
!
  zz(1:m) = b(1:m)
  CALL solve_triangular(zz)
  GO TO 210
  !                      ******  END OF SECONDARY LOOP  ******
!
  330 DO ip = 1,nsetp
    i = indx(ip)
    x(i) = zz(ip)
  END DO
  !        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.
  GO TO 30
!
  !                        ******  END OF MAIN LOOP  ******
!
  !                        COME TO HERE FOR TERMINATION.
  !                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.
!
  350 sm = zero
  IF (npp1 <= m) THEN
    sm = SUM( b(npp1:m)**2 )
  ELSE
    w(1:n) = zero
  END IF
  rnorm = SQRT(sm)
  RETURN
!
  CONTAINS
!
  SUBROUTINE solve_triangular(zz)
!
  !     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE
  !     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().
!
  REAL (dp), INTENT(IN OUT) :: zz(:)
!
  DO l = 1,nsetp
    ip = nsetp+1-l
    IF (l  /=  1) zz(1:ip) = zz(1:ip) - a(1:ip,jj)*zz(ip+1)
    jj = indx(ip)
    zz(ip) = zz(ip) / a(ip,jj)
  END DO
!
  RETURN
  END SUBROUTINE solve_triangular
!
END SUBROUTINE nnls
!
!
!
!
!
!
!
