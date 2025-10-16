!
SUBROUTINE h12(mode, lpivot, l1, m, u, up, c, ice, icv, ncv)
!
  !     SUBROUTINE h12 (mode, lpivot, l1, m, u, up, c, ice, icv, ncv)
!
  !  CONSTRUCTION AND/OR APPLICATION OF A SINGLE
  !  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B
!
  !  The original version of this code was developed by
  !  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
  !  1973 JUN 12, and published in the book
  !  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
  !  Revised FEB 1995 to accompany reprinting of the book by SIAM.
  !     ------------------------------------------------------------------
  !                     Subroutine Arguments
!
  !     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
  !            Householder transformation, or Algorithm H2 to apply a
  !            previously constructed transformation.
  !     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
  !     L1,M   IF L1  <=  M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
  !            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M
  !            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
  !     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
  !            vector.  IUE is the storage increment between elements.
  !            On exit when MODE = 1, U() and UP contain quantities
  !            defining the vector U of the Householder transformation.
  !            on entry with MODE = 2, U() and UP should contain
  !            quantities previously computed with MODE = 1.  These will
  !            not be modified during the entry with MODE = 2.
  !     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
  !            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
  !            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
  !            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
  !     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
  !     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
  !     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV  <=  0
  !            NO OPERATIONS WILL BE DONE ON C().
  !     ------------------------------------------------------------------
  !     ------------------------------------------------------------------
!
  USE precision
  IMPLICIT NONE
  INTEGER, INTENT(IN)                     :: mode, lpivot, l1, m, ice, icv, ncv
  REAL (dp), DIMENSION(:), INTENT(IN OUT) :: u, c
  REAL (dp), INTENT(IN OUT)               :: up
!
  !  Local variables
  INTEGER          :: i, i2, i3, i4, incr, j
  REAL (dp)        :: b, cl, clinv, one = 1.0D0, sm
  !     ------------------------------------------------------------------
  IF (0 >= lpivot .OR. lpivot >= l1 .OR. l1 > m) RETURN
  cl = ABS(u(lpivot))
  IF (mode /= 2) THEN
  !                            ****** CONSTRUCT THE TRANSFORMATION. ******
    DO j = l1, m
      cl = MAX(ABS(u(j)),cl)
    END DO
    IF (cl <= 0) RETURN
    clinv = one / cl
    sm = (u(lpivot)*clinv) ** 2 + SUM( (u(l1:m)*clinv)**2 )
    cl = cl * SQRT(sm)
    IF (u(lpivot) > 0) THEN
      cl = -cl
    END IF
    up = u(lpivot) - cl
    u(lpivot) = cl
  ELSE
  !            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
!
    IF (cl <= 0) RETURN
  END IF
  IF (ncv <= 0) RETURN
!
  b = up * u(lpivot)
  !                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
!
  IF (b < 0) THEN
    b = one / b
    i2 = 1 - icv + ice * (lpivot-1)
    incr = ice * (l1-lpivot)
    DO j = 1, ncv
      i2 = i2 + icv
      i3 = i2 + incr
      i4 = i3
      sm = c(i2) * up
      DO i = l1, m
        sm = sm + c(i3) * u(i)
        i3 = i3 + ice
      END DO
      IF (sm /= 0) THEN
        sm = sm * b
        c(i2) = c(i2) + sm * up
        DO i = l1, m
          c(i4) = c(i4) + sm * u(i)
          i4 = i4 + ice
        END DO
      END IF
    END DO ! j = 1, ncv
  END IF
!
RETURN
END SUBROUTINE h12
!
!
