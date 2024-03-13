SUBROUTINE g1(a, b, cterm, sterm, sig)
!
  !     COMPUTE ORTHOGONAL ROTATION MATRIX..
!
  !  The original version of this code was developed by
  !  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
  !  1973 JUN 12, and published in the book
  !  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
  !  Revised FEB 1995 to accompany reprinting of the book by SIAM.
!
  !     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
  !                        (-S,C)         (-S,C)(B)   (   0          )
  !     COMPUTE SIG = SQRT(A**2+B**2)
  !        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
  !        SIG MAY BE IN THE SAME LOCATION AS A OR B .
  !     ------------------------------------------------------------------
  USE precision
  IMPLICIT NONE
  REAL (dp), INTENT(IN)  :: a, b
  REAL (dp), INTENT(OUT) :: cterm, sterm, sig
!
  !     Local variables
  REAL (dp) :: one = 1.0D0, xr, yr, zero = 0.0D0
  !     ------------------------------------------------------------------
  IF (ABS(a) > ABS(b)) THEN
    xr = b / a
    yr = SQRT(one + xr**2)
    cterm = SIGN(one/yr, a)
    sterm = cterm * xr
    sig = ABS(a) * yr
    RETURN
  END IF
!
  IF (b /= zero) THEN
    xr = a / b
    yr = SQRT(one + xr**2)
    sterm = SIGN(one/yr, b)
    cterm = sterm * xr
    sig = ABS(b) * yr
    RETURN
  END IF
!
  !      SIG = ZERO
  cterm = zero
  sterm = one
  RETURN
END SUBROUTINE g1
!
!
