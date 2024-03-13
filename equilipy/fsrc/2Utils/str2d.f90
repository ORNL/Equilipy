!
SUBROUTINE str2d(varname,str,n,m,l)
  use ModuleParseCS
  IMPLICIT NONE
  INTEGER, INTENT(IN)::n,m,l
  character(30), intent(IN)::varname
  CHARACTER, DIMENSION(n,m,l)::str
  INTEGER:: i, j,k
  IF(varname=='cPairNameCS') THEN
    IF(ALLOCATED(cPairNameCS)) DEALLOCATE(cPairNameCS)
    ALLOCATE(cPairNameCS(n,m))
    cPairNameCS=' '
    DO i = 1, n
      DO j = 1, m
        DO k = 1, l
          cPairNameCS(i,j)=trim(cPairNameCS(i,j))//str(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    print*,'Variable name is not supported'
    return
  END IF
END SUBROUTINE str2d
!
!
