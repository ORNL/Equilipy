!
SUBROUTINE str3d(varname,str,n,m,o,l)
  use ModuleParseCS
  IMPLICIT NONE
  INTEGER, INTENT(IN)::n,m,o,l
  character(30), intent(IN)::varname
  CHARACTER, DIMENSION(n,m,o,l)::str
  INTEGER:: i,j,k,y
  IF(varname=='cConstituentNameSUBCS') THEN
    IF(ALLOCATED(cConstituentNameSUBCS)) DEALLOCATE(cConstituentNameSUBCS)
    ALLOCATE(cConstituentNameSUBCS(n,m,o))
    cConstituentNameSUBCS=' '
    DO i = 1, n
      DO j = 1, m
        DO k = 1, o
          DO y = 1, l
            cConstituentNameSUBCS(i,j,k)=&
            trim(cConstituentNameSUBCS(i,j,k))//str(i,j,k,y)
          END DO
        END DO
      END DO
    END DO
  ELSE
    print*,'Variable name is not supported'
    return
  END IF
END SUBROUTINE str3d
!
!
