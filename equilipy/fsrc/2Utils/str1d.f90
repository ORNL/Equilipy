!
SUBROUTINE str1d(varname,str,n,l)
  use ModuleParseCS
  IMPLICIT NONE
  INTEGER, INTENT(IN)::n,l
  character(30), intent(IN)::varname
  CHARACTER, DIMENSION(n,l)::str
  CHARACTER(l), DIMENSION(n):: strTemp
  INTEGER:: i, j
  strTemp=' '
  DO i = 1, n
    DO j = 1, l
      strTemp(i)=trim(strTemp(i))//str(i,j)
    END DO
  END DO
!
  IF(varname=='cElementNameCS') THEN
    IF(ALLOCATED(cElementNameCS)) DEALLOCATE(cElementNameCS)
    ALLOCATE(cElementNameCS(n))
    cElementNameCS=strTemp
!
  ELSE IF(varname=='cSolnPhaseTypeCS') THEN
    IF(ALLOCATED(cSolnPhaseTypeCS)) DEALLOCATE(cSolnPhaseTypeCS)
    ALLOCATE(cSolnPhaseTypeCS(n))
    cSolnPhaseTypeCS=strTemp
!
  ELSE IF(varname=='cSolnPhaseNameCS') THEN
    IF(ALLOCATED(cSolnPhaseNameCS)) DEALLOCATE(cSolnPhaseNameCS)
    cSolnPhaseNameCS=strTemp
!
  ELSE IF(varname=='cSpeciesNameCS') THEN
    IF(ALLOCATED(cSpeciesNameCS)) DEALLOCATE(cSpeciesNameCS)
    ALLOCATE(cSpeciesNameCS(n))
    cSpeciesNameCS=strTemp
!
  ELSE IF(varname=='cRegularParamCS') THEN
    IF(ALLOCATED(cRegularParamCS)) DEALLOCATE(cRegularParamCS)
    ALLOCATE(cRegularParamCS(n))
    cRegularParamCS=strTemp
  ELSE
    print*,'Variable name is not supported'
    return
  END IF
END SUBROUTINE str1d
!
!
