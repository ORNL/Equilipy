!
subroutine SubMinCheckDuplicate(lDuplicate)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check whether this yields a duplicate
    ! set of mole fractions as the "other solution phase" that corresponds to
    ! the miscibility gap.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleSubMin
!
    implicit none
!
    integer :: i, j, k, iFirstSUBOther, iLastSUBOther
    real(8) :: dTemp
    logical :: lDuplicate
!
!
    ! Initialize variables:
    iFirstSUBOther = nSpeciesPhase(iSolnPhaseIndexOther-1) + 1
    iLastSUBOther  = nSpeciesPhase(iSolnPhaseIndexOther)
    dTemp       = 0D0
    lDuplicate  = .FALSE.
!
    ! Double check to make sure that they have the same number of constituents:
    if (iLastSUBOther - iFirstSUBOther + 1 /= nVar) then
        INFOThermo = 29
        return
    end if
!
    ! Compute the Euclidean norm between the two mole fraction vectors:
    do k = 1, nVar
!
        ! Absolute species index:
        i = iFirstSUB      + k - 1    ! Absolute index for first species in iSolnPhaseIndex
        j = iFirstSUBOther + k - 1    ! Absolute index for first species in iSolnPhaseIndexOther
!
        dTemp = dTemp + DABS(dMolfraction(i) - dMolfraction(j))
!
    end do
!
    ! Compute the normalized Euclidean norm between the mole fraction vectors between the two "phases":
    dTemp = ( dTemp ) / DFLOAT(nVar)
!
!
    ! Check if the normalized Euclidean norm is less than a specified tolerance:
    if (dTemp < dTolEuclideanNorm) lDuplicate = .TRUE.
!
end subroutine SubMinCheckDuplicate
!
!
