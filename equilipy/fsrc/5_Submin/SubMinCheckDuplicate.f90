!
subroutine SubMinCheckDuplicate(lDuplicate)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SubMinCheckDuplicate.f90
    !> \brief   Check whether a miscibility-gap subminimization result duplicates another phase copy.
    !> \author  M.H.A. Piro
    !> \date    Aug. 30, 2012
    !> \sa      Subminimization.f90
    !> \sa      ModuleSubMin.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   08/30/2012      M.H.A. Piro         Original duplicate-composition check
    !   07/20/2026      S.Y. Kwon           Used the shared Euclidean composition tolerance to identify duplicate subminimization solutions.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check whether the current
    !! subminimized composition duplicates the corresponding miscibility-gap
    !! copy of the same solution phase.  Duplicate minima are rejected so
    !! CheckPhaseAssemblage/PEA does not add the same local minimum twice.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! iSolnPhaseIndexOther Matching miscibility-gap solution phase index.
    ! iFirstSUB            First absolute species index for the active solution phase.
    ! iLastSUB             Last absolute species index for the active solution phase.
    ! nVar                 Number of species/endmembers in the active solution phase.
    ! dMolFraction         Current subminimized mole fractions and the other miscibility-copy fractions.
    ! dStoichSpecies       Endmember stoichiometry used to compare phase compositions.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out]  lDuplicate  True when the composition difference is below dTolEuclideanNorm.
    !
    ! INFOThermo            Set to 29 if the paired miscibility-gap phases do not have matching sizes.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! Subminimization       Calls this after a miscibility-gap subminimization trial.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - The active and paired miscibility-gap phases have the same endmember count.
    ! - Duplicate detection compares composition-space stoichiometry differences,
    !   not raw site-fraction differences.
    ! - dTolEuclideanNorm is the normalized composition-difference threshold.
    !
    !-------------------------------------------------------------------------------------------------------------
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

    dTemp = MAXVAL(DABS(MATMUL(dMolfraction(iFirstSUB:iLastSUB),dStoichSpecies(iFirstSUB:iLastSUB,:)) -&
    MATMUL(dMolfraction(iFirstSUBOther:iLastSUBOther),dStoichSpecies(iFirstSUBOther:iLastSUBOther,:))))

!
    ! Compute the normalized Euclidean norm between the mole fraction vectors between the two "phases":
    dTemp = ( dTemp ) / DFLOAT(nVar)
!
!
    ! Check if the normalized composition difference is less than the duplicate tolerance:
    if (dTemp < dTolEuclideanNorm) lDuplicate = .TRUE.
!
end subroutine SubMinCheckDuplicate
!
!
