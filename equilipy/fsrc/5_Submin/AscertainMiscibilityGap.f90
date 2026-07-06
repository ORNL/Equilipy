!> \brief Decide whether two same-name solution candidates represent a distinct miscibility gap.
!!
!! \details Evaluates the solution Gibbs energy at two candidate compositions
!! and at their midpoint.  The candidate compositions must be applied to both
!! dMolFraction and dMolesSpecies because several Gibbs models, including
!! SUBG/SUBQ, derive pair/site variables from dMolFraction.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    AscertainMiscibilityGap.f90
!> \brief   Decide whether two solution candidates are distinct miscibility-gap branches.
!> \author  S.Y. Kwon
!> \date    Jun. 27, 2026
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/27/2026      S.Y. Kwon           Restored phase-local state and evaluated duplicate checks with
!                                       both dMolFraction and dMolesSpecies synchronized.
!   06/27/2026      S.Y. Kwon           Preserved separate candidates when duplicate-check Gibbs evaluations
!                                       are invalid instead of printing diagnostics.
!
!
! Purpose:
! ========
!
!> \details Leveling-to-Lagrangian translation calls this routine before
!! merging same-name solution candidates.  A convex midpoint means the two
!! candidates are not separated by a miscibility gap and can be merged; a
!! concave midpoint means they should remain as separate solution phases.
!
!
! Required input variables:
! =========================
!
!> \param[in]  iSolnPhaseIndex       Solution phase index used for Gibbs evaluation.
!> \param[in]  iFirst                First species row of the candidate phase.
!> \param[in]  iLast                 Last species row of the candidate phase.
!> \param[in]  nConstituents         Number of phase-local constituents.
!> \param[in]  dMolFractionMiscible  New candidate composition.
!> \param[in]  dMolFractionImmiscible Existing candidate composition.
!
! dMolFraction          Phase-local mole fractions used by SUBG/SUBQ and CEF Gibbs models.
! dMolesSpecies         Phase-local species moles used to accumulate the solution Gibbs energy.
!
!
! Output/updated variables:
! =========================
!
!> \param[out] lDuplicate            True when the two candidates should be merged.
!
! The routine restores dMolFraction, dMolesSpecies, dChemicalPotential,
! phase-local property arrays, and dGibbsSolnPhase before returning.
!
!
! Called subroutines/functions:
! =============================
!
! CompExcessGibbsEnergy              Evaluates the phase Gibbs energy at each candidate composition.
!
!
! Primary callers:
! ================
!
! Subminimization                    Rejects duplicate miscibility minima during phase-local searches.
! Level2Lagrange                     Prevents accidental merging of true miscibility-gap candidates.
!
!
! Numerical assumptions:
! ======================
!
! - Input candidate vectors are phase-local mole fractions and are expected to
!   be normalized by the caller.
! - A midpoint no higher than the average endpoint Gibbs energy is treated as
!   duplicate/mergeable on the current composition line.
!
!-------------------------------------------------------------------------------------------------------------



subroutine AscertainMiscibilityGap(iSolnPhaseIndex,iFirst,iLast,nConstituents,dMolFractionMiscible,dMolFractionImmiscible,lDuplicate)

    USE ModuleThermo 
    USE ModuleGEMSolver
!
    implicit none
!
    integer, intent(in)::  iSolnPhaseIndex,iFirst, iLast, nConstituents
    real(8), intent(in):: dMolFractionMiscible(nConstituents),dMolFractionImmiscible(nConstituents)
    logical, intent(out)::  lDuplicate
    real(8):: dGMiscible,dGImmiscible,dGMiddle,dGibbsSolnPhaseDummy
    real(8),dimension(nConstituents)    :: dMolesSpeciesDummy, dMolFractionDummy
    real(8),dimension(nConstituents)    :: dChemicalPotentialDummy
    real(8),dimension(nConstituents)    :: dPartialEnthalpyDummy, dPartialEntropyDummy, dPartialHeatCapacityDummy
    real(8),dimension(nConstituents)    :: dPartialExcessGibbsDummy
    real(8),dimension(nConstituents)    :: dPartialEnthalpyXSDummy, dPartialEntropyXSDummy, dPartialHeatCapacityXSDummy
!
    ! Initialize variables:
    lDuplicate        = .FALSE.
    dMolesSpeciesDummy = dMolesSpecies(iFirst:iLast)
    dMolFractionDummy = dMolFraction(iFirst:iLast)
    dChemicalPotentialDummy = dChemicalPotential(iFirst:iLast)
    dPartialEnthalpyDummy = dPartialEnthalpy(iFirst:iLast)
    dPartialEntropyDummy = dPartialEntropy(iFirst:iLast)
    dPartialHeatCapacityDummy = dPartialHeatCapacity(iFirst:iLast)
    dPartialExcessGibbsDummy = dPartialExcessGibbs(iFirst:iLast)
    dPartialEnthalpyXSDummy = dPartialEnthalpyXS(iFirst:iLast)
    dPartialEntropyXSDummy = dPartialEntropyXS(iFirst:iLast)
    dPartialHeatCapacityXSDummy = dPartialHeatCapacityXS(iFirst:iLast)
    dGibbsSolnPhaseDummy = dGibbsSolnPhase(iSolnPhaseIndex)

    dGMiscible = 0D0
    dGImmiscible = 0D0
    dGMiddle = 0D0

    dMolFraction(iFirst:iLast) = dMolFractionMiscible
    dMolesSpecies(iFirst:iLast) = dMolFractionMiscible
    dGibbsSolnPhase(iSolnPhaseIndex) = 0D0
    call CompExcessGibbsEnergy(iSolnPhaseIndex)
    dGMiscible = dGibbsSolnPhase(iSolnPhaseIndex)

    dMolFraction(iFirst:iLast) = dMolFractionImmiscible
    dMolesSpecies(iFirst:iLast) = dMolFractionImmiscible
    dGibbsSolnPhase(iSolnPhaseIndex) = 0D0
    call CompExcessGibbsEnergy(iSolnPhaseIndex)
    dGImmiscible = dGibbsSolnPhase(iSolnPhaseIndex)

    dMolFraction(iFirst:iLast) = 0.5D0*(dMolFractionImmiscible+dMolFractionMiscible)
    dMolesSpecies(iFirst:iLast) = dMolFraction(iFirst:iLast)
    dGibbsSolnPhase(iSolnPhaseIndex) = 0D0
    call CompExcessGibbsEnergy(iSolnPhaseIndex)
    dGMiddle = dGibbsSolnPhase(iSolnPhaseIndex)

    if ((dGMiscible /= dGMiscible).OR.(dGImmiscible /= dGImmiscible).OR.(dGMiddle /= dGMiddle)) then
        lDuplicate = .FALSE.
    else if (dGMiddle-0.5D0*(dGMiscible +dGImmiscible)<= 1D-12) then
        lDuplicate = .TRUE.
    end if
    dMolesSpecies(iFirst:iLast) = dMolesSpeciesDummy
    dMolFraction(iFirst:iLast) = dMolFractionDummy
    dChemicalPotential(iFirst:iLast) = dChemicalPotentialDummy
    dPartialEnthalpy(iFirst:iLast) = dPartialEnthalpyDummy
    dPartialEntropy(iFirst:iLast) = dPartialEntropyDummy
    dPartialHeatCapacity(iFirst:iLast) = dPartialHeatCapacityDummy
    dPartialExcessGibbs(iFirst:iLast) = dPartialExcessGibbsDummy
    dPartialEnthalpyXS(iFirst:iLast) = dPartialEnthalpyXSDummy
    dPartialEntropyXS(iFirst:iLast) = dPartialEntropyXSDummy
    dPartialHeatCapacityXS(iFirst:iLast) = dPartialHeatCapacityXSDummy
    dGibbsSolnPhase(iSolnPhaseIndex) = dGibbsSolnPhaseDummy
    ! lDuplicate = .TRUE.

!
!
    return
!
end subroutine AscertainMiscibilityGap
!
!
