!-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergySUBL.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of solution phase constituents in a SUBL
    !!           or SUBLM solution phase.
    !> \author  M.H.A. Piro
    !> \date    January 17, 2013
    !> \sa      CompExcessGibbsEnergy.f90
    !> \sa      CompExcessGibbsEnergyRKMP.f90
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/17/2013      M.H.A. Piro         Original code.
    !   02/11/2013      M.H.A. Piro         Added capability of handling non-ideal terms.
    !   02/14/2013      M.H.A. Piro         Corrected handling of non-ideal terms with higher order
    !                                       mixing parameters (happy valentine's day).
    !   02/15/2013      M.H.A. Piro         Fix bug in handling higher order terms from yesterday.
    !   02/27/2013      M.H.A. Piro         The sum of site fractions may not necessarily equal unity. Apply
    !                                       an appropriate correction.
    !   07/20/2026      S.Y. Kwon           Replaced legacy CEF partial properties with site-fraction derivatives and retained raw ordered SUBOM energies.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the partial molar excess Gibbs energy of mixing
    !! (dPartialExcessGibbs) of all constituents in a non-ideal solution phase designated as 'SUBL'
    !! (Compound Energy Formalism) or 'SUBLM' (SUBL-Magnetic) phases.  For more information on the
    !! SUBL model and the derivation of equations, the reader is referred to the following literature:
    !!
    !!      B. Sundman and J. Agren, "A Regular Solution Model for Phases with Several Components
    !!      and Sublattices, Suitable for Computer Applications," Journal of Physical Chemistry of
    !!      Solids, 42 (1981) 297-301.
    !!
    !!      M. Hillert, "The Compound Energy Formalism," Journal of Alloys and Compounds, 320 (2001)
    !!      161-176.
    !!
    !! An important distinction between a phase represented by the SUBL model and any other model
    !! is that the phase components and the constituents are not the same entity (refer to above
    !! literature for a theoretical discussion).  As a result of this, the contribution to the chemical
    !! potential term from the reference molar Gibbs energy is NOT \f$ g_i^{\circ} \f$, but rather a
    !! more complicated expression involved all components in this phase.
    !!
    !! The chemical potential of a component in a SUBL phase is defined by the following equation:
    !!
    !! 	\f$ \mu_{i(\lambda)} = \sum_{j=1}^{N_{\lambda}} x_{j(\lambda)} g_{j(\lambda)}^{\circ}
    !!      \left(1 - N_s + \sum_{s=1}^{N_s} \frac{\delta_{i,j}}{y_{j(s)}} \right)
    !!      + RT \sum_{s=1}^{N_s} a_s \mathrm{ln} (y_{i(s)}) + g_{i(\lambda)}^{ex} \f$
    !!
    !! The partial molar excess Gibbs energy of mixing of a component of a SUBL phase is:
    !!
    !! 	\f$ g_{i(\lambda)}^{ex} = \sum_{p=1}^{N_p} \left( \prod_{m=1} y_{m(s)}  \right)
    !!      \sum_{z=0} {^zL_{j,k}} \left(  1 - (N_s + z) + \sum_{s=1}^{N_s}
    !!      \frac{\delta _{i,p}}{y_{i(s)}}   \right) \f$
    !!
    !! The mole fraction of any component in a SUBL phase is related to the site fractions of the
    !! constituents through the following multiplicative relationship:
    !!
    !! \f$ x_{i(\lambda)} = \prod_{s=1}^{N_s} y_{i(s)} \f$
    !!
    !! Similarly, the site fraction of any constituent is related to the mole fraction through the following
    !! summation:
    !!
    !! \f$ y_{c(s)} = \sum_{i=1}^{N_{\lambda}} x_{i(\lambda)} \delta_{i,c(s)} \f$
    !!
    !! Since the molar Gibbs energy terms are defined by the phase components and NOT the constituents, the
    !! numerics necessarily works with mole fractions (corresponding to components) rather than site fractions
    !! (corresponding to constituents).  Therefore, the mole fractions are modified at any iteration and the
    !! site fractions are computed from the mole fractions.  Note that the computation of the site fractions
    !! is a summation operation and it is thus insensitive to components with relatively small mole fractions.
    !! However, the calculation of a mole fraction from the site fractions is sensitive because it is a
    !! multiplicative function. A correction is applied to the mole fractions of all components after the site
    !! fractions are computed by defining the mole fraction as a function of the site fractions.
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iSolnIndex    Absolute index of a solution phase
    !
    ! nSpeciesPhase             An integec vector representing the last index of a species in a particular
    !                            solution phase
    ! nParamPhase               An integer vector representing the number of parameters for a particular phase.
    ! nSublattice               An integer scalar representing the number of sublattices associated with this
    !                            particular phase (used for convenience).
    ! nSublatticePhase          An integer vector representing the number of sublattices for each charged phase.
    ! iChargedPhaseID           An integer scalar representing the relative index of the charged phase.
    ! iParam                    An integer scalar representing the index number of a parameter.
    ! iPhaseSublattice          An integer vector representing the relative index of a charged phase corresponding
    !                            to each solution phase in the system.  Thus, it is equal to zero for a phase
    !                            that does not contain any sublattices.
    ! iRegularParam             An integer matrix representing the relative indices of species mixing for a
    !                            particular parmater.
    ! iFirst             An integer scalar representing the absolute index of the first species in a phase.
    ! iLast              An integer scalar representing the absolute index of the last species in a phase.
    ! dChemicalPotential        A double real vector representing the chemical potential for every species in the
    !                            system.
    ! dMolFraction              A double real vector representing hte mole fraction of all species in the system.
    ! dSiteFraction             A double real array representing the site fraction of each constituent on each
    !                            sublattice.  The first dimension corresponds to the charged phase index, the
    !                            second dimension corresponds to the sublattice index and the third dimension
    !                            corresponds to the constituent index.
    ! dExcessGibbsParam         A double real vector representing the molar excess Gibbs energy of mixing for
    !                            each subsystem.
    ! dPartialExcessGibbs       Partial molar excess Gibbs energy of mixing of species.
    ! dMolFraction              Current estimated mole fraction.
    ! cSolnPhaseType            A character vector representing the solution phase type.
    ! iConstituentSublattice    An integer array representing the constituent index for a particular sublattice
    !                            phase.  The first dimension refers to the relative phase index for charged phases,
    !                            the second dimension represents the sublattice index and the third dimension
    !                            corresponds to the relative component index (not absolute).
    ! iFirstParam               An integer scalar representing the constituent index of the first constituent
    !                            that is being mixed for a mixing parameter.
    ! iSeconParam               An integer scalar representing the constituent index of the second constituent
    !                            that is being mixed for a mixing parameter.
    ! iSubParam                 An integer scalar representing the sublattice index corresponding to the mixing
    !                            parameter.
    !
    !-------------------------------------------------------------------------------------------------------------
!
!
subroutine CompExcessGibbsEnergySUBL(iSolnIndex)
!
    USE ModuleThermo
    USE ModuleThermoIO
!
    implicit none
!
    integer, intent(in) :: iSolnIndex

    integer :: iFirst, iLast, iPhaseID, nSiteCapacity, nSiteOut, iInfo
    real(8) :: dScalarGibbs, dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
    real(8), dimension(:), allocatable :: dSiteGradientG, dSiteGradientH
    real(8), dimension(:), allocatable :: dSiteGradientS, dSiteGradientCp
!
!
    if ((cSolnPhaseType(iSolnIndex) /= 'SUBL').AND.(cSolnPhaseType(iSolnIndex) /= 'SUBLM').AND. &
        (cSolnPhaseType(iSolnIndex) /= 'SUBOM')) return

    iPhaseID = iPhaseSublattice(iSolnIndex)
    if (iPhaseID <= 0) then
        INFOThermo = 36
        return
    end if

    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnIndex)
    nSiteCapacity = nMaxSublatticeSys * nMaxConstituentSys

    allocate(dSiteGradientG(nSiteCapacity), dSiteGradientH(nSiteCapacity), &
        dSiteGradientS(nSiteCapacity), dSiteGradientCp(nSiteCapacity))

    call CompGradientSUBL(iSolnIndex, nSiteCapacity, dSiteGradientG, dScalarGibbs, &
        dSiteGradientH, dSiteGradientS, dSiteGradientCp, dScalarEnthalpy, dScalarEntropy, &
        dScalarHeatCapacity, nSiteOut, iInfo)

    if (iInfo /= 0) then
        if (iInfo == 4) print *, 'Unrecognized excess mixing term in SUBL phase ', cSolnPhaseName(iSolnIndex)
        INFOThermo = 36
        if (allocated(dSiteGradientG)) deallocate(dSiteGradientG)
        if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
        if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
        if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)
        return
    end if

    call CorrectSUBLMoleFractionsFromSites(iSolnIndex, iPhaseID, iFirst, iLast)
    call ProjectSUBLSiteGradientsToEndmembers(iSolnIndex, iPhaseID, iFirst, iLast, nSiteOut, &
        dScalarGibbs, dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity, &
        dSiteGradientG, dSiteGradientH, dSiteGradientS, dSiteGradientCp)

    dPartialExcessGibbs(iFirst:iLast) = 0D0
    dPartialEnthalpyXS(iFirst:iLast) = 0D0
    dPartialEntropyXS(iFirst:iLast) = 0D0
    dPartialHeatCapacityXS(iFirst:iLast) = 0D0
    dMagGibbsEnergy(iFirst:iLast) = 0D0
    dMagEnthalpy(iFirst:iLast) = 0D0
    dMagEntropy(iFirst:iLast) = 0D0
    dMagHeatCapacity(iFirst:iLast) = 0D0

    if (allocated(dSiteGradientG)) deallocate(dSiteGradientG)
    if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
    if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
    if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)

    return
!
end subroutine CompExcessGibbsEnergySUBL
!
!
subroutine CorrectSUBLMoleFractionsFromSites(iSolnIndex, iPhaseID, iFirst, iLast)
!
    USE ModuleThermo
!
    implicit none
!
    integer, intent(in) :: iSolnIndex, iPhaseID, iFirst, iLast
    integer :: i, m, s, c
    real(8) :: dProduct
!
    do i = iFirst, iLast
        m = i - iFirst + 1
        dProduct = 1D0
        do s = 1, nSublatticePhase(iPhaseID)
            c = iConstituentSublattice(iPhaseID,s,m)
            dProduct = dProduct * dSiteFraction(iPhaseID,s,c)
        end do
        dMolFraction(i) = dProduct
    end do
!
    return
!
end subroutine CorrectSUBLMoleFractionsFromSites
!
!
subroutine ProjectSUBLSiteGradientsToEndmembers(iSolnIndex, iPhaseID, iFirst, iLast, nSiteOut, &
    dScalarGibbs, dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity, &
    dSiteGradientG, dSiteGradientH, dSiteGradientS, dSiteGradientCp)
!
    USE ModuleThermo
!
    implicit none
!
    integer, intent(in) :: iSolnIndex, iPhaseID, iFirst, iLast, nSiteOut
    real(8), intent(in) :: dScalarGibbs, dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
    real(8), intent(in) :: dSiteGradientG(nSiteOut), dSiteGradientH(nSiteOut)
    real(8), intent(in) :: dSiteGradientS(nSiteOut), dSiteGradientCp(nSiteOut)

    integer :: i, j, m, s, c, a, nSublattice, iCursor
    integer, dimension(nMaxSublatticeSys) :: iSiteOffset
    real(8) :: dStoichSum
    real(8), dimension(nMaxSublatticeSys) :: dAverageG, dAverageH, dAverageS, dAverageCp
!
    nSublattice = nSublatticePhase(iPhaseID)
    iSiteOffset = 0
    iCursor = 0
    do s = 1, nSublattice
        iSiteOffset(s) = iCursor
        iCursor = iCursor + nConstituentSublattice(iPhaseID,s)
    end do
!
    dAverageG = 0D0
    dAverageH = 0D0
    dAverageS = 0D0
    dAverageCp = 0D0
    do s = 1, nSublattice
        do c = 1, nConstituentSublattice(iPhaseID,s)
            a = iSiteOffset(s) + c
            dAverageG(s) = dAverageG(s) + dSiteFraction(iPhaseID,s,c) * dSiteGradientG(a)
            dAverageH(s) = dAverageH(s) + dSiteFraction(iPhaseID,s,c) * dSiteGradientH(a)
            dAverageS(s) = dAverageS(s) + dSiteFraction(iPhaseID,s,c) * dSiteGradientS(a)
            dAverageCp(s) = dAverageCp(s) + dSiteFraction(iPhaseID,s,c) * dSiteGradientCp(a)
        end do
    end do
!
    do i = iFirst, iLast
        m = i - iFirst + 1
        dChemicalPotential(i) = dScalarGibbs
        dPartialEnthalpy(i) = dScalarEnthalpy
        dPartialEntropy(i) = dScalarEntropy
        dPartialHeatCapacity(i) = dScalarHeatCapacity
!
        do s = 1, nSublattice
            c = iConstituentSublattice(iPhaseID,s,m)
            a = iSiteOffset(s) + c
            dChemicalPotential(i) = dChemicalPotential(i) + dSiteGradientG(a) - dAverageG(s)
            dPartialEnthalpy(i) = dPartialEnthalpy(i) + dSiteGradientH(a) - dAverageH(s)
            dPartialEntropy(i) = dPartialEntropy(i) + dSiteGradientS(a) - dAverageS(s)
            dPartialHeatCapacity(i) = dPartialHeatCapacity(i) + dSiteGradientCp(a) - dAverageCp(s)
        end do
!
        dStoichSum = 0D0
        do j = 1, nElements
            dStoichSum = dStoichSum + DABS(dStoichSpecies(i,j))
        end do
        if ((dStoichSum == 0D0).AND.(dMolFraction(i) > 0.9D0)) dChemicalPotential(i) = dChemicalPotential(i) + 1D3
    end do
!
    return
!
end subroutine ProjectSUBLSiteGradientsToEndmembers
!
!
subroutine ApplyFullOrderDisorderCorrectionSUBL(iSolnIndex)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergySUBL.f90
    !> \brief   Apply the full SUBOM order/disorder correction for public CEF evaluations.
    !> \author  S.Y. Kwon
    !> \date    Jun. 24, 2026
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Kept raw ordered SUBOM energy evaluation internal to the sublattice model.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details Public thermodynamic callers should see full SUBOM Gibbs energy:
    !! raw ordered CEF plus the disorder correction.  The raw ordered expression
    !! is used only inside ApplyOrderDisorderCorrectionSUBL when it evaluates
    !! collapsed/random states.  The lOrderDisorderEvaluation guard identifies
    !! that internal recursion and prevents applying the correction twice.
    !
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in] iSolnIndex  Absolute solution phase index.
    !
    ! cSolnPhaseType            Solution model type for iSolnIndex.
    ! lOrderDisorderEvaluation  True only during internal order/disorder raw-state evaluation.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! dPartialExcessGibbs       Receives the SUBOM disorder correction when applicable.
    ! dPartialEnthalpyXS        Receives the SUBOM disorder correction when applicable.
    ! dPartialEntropyXS         Receives the SUBOM disorder correction when applicable.
    ! dPartialHeatCapacityXS    Receives the SUBOM disorder correction when applicable.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! ApplyOrderDisorderCorrectionSUBL  Computes and applies the full SUBOM correction.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! CompExcessGibbsEnergySUBL         Calls this before returning from the public SUBL/SUBOM evaluator.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - This helper is intentionally a no-op for SUBL/SUBLM phases.
    ! - This helper is intentionally a no-op while ApplyOrderDisorderCorrectionSUBL is evaluating raw
    !   ordered or disordered states.  That is the only supported raw SUBOM path.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
!
    implicit none
!
    integer, intent(in) :: iSolnIndex
!
    if ((TRIM(cSolnPhaseType(iSolnIndex)) == 'SUBOM').AND.(.NOT.lOrderDisorderEvaluation)) then
        call ApplyOrderDisorderCorrectionSUBL(iSolnIndex)
    end if
!
    return
!
end subroutine ApplyFullOrderDisorderCorrectionSUBL
!
!
subroutine ApplyOrderDisorderCorrectionSUBL(iSolnIndex)
!
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer, intent(in) :: iSolnIndex
    integer :: i, j, m, s, c, d, sd, so, mo, md
    integer :: iDisIndex, iOrdPhaseID, iDisPhaseID, nOrdSub, nDisSub
    integer :: iFirstOrd, iLastOrd, iFirstDis, iLastDis, nOrdSpecies, nDisSpecies
    integer :: iBestDisSub, iBestCount, iCountDis, iConOrd, iConDis
    real(8) :: dTemp, dSum, dNorm, dWeight, dMappedG, dMappedH, dMappedS, dMappedCp
    real(8) :: dDisScalarG, dDisScalarH, dDisScalarS, dDisScalarCp
    real(8) :: dRandomScalarG, dRandomScalarH, dRandomScalarS, dRandomScalarCp
    real(8) :: dCorrG, dCorrH, dCorrS, dCorrCp
    real(8) :: dOrdTCrit, dOrdB, dDisTCrit, dDisB, dRandomTCrit, dRandomB
    real(8) :: dPartTCrit, dPartB
    real(8) :: dMappedTCrit, dMappedB, dMappedMagG, dMappedMagH, dMappedMagS, dMappedMagCp
    real(8) :: dOrdMagG, dOrdMagH, dOrdMagS, dOrdMagCp
    real(8) :: dDisMagG, dDisMagH, dDisMagS, dDisMagCp
    real(8) :: dRandomMagG, dRandomMagH, dRandomMagS, dRandomMagCp
    real(8) :: dPartMagG, dPartMagH, dPartMagS, dPartMagCp
    real(8) :: dAvgRelG, dAvgRelH, dAvgRelS, dAvgRelCp, dAvgRelTCrit, dAvgRelB
    real(8) :: dConditionalG, dConditionalH, dConditionalS, dConditionalCp, dConditionalNorm
    real(8) :: dConditionalTCrit, dConditionalB
    logical :: lSubset, lFound, lCorrectionApplied
!
    integer, dimension(nMaxSublatticeSys) :: iOrdToDis
    integer, dimension(nMaxSublatticeSys,nMaxConstituentSys) :: iOrdConToDisCon
    real(8), dimension(nMaxSublatticeSys) :: dGroupStoich
    real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys) :: dCollapsedSiteFraction
    real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys) :: dRelCorrG, dRelCorrH, dRelCorrS, dRelCorrCp
    real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys) :: dRelCorrTCrit, dRelCorrB
!
    real(8), dimension(:), allocatable :: dSaveMolFraction, dSaveChemicalPotential, dSaveGibbsSolnPhase
    real(8), dimension(:), allocatable :: dSavePartialEnthalpy, dSavePartialEntropy, dSavePartialHeatCapacity
    real(8), dimension(:), allocatable :: dSavePartialExcessGibbs, dSavePartialEnthalpyXS
    real(8), dimension(:), allocatable :: dSavePartialEntropyXS, dSavePartialHeatCapacityXS
    real(8), dimension(:), allocatable :: dSaveMagGibbsEnergy, dSaveMagEnthalpy
    real(8), dimension(:), allocatable :: dSaveMagEntropy, dSaveMagHeatCapacity
    real(8), dimension(:), allocatable :: dDisG, dDisH, dDisS, dDisCp
    real(8), dimension(:), allocatable :: dRandomG, dRandomH, dRandomS, dRandomCp
    real(8), dimension(:), allocatable :: dOrdTCritPartial, dOrdBPartial
    real(8), dimension(:), allocatable :: dDisTCritPartial, dDisBPartial
    real(8), dimension(:), allocatable :: dRandomTCritPartial, dRandomBPartial
    real(8), dimension(:,:,:), allocatable :: dSaveSiteFraction
!
    iDisIndex = 0
    if (allocated(iDisorderedPhase)) iDisIndex = iDisorderedPhase(iSolnIndex)
    if (iDisIndex <= 0) return
    if (iDisIndex > nSolnPhasesSys) return
    if (.NOT.allocated(iPhaseSublattice)) return
!
    iOrdPhaseID = iPhaseSublattice(iSolnIndex)
    iDisPhaseID = iPhaseSublattice(iDisIndex)
    if ((iOrdPhaseID <= 0).OR.(iDisPhaseID <= 0)) return
!
    nOrdSub = nSublatticePhase(iOrdPhaseID)
    nDisSub = nSublatticePhase(iDisPhaseID)
    if ((nOrdSub <= 0).OR.(nDisSub <= 0)) return
!
    iFirstOrd = nSpeciesPhase(iSolnIndex-1) + 1
    iLastOrd  = nSpeciesPhase(iSolnIndex)
    iFirstDis = nSpeciesPhase(iDisIndex-1) + 1
    iLastDis  = nSpeciesPhase(iDisIndex)
    nOrdSpecies = iLastOrd - iFirstOrd + 1
    nDisSpecies = iLastDis - iFirstDis + 1
    if ((nOrdSpecies <= 0).OR.(nDisSpecies <= 0)) return
    lCorrectionApplied = .FALSE.
!
    iOrdToDis = 0
    iOrdConToDisCon = 0
    dGroupStoich = 0D0
    dCollapsedSiteFraction = 0D0
    dRelCorrG = 0D0
    dRelCorrH = 0D0
    dRelCorrS = 0D0
    dRelCorrCp = 0D0
    dRelCorrTCrit = 0D0
    dRelCorrB = 0D0
    dDisScalarG = 0D0
    dDisScalarH = 0D0
    dDisScalarS = 0D0
    dDisScalarCp = 0D0
    dRandomScalarG = 0D0
    dRandomScalarH = 0D0
    dRandomScalarS = 0D0
    dRandomScalarCp = 0D0
    dOrdTCrit = 0D0
    dOrdB = 0D0
    dDisTCrit = 0D0
    dDisB = 0D0
    dRandomTCrit = 0D0
    dRandomB = 0D0
    dPartTCrit = 0D0
    dPartB = 0D0
    dOrdMagG = 0D0
    dOrdMagH = 0D0
    dOrdMagS = 0D0
    dOrdMagCp = 0D0
    dDisMagG = 0D0
    dDisMagH = 0D0
    dDisMagS = 0D0
    dDisMagCp = 0D0
    dRandomMagG = 0D0
    dRandomMagH = 0D0
    dRandomMagS = 0D0
    dRandomMagCp = 0D0
    dPartMagG = 0D0
    dPartMagH = 0D0
    dPartMagS = 0D0
    dPartMagCp = 0D0
!
    ! Map every ordered sublattice to the smallest compatible disordered
    ! sublattice.  This resolves cases where Va appears on both a
    ! substitutional sublattice and a pure vacancy sublattice.
    LOOP_MapOrdSub: do so = 1, nOrdSub
        iBestDisSub = 0
        iBestCount = nMaxConstituentSys + 1
!
        do sd = 1, nDisSub
            lSubset = .TRUE.
            do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                lFound = .FALSE.
                do d = 1, nConstituentSublattice(iDisPhaseID,sd)
                    if (TRIM(cConstituentNameSUB(iOrdPhaseID,so,c)) == &
                        TRIM(cConstituentNameSUB(iDisPhaseID,sd,d))) then
                        lFound = .TRUE.
                        exit
                    end if
                end do
                if (.NOT.lFound) then
                    lSubset = .FALSE.
                    exit
                end if
            end do
!
            iCountDis = nConstituentSublattice(iDisPhaseID,sd)
            if (lSubset.AND.(iCountDis < iBestCount)) then
                iBestDisSub = sd
                iBestCount = iCountDis
            end if
        end do
!
        if (iBestDisSub == 0) return
        iOrdToDis(so) = iBestDisSub
        dGroupStoich(iBestDisSub) = dGroupStoich(iBestDisSub) + dStoichSublattice(iOrdPhaseID,so)
!
        do c = 1, nConstituentSublattice(iOrdPhaseID,so)
            do d = 1, nConstituentSublattice(iDisPhaseID,iBestDisSub)
                if (TRIM(cConstituentNameSUB(iOrdPhaseID,so,c)) == &
                    TRIM(cConstituentNameSUB(iDisPhaseID,iBestDisSub,d))) then
                    iOrdConToDisCon(so,c) = d
                    exit
                end if
            end do
            if (iOrdConToDisCon(so,c) == 0) return
        end do
    end do LOOP_MapOrdSub
!
    do sd = 1, nDisSub
        if (dGroupStoich(sd) <= 0D0) return
    end do
!
    ! Collapse the ordered site fractions into the disordered phase site
    ! fractions using the ordered sublattice stoichiometries.
    do so = 1, nOrdSub
        dSum = 0D0
        do c = 1, nConstituentSublattice(iOrdPhaseID,so)
            dSum = dSum + dSiteFraction(iOrdPhaseID,so,c)
        end do
        if (dSum <= 0D0) return
!
        sd = iOrdToDis(so)
        dWeight = dStoichSublattice(iOrdPhaseID,so) / dGroupStoich(sd)
        do c = 1, nConstituentSublattice(iOrdPhaseID,so)
            d = iOrdConToDisCon(so,c)
            dCollapsedSiteFraction(sd,d) = dCollapsedSiteFraction(sd,d) + &
                dWeight * dSiteFraction(iOrdPhaseID,so,c) / dSum
        end do
    end do
!
    do sd = 1, nDisSub
        dNorm = 0D0
        do d = 1, nConstituentSublattice(iDisPhaseID,sd)
            dNorm = dNorm + dCollapsedSiteFraction(sd,d)
        end do
        if (dNorm <= 0D0) return
        do d = 1, nConstituentSublattice(iDisPhaseID,sd)
            dCollapsedSiteFraction(sd,d) = dCollapsedSiteFraction(sd,d) / dNorm
        end do
    end do
!
    allocate(dSaveMolFraction(nSpecies), dSaveChemicalPotential(nSpecies), dSaveGibbsSolnPhase(nSolnPhasesSys), &
        dSavePartialEnthalpy(nSpecies), dSavePartialEntropy(nSpecies), dSavePartialHeatCapacity(nSpecies), &
        dSavePartialExcessGibbs(nSpecies), dSavePartialEnthalpyXS(nSpecies), dSavePartialEntropyXS(nSpecies), &
        dSavePartialHeatCapacityXS(nSpecies), dSaveMagGibbsEnergy(nSpecies), dSaveMagEnthalpy(nSpecies), &
        dSaveMagEntropy(nSpecies), dSaveMagHeatCapacity(nSpecies))
    allocate(dDisG(nDisSpecies), dDisH(nDisSpecies), dDisS(nDisSpecies), dDisCp(nDisSpecies))
    allocate(dRandomG(nOrdSpecies), dRandomH(nOrdSpecies), dRandomS(nOrdSpecies), dRandomCp(nOrdSpecies))
    allocate(dOrdTCritPartial(nOrdSpecies), dOrdBPartial(nOrdSpecies), &
        dDisTCritPartial(nDisSpecies), dDisBPartial(nDisSpecies), &
        dRandomTCritPartial(nOrdSpecies), dRandomBPartial(nOrdSpecies))
    allocate(dSaveSiteFraction(SIZE(dSiteFraction,1),SIZE(dSiteFraction,2),SIZE(dSiteFraction,3)))
!
    dDisG = 0D0
    dDisH = 0D0
    dDisS = 0D0
    dDisCp = 0D0
    dRandomG = 0D0
    dRandomH = 0D0
    dRandomS = 0D0
    dRandomCp = 0D0
    dOrdTCritPartial = 0D0
    dOrdBPartial = 0D0
    dDisTCritPartial = 0D0
    dDisBPartial = 0D0
    dRandomTCritPartial = 0D0
    dRandomBPartial = 0D0
!
    dSaveMolFraction = dMolFraction
    dSaveChemicalPotential = dChemicalPotential
    dSaveGibbsSolnPhase = dGibbsSolnPhase
    dSavePartialEnthalpy = dPartialEnthalpy
    dSavePartialEntropy = dPartialEntropy
    dSavePartialHeatCapacity = dPartialHeatCapacity
    dSavePartialExcessGibbs = dPartialExcessGibbs
    dSavePartialEnthalpyXS = dPartialEnthalpyXS
    dSavePartialEntropyXS = dPartialEntropyXS
    dSavePartialHeatCapacityXS = dPartialHeatCapacityXS
    dSaveMagGibbsEnergy = dMagGibbsEnergy
    dSaveMagEnthalpy = dMagEnthalpy
    dSaveMagEntropy = dMagEntropy
    dSaveMagHeatCapacity = dMagHeatCapacity
    dSaveSiteFraction = dSiteFraction
    call CompMagneticTemperatureMoment(iSolnIndex, dOrdTCrit, dOrdB)
    call CompMagneticVariablePartialsSUBL(iSolnIndex, dOrdTCritPartial, dOrdBPartial)
    call EvaluateMagneticScalarSUBL(iSolnIndex, dOrdTCrit, dOrdB, &
        dOrdMagG, dOrdMagH, dOrdMagS, dOrdMagCp)
    lOrderDisorderEvaluation = .TRUE.
!
    ! Evaluate the disordered phase at the collapsed state.  The partitioned
    ! CEF correction is G_dis(collapsed) - G_ord(random); after this correction
    ! is added to the raw ordered expression, the ordering contribution vanishes
    ! at the random state and the ordered phase recovers the disordered Gibbs
    ! energy.
    !
    ! Magnetic SUBOM phases partition TC/BMAGN before evaluating the
    ! Hillert-Jarl magnetic function.  Value-level magnetic Gibbs subtraction is
    ! not equivalent because the magnetic Gibbs function is nonlinear in those
    ! variables.
    dMolFraction(iFirstDis:iLastDis) = 0D0
    do i = iFirstDis, iLastDis
        m = i - iFirstDis + 1
        dTemp = 1D0
        do sd = 1, nDisSub
            d = iConstituentSublattice(iDisPhaseID,sd,m)
            dTemp = dTemp * dCollapsedSiteFraction(sd,d)
        end do
        dMolFraction(i) = dTemp
    end do
    dNorm = SUM(dMolFraction(iFirstDis:iLastDis))
    if (dNorm > 0D0) dMolFraction(iFirstDis:iLastDis) = dMolFraction(iFirstDis:iLastDis) / dNorm
!
    dChemicalPotential(iFirstDis:iLastDis) = 0D0
    dPartialEnthalpy(iFirstDis:iLastDis) = 0D0
    dPartialEntropy(iFirstDis:iLastDis) = 0D0
    dPartialHeatCapacity(iFirstDis:iLastDis) = 0D0
    dPartialExcessGibbs(iFirstDis:iLastDis) = 0D0
    dPartialEnthalpyXS(iFirstDis:iLastDis) = 0D0
    dPartialEntropyXS(iFirstDis:iLastDis) = 0D0
    dPartialHeatCapacityXS(iFirstDis:iLastDis) = 0D0
    dMagGibbsEnergy(iFirstDis:iLastDis) = 0D0
    dMagEnthalpy(iFirstDis:iLastDis) = 0D0
    dMagEntropy(iFirstDis:iLastDis) = 0D0
    dMagHeatCapacity(iFirstDis:iLastDis) = 0D0
    dGibbsSolnPhase(iDisIndex) = 0D0
!
    call CompExcessGibbsEnergySUBL(iDisIndex)
    if (INFOThermo /= 0) goto 900
    call CompMagneticTemperatureMoment(iDisIndex, dDisTCrit, dDisB)
    call CompMagneticVariablePartialsSUBL(iDisIndex, dDisTCritPartial, dDisBPartial)
    call EvaluateMagneticScalarSUBL(iDisIndex, dDisTCrit, dDisB, &
        dDisMagG, dDisMagH, dDisMagS, dDisMagCp)
!
    do i = iFirstDis, iLastDis
        m = i - iFirstDis + 1
        dDisG(m) = dChemicalPotential(i) + dPartialExcessGibbs(i)
        dDisH(m) = dPartialEnthalpy(i) + dPartialEnthalpyXS(i)
        dDisS(m) = dPartialEntropy(i) + dPartialEntropyXS(i)
        dDisCp(m) = dPartialHeatCapacity(i) + dPartialHeatCapacityXS(i)
        dDisScalarG = dDisScalarG + dMolFraction(i) * dDisG(m)
        dDisScalarH = dDisScalarH + dMolFraction(i) * dDisH(m)
        dDisScalarS = dDisScalarS + dMolFraction(i) * dDisS(m)
        dDisScalarCp = dDisScalarCp + dMolFraction(i) * dDisCp(m)
    end do
!
    ! Evaluate the ordered phase at the random disordered state.
    dMolFraction(iFirstOrd:iLastOrd) = 0D0
    do i = iFirstOrd, iLastOrd
        m = i - iFirstOrd + 1
        dTemp = 1D0
        do so = 1, nOrdSub
            c = iConstituentSublattice(iOrdPhaseID,so,m)
            sd = iOrdToDis(so)
            d = iOrdConToDisCon(so,c)
            dTemp = dTemp * dCollapsedSiteFraction(sd,d)
        end do
        dMolFraction(i) = dTemp
    end do
    dNorm = SUM(dMolFraction(iFirstOrd:iLastOrd))
    if (dNorm > 0D0) dMolFraction(iFirstOrd:iLastOrd) = dMolFraction(iFirstOrd:iLastOrd) / dNorm
!
    dChemicalPotential(iFirstOrd:iLastOrd) = 0D0
    dPartialEnthalpy(iFirstOrd:iLastOrd) = 0D0
    dPartialEntropy(iFirstOrd:iLastOrd) = 0D0
    dPartialHeatCapacity(iFirstOrd:iLastOrd) = 0D0
    dPartialExcessGibbs(iFirstOrd:iLastOrd) = 0D0
    dPartialEnthalpyXS(iFirstOrd:iLastOrd) = 0D0
    dPartialEntropyXS(iFirstOrd:iLastOrd) = 0D0
    dPartialHeatCapacityXS(iFirstOrd:iLastOrd) = 0D0
    dMagGibbsEnergy(iFirstOrd:iLastOrd) = 0D0
    dMagEnthalpy(iFirstOrd:iLastOrd) = 0D0
    dMagEntropy(iFirstOrd:iLastOrd) = 0D0
    dMagHeatCapacity(iFirstOrd:iLastOrd) = 0D0
    dGibbsSolnPhase(iSolnIndex) = 0D0
!
    call CompExcessGibbsEnergySUBL(iSolnIndex)
    if (INFOThermo /= 0) goto 900
    call CompMagneticTemperatureMoment(iSolnIndex, dRandomTCrit, dRandomB)
    call CompMagneticVariablePartialsSUBL(iSolnIndex, dRandomTCritPartial, dRandomBPartial)
    call EvaluateMagneticScalarSUBL(iSolnIndex, dRandomTCrit, dRandomB, &
        dRandomMagG, dRandomMagH, dRandomMagS, dRandomMagCp)
!
    do i = iFirstOrd, iLastOrd
        m = i - iFirstOrd + 1
        dRandomG(m) = dChemicalPotential(i) + dPartialExcessGibbs(i)
        dRandomH(m) = dPartialEnthalpy(i) + dPartialEnthalpyXS(i)
        dRandomS(m) = dPartialEntropy(i) + dPartialEntropyXS(i)
        dRandomCp(m) = dPartialHeatCapacity(i) + dPartialHeatCapacityXS(i)
        dRandomScalarG = dRandomScalarG + dMolFraction(i) * dRandomG(m)
        dRandomScalarH = dRandomScalarH + dMolFraction(i) * dRandomH(m)
        dRandomScalarS = dRandomScalarS + dMolFraction(i) * dRandomS(m)
        dRandomScalarCp = dRandomScalarCp + dMolFraction(i) * dRandomCp(m)
    end do
!
    dPartTCrit = dOrdTCrit + dDisTCrit - dRandomTCrit
    dPartB = dOrdB + dDisB - dRandomB
    call EvaluateMagneticScalarSUBL(iSolnIndex, dPartTCrit, dPartB, &
        dPartMagG, dPartMagH, dPartMagS, dPartMagCp)
!
    dCorrG = dDisScalarG - dRandomScalarG + dPartMagG - dOrdMagG
    dCorrH = dDisScalarH - dRandomScalarH + dPartMagH - dOrdMagH
    dCorrS = dDisScalarS - dRandomScalarS + dPartMagS - dOrdMagS
    dCorrCp = dDisScalarCp - dRandomScalarCp + dPartMagCp - dOrdMagCp
!
    ! Convert the disordered-minus-random-ordered scalar correction into
    ! relative derivatives with respect to the collapsed disordered site
    ! fractions.  The relative form is invariant to the arbitrary additive
    ! constant in site-fraction gradients.
    do sd = 1, nDisSub
        if (nConstituentSublattice(iDisPhaseID,sd) <= 1) cycle
!
        do d = 1, nConstituentSublattice(iDisPhaseID,sd)
            dConditionalG = 0D0
            dConditionalH = 0D0
            dConditionalS = 0D0
            dConditionalCp = 0D0
            dConditionalTCrit = 0D0
            dConditionalB = 0D0
            dConditionalNorm = 0D0
!
            do j = iFirstDis, iLastDis
                md = j - iFirstDis + 1
                if (iConstituentSublattice(iDisPhaseID,sd,md) /= d) cycle
!
                dWeight = 1D0
                do s = 1, nDisSub
                    if (s == sd) cycle
                    c = iConstituentSublattice(iDisPhaseID,s,md)
                    dWeight = dWeight * dCollapsedSiteFraction(s,c)
                end do
!
                dConditionalG = dConditionalG + dWeight * dDisG(md)
                dConditionalH = dConditionalH + dWeight * dDisH(md)
                dConditionalS = dConditionalS + dWeight * dDisS(md)
                dConditionalCp = dConditionalCp + dWeight * dDisCp(md)
                dConditionalTCrit = dConditionalTCrit + dWeight * dDisTCritPartial(md)
                dConditionalB = dConditionalB + dWeight * dDisBPartial(md)
                dConditionalNorm = dConditionalNorm + dWeight
            end do
!
            if (dConditionalNorm > 0D0) then
                dRelCorrG(sd,d) = dRelCorrG(sd,d) + dConditionalG / dConditionalNorm - dDisScalarG
                dRelCorrH(sd,d) = dRelCorrH(sd,d) + dConditionalH / dConditionalNorm - dDisScalarH
                dRelCorrS(sd,d) = dRelCorrS(sd,d) + dConditionalS / dConditionalNorm - dDisScalarS
                dRelCorrCp(sd,d) = dRelCorrCp(sd,d) + dConditionalCp / dConditionalNorm - dDisScalarCp
                dRelCorrTCrit(sd,d) = dRelCorrTCrit(sd,d) + dConditionalTCrit / dConditionalNorm - dDisTCrit
                dRelCorrB(sd,d) = dRelCorrB(sd,d) + dConditionalB / dConditionalNorm - dDisB
            end if
        end do
!
        do so = 1, nOrdSub
            if (iOrdToDis(so) /= sd) cycle
!
            do d = 1, nConstituentSublattice(iDisPhaseID,sd)
                iConOrd = 0
                do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                    if (iOrdConToDisCon(so,c) == d) then
                        iConOrd = c
                        exit
                    end if
                end do
                if (iConOrd == 0) cycle
!
                dConditionalG = 0D0
                dConditionalH = 0D0
                dConditionalS = 0D0
                dConditionalCp = 0D0
                dConditionalTCrit = 0D0
                dConditionalB = 0D0
                dConditionalNorm = 0D0
!
                do j = iFirstOrd, iLastOrd
                    md = j - iFirstOrd + 1
                    c = iConstituentSublattice(iOrdPhaseID,so,md)
                    if (c /= iConOrd) cycle
!
                    dWeight = 1D0
                    do s = 1, nOrdSub
                        if (s == so) cycle
                        c = iConstituentSublattice(iOrdPhaseID,s,md)
                        dWeight = dWeight * dCollapsedSiteFraction(iOrdToDis(s),iOrdConToDisCon(s,c))
                    end do
!
                    dConditionalG = dConditionalG + dWeight * dRandomG(md)
                    dConditionalH = dConditionalH + dWeight * dRandomH(md)
                    dConditionalS = dConditionalS + dWeight * dRandomS(md)
                    dConditionalCp = dConditionalCp + dWeight * dRandomCp(md)
                    dConditionalTCrit = dConditionalTCrit + dWeight * dRandomTCritPartial(md)
                    dConditionalB = dConditionalB + dWeight * dRandomBPartial(md)
                    dConditionalNorm = dConditionalNorm + dWeight
                end do
!
                if (dConditionalNorm > 0D0) then
                    dRelCorrG(sd,d) = dRelCorrG(sd,d) - (dConditionalG / dConditionalNorm - dRandomScalarG)
                    dRelCorrH(sd,d) = dRelCorrH(sd,d) - (dConditionalH / dConditionalNorm - dRandomScalarH)
                    dRelCorrS(sd,d) = dRelCorrS(sd,d) - (dConditionalS / dConditionalNorm - dRandomScalarS)
                    dRelCorrCp(sd,d) = dRelCorrCp(sd,d) - (dConditionalCp / dConditionalNorm - dRandomScalarCp)
                    dRelCorrTCrit(sd,d) = dRelCorrTCrit(sd,d) - &
                        (dConditionalTCrit / dConditionalNorm - dRandomTCrit)
                    dRelCorrB(sd,d) = dRelCorrB(sd,d) - (dConditionalB / dConditionalNorm - dRandomB)
                end if
            end do
        end do
    end do
!
    ! Restore the caller's ordered-phase state and add the CEF partial-molar
    ! transform of the analytic disordered-minus-random-ordered correction.
    dMolFraction = dSaveMolFraction
    dSiteFraction = dSaveSiteFraction
    dChemicalPotential = dSaveChemicalPotential
    dGibbsSolnPhase = dSaveGibbsSolnPhase
    dPartialEnthalpy = dSavePartialEnthalpy
    dPartialEntropy = dSavePartialEntropy
    dPartialHeatCapacity = dSavePartialHeatCapacity
    dPartialExcessGibbs = dSavePartialExcessGibbs
    dPartialEnthalpyXS = dSavePartialEnthalpyXS
    dPartialEntropyXS = dSavePartialEntropyXS
    dPartialHeatCapacityXS = dSavePartialHeatCapacityXS
    dMagGibbsEnergy = dSaveMagGibbsEnergy
    dMagEnthalpy = dSaveMagEnthalpy
    dMagEntropy = dSaveMagEntropy
    dMagHeatCapacity = dSaveMagHeatCapacity
!
    do i = iFirstOrd, iLastOrd
        mo = i - iFirstOrd + 1
        dMappedG = dCorrG
        dMappedH = dCorrH
        dMappedS = dCorrS
        dMappedCp = dCorrCp
        dMappedTCrit = dPartTCrit + dOrdTCritPartial(mo) - dOrdTCrit
        dMappedB = dPartB + dOrdBPartial(mo) - dOrdB
!
        do so = 1, nOrdSub
            sd = iOrdToDis(so)
            if (nConstituentSublattice(iDisPhaseID,sd) <= 1) cycle
!
            iConOrd = iConstituentSublattice(iOrdPhaseID,so,mo)
            iConDis = iOrdConToDisCon(so,iConOrd)
            dTemp = dStoichSublattice(iOrdPhaseID,so) / dGroupStoich(sd)
!
            dAvgRelG = 0D0
            dAvgRelH = 0D0
            dAvgRelS = 0D0
            dAvgRelCp = 0D0
            dAvgRelTCrit = 0D0
            dAvgRelB = 0D0
            dNorm = 0D0
!
            do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                d = iOrdConToDisCon(so,c)
                dWeight = dSiteFraction(iOrdPhaseID,so,c)
                dAvgRelG = dAvgRelG + dWeight * dRelCorrG(sd,d)
                dAvgRelH = dAvgRelH + dWeight * dRelCorrH(sd,d)
                dAvgRelS = dAvgRelS + dWeight * dRelCorrS(sd,d)
                dAvgRelCp = dAvgRelCp + dWeight * dRelCorrCp(sd,d)
                dAvgRelTCrit = dAvgRelTCrit + dWeight * dRelCorrTCrit(sd,d)
                dAvgRelB = dAvgRelB + dWeight * dRelCorrB(sd,d)
                dNorm = dNorm + dWeight
            end do
!
            if (dNorm > 0D0) then
                dAvgRelG = dAvgRelG / dNorm
                dAvgRelH = dAvgRelH / dNorm
                dAvgRelS = dAvgRelS / dNorm
                dAvgRelCp = dAvgRelCp / dNorm
                dAvgRelTCrit = dAvgRelTCrit / dNorm
                dAvgRelB = dAvgRelB / dNorm
            end if
!
            dMappedG = dMappedG + dTemp * (dRelCorrG(sd,iConDis) - dAvgRelG)
            dMappedH = dMappedH + dTemp * (dRelCorrH(sd,iConDis) - dAvgRelH)
            dMappedS = dMappedS + dTemp * (dRelCorrS(sd,iConDis) - dAvgRelS)
            dMappedCp = dMappedCp + dTemp * (dRelCorrCp(sd,iConDis) - dAvgRelCp)
            dMappedTCrit = dMappedTCrit + dTemp * (dRelCorrTCrit(sd,iConDis) - dAvgRelTCrit)
            dMappedB = dMappedB + dTemp * (dRelCorrB(sd,iConDis) - dAvgRelB)
        end do
!
        call EvaluateMagneticPartialSUBL(iSolnIndex, dPartTCrit, dPartB, dMappedTCrit, dMappedB, &
            dMappedMagG, dMappedMagH, dMappedMagS, dMappedMagCp)
        dMappedG = dMappedG + (dMappedMagG - dPartMagG) - (dSaveMagGibbsEnergy(i) - dOrdMagG)
        dMappedH = dMappedH + (dMappedMagH - dPartMagH) - (dSaveMagEnthalpy(i) - dOrdMagH)
        dMappedS = dMappedS + (dMappedMagS - dPartMagS) - (dSaveMagEntropy(i) - dOrdMagS)
        dMappedCp = dMappedCp + (dMappedMagCp - dPartMagCp) - &
            (dSaveMagHeatCapacity(i) - dOrdMagCp)
        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) + dMappedG
        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i) + dMappedH
        dPartialEntropyXS(i) = dPartialEntropyXS(i) + dMappedS
        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i) + dMappedCp
    end do
    lCorrectionApplied = .TRUE.
!
900 continue
    lOrderDisorderEvaluation = .FALSE.
    dMolFraction = dSaveMolFraction
    dSiteFraction = dSaveSiteFraction
    dChemicalPotential = dSaveChemicalPotential
    dGibbsSolnPhase = dSaveGibbsSolnPhase
    dPartialEnthalpy = dSavePartialEnthalpy
    dPartialEntropy = dSavePartialEntropy
    dPartialHeatCapacity = dSavePartialHeatCapacity
    if (.NOT.lCorrectionApplied) then
        dPartialExcessGibbs = dSavePartialExcessGibbs
        dPartialEnthalpyXS = dSavePartialEnthalpyXS
        dPartialEntropyXS = dSavePartialEntropyXS
        dPartialHeatCapacityXS = dSavePartialHeatCapacityXS
    end if
    dMagGibbsEnergy = dSaveMagGibbsEnergy
    dMagEnthalpy = dSaveMagEnthalpy
    dMagEntropy = dSaveMagEntropy
    dMagHeatCapacity = dSaveMagHeatCapacity
!
    deallocate(dSaveMolFraction, dSaveChemicalPotential, dSaveGibbsSolnPhase, dSavePartialEnthalpy, dSavePartialEntropy, &
        dSavePartialHeatCapacity, dSavePartialExcessGibbs, dSavePartialEnthalpyXS, dSavePartialEntropyXS, &
        dSavePartialHeatCapacityXS, dSaveMagGibbsEnergy, dSaveMagEnthalpy, dSaveMagEntropy, &
        dSaveMagHeatCapacity, dDisG, dDisH, dDisS, dDisCp, dRandomG, dRandomH, dRandomS, dRandomCp, &
        dOrdTCritPartial, dOrdBPartial, dDisTCritPartial, dDisBPartial, dRandomTCritPartial, &
        dRandomBPartial, dSaveSiteFraction)
!
    return
!
contains
!
    subroutine EvaluateMagneticScalarSUBL(iSolnIndexIn, dTCritIn, dBIn, dMagG, dMagH, dMagS, dMagCp)
!
    implicit none
!
    integer, intent(in) :: iSolnIndexIn
    real(8), intent(in) :: dTCritIn, dBIn
    real(8), intent(out) :: dMagG, dMagH, dMagS, dMagCp
!
    integer :: iFirstLocal
    real(8) :: dTCritLocal, dBLocal, dStructureFactor, dP, dInvPMinusOne
    real(8) :: dTau, dD, dA, dBTemp, dC, dF, dFp, dFdp
    real(8) :: dLogMoment
!
    dMagG = 0D0
    dMagH = 0D0
    dMagS = 0D0
    dMagCp = 0D0
    if (iSolnIndexIn <= 0) return
!
    iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
    dStructureFactor = dCoeffGibbsMagnetic(iFirstLocal,3)
    dP = dCoeffGibbsMagnetic(iFirstLocal,4)
    if (dP == 0D0) return
!
    dTCritLocal = dTCritIn
    dBLocal = dBIn
    if (dBLocal < 0D0) dBLocal = -dBLocal * dStructureFactor
    if (dTCritLocal < 0D0) dTCritLocal = -dTCritLocal * dStructureFactor
    if ((dTCritLocal == 0D0).OR.(dBLocal <= -1D0)) return
!
    dTau = dTemperature / dTCritLocal
    dInvPMinusOne = 1D0 / dP - 1D0
    dD = (518D0/1125D0) + (11692D0/15975D0) * dInvPMinusOne
!
    if (dTau > 1D0) then
        dA = dTau**(-5)
        dBTemp = dA**3
        dC = dA * dA * dBTemp
!
        dF = -(dA/10D0 + dBTemp/315D0 + dC/1500D0) / dD
        dFp = (1D0 / (dD * dTCritLocal)) * (dA/2D0 + dBTemp/21D0 + dC/60D0)
        dFdp = -(1D0 / (dD * dTCritLocal)) * (3D0*dA + 16D0*dBTemp/21D0 + 13D0*dC/30D0)
    else
        dA = dTau**3
        dBTemp = dA**3
        dC = dA * dA * dBTemp
!
        dF = 1D0 - (79D0/(140D0*dP*dTau) + &
            (474D0/497D0)*dInvPMinusOne*(dA/6D0 + dBTemp/135D0 + dC/600D0)) / dD
        dFp = (79D0/(140D0*dP*dTemperature) - &
            (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA/2D0 + dBTemp/15D0 + dC/40D0)) / dD
        dFdp = -(79D0/(70D0*dP*dTemperature) + &
            (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA + 8D0*dBTemp/15D0 + 7D0*dC/20D0)) / dD
    end if
!
    dLogMoment = DLOG(1D0 + dBLocal)
    dMagG = dLogMoment * dF
    dMagH = -dLogMoment * dTCritLocal * dFp
    dMagS = -dLogMoment * (dF + dTCritLocal * dFp)
    dMagCp = -dLogMoment * dTCritLocal * (2D0*dFp + dFdp)
!
    return
!
    end subroutine EvaluateMagneticScalarSUBL
!
    subroutine EvaluateMagneticPartialSUBL(iSolnIndexIn, dTCritIn, dBIn, &
        dTCritPartialIn, dBPartialIn, dMagG, dMagH, dMagS, dMagCp)
!
    implicit none
!
    integer, intent(in) :: iSolnIndexIn
    real(8), intent(in) :: dTCritIn, dBIn, dTCritPartialIn, dBPartialIn
    real(8), intent(out) :: dMagG, dMagH, dMagS, dMagCp
!
    integer :: iFirstLocal
    real(8) :: dTCritLocal, dBLocal, dTCritPartialLocal, dBPartialLocal
    real(8) :: dStructureFactor, dP, dInvPMinusOne
    real(8) :: dTau, dD, dA, dBTemp, dC, dF, dFp, dFdp, dFtp
    real(8) :: dQp, dQdp, dQtp, dParam1, dParam2, dParam3
!
    dMagG = 0D0
    dMagH = 0D0
    dMagS = 0D0
    dMagCp = 0D0
    if (iSolnIndexIn <= 0) return
!
    iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
    dStructureFactor = dCoeffGibbsMagnetic(iFirstLocal,3)
    dP = dCoeffGibbsMagnetic(iFirstLocal,4)
    if (dP == 0D0) return
!
    dTCritLocal = dTCritIn
    dBLocal = dBIn
    dTCritPartialLocal = dTCritPartialIn
    dBPartialLocal = dBPartialIn
    if (dBLocal < 0D0) then
        dBLocal = -dBLocal * dStructureFactor
        dBPartialLocal = -dBPartialLocal * dStructureFactor
    end if
    if (dTCritLocal < 0D0) then
        dTCritLocal = -dTCritLocal * dStructureFactor
        dTCritPartialLocal = -dTCritPartialLocal * dStructureFactor
    end if
    if ((dTCritLocal == 0D0).OR.(dBLocal <= -1D0)) return
!
    dTau = dTemperature / dTCritLocal
    dInvPMinusOne = 1D0 / dP - 1D0
    dD = (518D0/1125D0) + (11692D0/15975D0) * dInvPMinusOne
!
    if (dTau > 1D0) then
        dA = dTau**(-5)
        dBTemp = dA**3
        dC = dA * dA * dBTemp
!
        dF = -(dA/10D0 + dBTemp/315D0 + dC/1500D0) / dD
        dFp = (1D0 / (dD * dTCritLocal)) * (dA/2D0 + dBTemp/21D0 + dC/60D0)
        dFdp = -(1D0 / (dD * dTCritLocal)) * &
            (3D0*dA + 16D0*dBTemp/21D0 + 13D0*dC/30D0)
        dFtp = +(1D0 / (dD * dTCritLocal)) * &
            (21D0*dA + 272D0*dBTemp/21D0 + 117D0*dC/10D0)
        dQp = dFp
        dQdp = dFdp
        dQtp = dFtp
    else
        dA = dTau**3
        dBTemp = dA**3
        dC = dA * dA * dBTemp
!
        dF = 1D0 - (79D0/(140D0*dP*dTau) + &
            (474D0/497D0)*dInvPMinusOne*(dA/6D0 + dBTemp/135D0 + dC/600D0)) / dD
        dFp = (79D0/(140D0*dP*dTemperature) - &
            (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA/2D0 + dBTemp/15D0 + dC/40D0)) / dD
        dFdp = -(79D0/(70D0*dP*dTemperature) + &
            (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA + 8D0*dBTemp/15D0 + 7D0*dC/20D0)) / dD
        dFtp = (237D0/(70D0*dP*dTemperature) - &
            (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA + 56D0*dBTemp/15D0 + &
            91D0*dC/20D0)) / dD
        dQp = -dFp
        dQdp = -dFdp
        dQtp = -dFtp
    end if
!
    dParam1 = (dBPartialLocal - dBLocal) / (1D0 + dBLocal)
    dParam2 = DLOG(1D0 + dBLocal)
    dParam3 = -(dTCritLocal - dTCritPartialLocal)
!
    dMagG = dF * dParam1 + dParam2 * (dParam3 * dQp + dF)
    dMagH = -(dTCritLocal * dParam1 * dFp + &
        dParam2 * (dTCritLocal * dFp + dParam3 * (dQp + dQdp)))
    dMagS = -dParam1 * (dF + dTCritLocal * dFp) - &
        dParam2 * (dF + dTCritLocal * dFp + dParam3 * (2D0*dQp + dQdp))
    dMagCp = -dParam1 * dTCritLocal * (2D0*dFp + dFdp) - &
        dParam2 * (dTCritLocal * (2D0*dFp + dFdp) + &
        dParam3 * (2D0*dQp + 4D0*dQdp + dQtp))
!
    return
!
    end subroutine EvaluateMagneticPartialSUBL
!
    subroutine CompMagneticVariablePartialsSUBL(iSolnIndexIn, dTcritPartial, dBPartial)
!
    implicit none
!
    integer, intent(in) :: iSolnIndexIn
    real(8), dimension(:), intent(out) :: dTcritPartial, dBPartial
!
    integer :: ii, jj, iFirstLocal, iLastLocal, iParam, iExponent
    integer :: iPhaseIDLocal, nSublatticeLocal, mm, nn, ss, cc, kk, ll, nParamCon
    integer :: iFirstParamLocal, iSecondParamLocal, iSubParamLocal, iSubLocal
    integer :: iKD
    real(8) :: dPreFactor, dPreFactorT, dPreFactorB, dTempLocal, dX1, dX2
!
    dTcritPartial = 0D0
    dBPartial = 0D0
    if (iSolnIndexIn <= 0) return
!
    iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
    iLastLocal = nSpeciesPhase(iSolnIndexIn)
    if (iLastLocal < iFirstLocal) return
    iPhaseIDLocal = iPhaseSublattice(iSolnIndexIn)
    if (iPhaseIDLocal <= 0) return
    nSublatticeLocal = nSublatticePhase(iPhaseIDLocal)
!
    do ii = iFirstLocal, iLastLocal
        mm = ii - iFirstLocal + 1
        do jj = iFirstLocal, iLastLocal
            nn = jj - iFirstLocal + 1
!
            dPreFactor = 1D0 - DFLOAT(nSublatticeLocal)
            do ss = 1, nSublatticeLocal
                kk = iConstituentSublattice(iPhaseIDLocal,ss,mm)
                ll = iConstituentSublattice(iPhaseIDLocal,ss,nn)
                if (kk == ll) dPreFactor = dPreFactor + &
                    1D0 / DMAX1(dSiteFraction(iPhaseIDLocal,ss,kk), 1D-75)
            end do
!
            dTcritPartial(nn) = dTcritPartial(nn) + &
                dPreFactor * dMolFraction(ii) * dCoeffGibbsMagnetic(ii,1)
            dBPartial(nn) = dBPartial(nn) + &
                dPreFactor * dMolFraction(ii) * dCoeffGibbsMagnetic(ii,2)
        end do
    end do
!
    LOOP_MagVar_Param: do iParam = nMagParamPhase(iSolnIndexIn-1)+1, nMagParamPhase(iSolnIndexIn)
        dPreFactor = 1D0
        nParamCon = iMagneticParam(iParam,1)
        iExponent = iMagneticParam(iParam,nParamCon+2)
        iFirstParamLocal = 0
        iSecondParamLocal = 0
        iSubParamLocal = 0
        dX1 = 0D0
        dX2 = 0D0
!
        do kk = 2, nParamCon + 1
            cc = MOD(iMagneticParam(iParam,kk), 10000)
            ss = (iMagneticParam(iParam,kk) - cc) / 10000
            dPreFactor = dPreFactor * dSiteFraction(iPhaseIDLocal,ss,cc)
            if (kk == 2) then
                dX1 = dSiteFraction(iPhaseIDLocal,ss,cc)
                iFirstParamLocal = cc
                iSubParamLocal = ss
            else if (kk == 3) then
                dX2 = dSiteFraction(iPhaseIDLocal,ss,cc)
                iSecondParamLocal = cc
            end if
        end do
!
        dPreFactorT = dPreFactor * dMagneticParam(iParam,1) * (dX1 - dX2)**iExponent
        dPreFactorB = dPreFactor * dMagneticParam(iParam,2) * (dX1 - dX2)**iExponent
!
        do ii = iFirstLocal, iLastLocal
            iKD = 0
            mm = ii - iFirstLocal + 1
            ! CEF partial-molar transform base term for this interaction arity.
            dTempLocal = 1D0 - DFLOAT(nParamCon + iExponent)
!
            do ss = 1, nSublatticeLocal
                cc = iConstituentSublattice(iPhaseIDLocal,ss,mm)
                if (ss == iSubParamLocal) then
                    if (cc == iFirstParamLocal) then
                        iKD = 1
                    else if (cc == iSecondParamLocal) then
                        iKD = -1
                    end if
                end if
!
                do jj = 2, nParamCon + 1
                    kk = MOD(iMagneticParam(iParam,jj), 10000)
                    iSubLocal = (iMagneticParam(iParam,jj) - kk) / 10000
                    if (iSubLocal /= ss) cycle
                    if (kk /= cc) cycle
                    dTempLocal = dTempLocal + &
                        1D0 / DMAX1(dSiteFraction(iPhaseIDLocal,ss,cc), 1D-75)
                end do
            end do
!
            if (dX1 /= dX2) then
                dTempLocal = dTempLocal + DFLOAT(iKD * iExponent) / (dX1 - dX2)
            end if
!
            nn = ii - iFirstLocal + 1
            dTcritPartial(nn) = dTcritPartial(nn) + dPreFactorT * dTempLocal
            dBPartial(nn) = dBPartial(nn) + dPreFactorB * dTempLocal
        end do
    end do LOOP_MagVar_Param
!
    return
!
    end subroutine CompMagneticVariablePartialsSUBL
!
end subroutine ApplyOrderDisorderCorrectionSUBL
!
!
