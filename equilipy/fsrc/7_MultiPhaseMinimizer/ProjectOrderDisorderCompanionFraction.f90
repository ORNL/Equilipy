!> \brief Project a disordered companion constitution into its ordered parent.
!!
!! \details Converts the companion's current element composition into the
!! product-fraction representation of the ordered SUBOM parent on its exact
!! disordered manifold.  This is a representation map, not a minimization or
!! an active-set decision.
!
subroutine ProjectOrderDisorderCompanionFraction(iCompanionPhase, iOrderedPhase, &
    nProjectedFraction, dProjectedFraction, lProjected, iProjectionStatus)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ProjectOrderDisorderCompanionFraction.f90
    !> \brief   Map a DIS_PART companion constitution into an ordered SUBOM parent.
    !> \author  S.Y. Kwon
    !> \date    Jul. 16, 2026
    !> \sa      RegisterSUBOMTwoSetCandidateRows.f90
    !> \sa      ReconcileOrderDisorderCandidateRows.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Projected declared order/disorder companion constitutions structurally and rejected unsupported topologies with typed status.
    !
    ! Purpose:
    ! ========
    !
    !> \details The disordered companion and ordered parent use different
    !! constitution coordinates.  This routine computes the companion's
    !! normalized element composition, assigns that composition identically to
    !! every substitutionally equivalent ordered sublattice, and returns the
    !! resulting ordered-parent product fractions.  Fixed vacancy sublattices
    !! remain fixed; mobile vacancies are rejected by a zero projected norm.
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in] iCompanionPhase  Active-system disordered companion index.
    !> \param[in] iOrderedPhase    Active-system ordered SUBOM parent index.
    ! dMolFraction                 Current companion endmember fractions.
    ! dStoichSpecies               Species stoichiometry on the active element basis.
    ! iConstituentSublattice       Ordered-parent endmember topology.
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[in] nProjectedFraction  Number of ordered-parent product fractions.
    !> \param[out] dProjectedFraction Ordered-parent product fractions.
    !> \param[out] lProjected         True only when the structural map succeeds.
    !> \param[out] iProjectionStatus  Typed success, invalid-input, empty-composition, or topology verdict.
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - The map is exact in the declared sublattice topology; it contains no
    !   composition-distance or ordering-degree tolerance.
    ! - The existing positive-normalization floor only guards division by an
    !   absent composition and is not a phase-identity decision.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: OD_PROJECTION_TOPOLOGY_UNSUPPORTED, &
        OD_PROJECTION_SUCCESS, OD_PROJECTION_INVALID_INPUT, &
        OD_PROJECTION_EMPTY_COMPOSITION, OD_PROJECTION_AMBIGUOUS_TOPOLOGY

    implicit none

    integer, intent(in) :: iCompanionPhase, iOrderedPhase, nProjectedFraction
    real(8), dimension(nProjectedFraction), intent(out) :: dProjectedFraction
    logical, intent(out) :: lProjected
    integer, intent(out) :: iProjectionStatus

    integer :: i, j, s, m, iCompanionFirst, iCompanionLast
    integer :: iOrderedFirst, iOrderedLast, iOrderedLocal, iSublPhase, iElement
    integer :: iChar, iCode
    real(8) :: dCompanionTotal, dNorm, dProduct, dX
    real(8), dimension(nElements) :: dCompanionElementFraction
    character(8) :: cConstituent, cElement
    integer :: OrderDisorderProjectionTopologyClass

    dProjectedFraction = 0D0
    lProjected = .FALSE.
    iProjectionStatus = OD_PROJECTION_INVALID_INPUT
    dCompanionElementFraction = 0D0
    if ((iCompanionPhase < 1).OR.(iCompanionPhase > nSolnPhasesSys)) return
    if ((iOrderedPhase < 1).OR.(iOrderedPhase > nSolnPhasesSys)) return
    if (OrderDisorderProjectionTopologyClass(iOrderedPhase) == &
        OD_PROJECTION_TOPOLOGY_UNSUPPORTED) then
        iProjectionStatus = OD_PROJECTION_AMBIGUOUS_TOPOLOGY
        return
    end if

    iCompanionFirst = nSpeciesPhase(iCompanionPhase-1) + 1
    iCompanionLast = nSpeciesPhase(iCompanionPhase)
    do i = iCompanionFirst, iCompanionLast
        dX = DMAX1(dMolFraction(i), 0D0)
        do j = 1, nElements
            dCompanionElementFraction(j) = dCompanionElementFraction(j) + &
                dX * DMAX1(dStoichSpecies(i,j), 0D0)
        end do
    end do

    dCompanionTotal = SUM(dCompanionElementFraction)
    if (dCompanionTotal <= 1D-300) then
        iProjectionStatus = OD_PROJECTION_EMPTY_COMPOSITION
        return
    end if
    dCompanionElementFraction = dCompanionElementFraction / dCompanionTotal

    iSublPhase = iPhaseSublattice(iOrderedPhase)
    if (iSublPhase <= 0) return
    iOrderedFirst = nSpeciesPhase(iOrderedPhase-1) + 1
    iOrderedLast = nSpeciesPhase(iOrderedPhase)
    if (nProjectedFraction /= (iOrderedLast - iOrderedFirst + 1)) return

    do i = iOrderedFirst, iOrderedLast
        iOrderedLocal = i - iOrderedFirst + 1
        dProduct = 1D0
        do s = 1, nSublatticePhase(iSublPhase)
            m = iConstituentSublattice(iSublPhase,s,iOrderedLocal)
            if (m <= 0) cycle
            cConstituent = cConstituentNameSUB(iSublPhase,s,m)
            do iChar = 1, LEN(cConstituent)
                iCode = IACHAR(cConstituent(iChar:iChar))
                if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                    cConstituent(iChar:iChar) = ACHAR(iCode - 32)
                end if
            end do
            if ((TRIM(cConstituent) == 'VA').OR.(TRIM(cConstituent) == 'VACANCY')) then
                if (nConstituentSublattice(iSublPhase,s) > 1) dProduct = 0D0
                cycle
            end if

            iElement = 0
            do j = 1, nElements
                cElement = cElementName(j)
                do iChar = 1, LEN(cElement)
                    iCode = IACHAR(cElement(iChar:iChar))
                    if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                        cElement(iChar:iChar) = ACHAR(iCode - 32)
                    end if
                end do
                if (TRIM(cConstituent) == TRIM(cElement)) then
                    iElement = j
                    exit
                end if
            end do
            if (iElement <= 0) then
                dProduct = 0D0
            else
                dProduct = dProduct * DMAX1(dCompanionElementFraction(iElement), 0D0)
            end if
        end do
        dProjectedFraction(iOrderedLocal) = dProduct
    end do

    dNorm = SUM(dProjectedFraction)
    if (dNorm <= 1D-300) then
        iProjectionStatus = OD_PROJECTION_EMPTY_COMPOSITION
        return
    end if
    dProjectedFraction = dProjectedFraction / dNorm
    lProjected = .TRUE.
    iProjectionStatus = OD_PROJECTION_SUCCESS

    return

end subroutine ProjectOrderDisorderCompanionFraction


integer function OrderDisorderProjectionTopologyClass(iOrderedPhase)
    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: OD_PROJECTION_TOPOLOGY_UNSUPPORTED, &
        OD_PROJECTION_TOPOLOGY_B2, OD_PROJECTION_TOPOLOGY_L1

    implicit none

    integer, intent(in) :: iOrderedPhase
    integer :: iSublPhase, nSublattice, iSeed, iOther, nMatching, nBest
    integer :: iBestSeed, nBestCount, nBestTies
    logical :: lOutsideSupported

    OrderDisorderProjectionTopologyClass = OD_PROJECTION_TOPOLOGY_UNSUPPORTED
    if ((iOrderedPhase < 1).OR.(iOrderedPhase > nSolnPhasesSys)) return
    if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') return
    iSublPhase = iPhaseSublattice(iOrderedPhase)
    if ((iSublPhase < 1).OR.(iSublPhase > nCountSublattice)) return
    nSublattice = nSublatticePhase(iSublPhase)

    ! The supported B2 and L1 families contain exactly one repeated group of
    ! two or four equal-weight, equal-constituent ordering sublattices.
    iBestSeed = 0
    nBestCount = 0
    nBestTies = 0
    do iSeed = 1, nSublattice
        if (nConstituentSublattice(iSublPhase,iSeed) <= 1) cycle
        nMatching = 0
        do iOther = 1, nSublattice
            if (SublatticeTopologyMatches(iSublPhase, iSeed, iOther)) then
                nMatching = nMatching + 1
            end if
        end do
        if (nMatching > nBestCount) then
            iBestSeed = iSeed
            nBestCount = nMatching
            nBestTies = 1
        else if ((nMatching == nBestCount).AND.(nMatching > 1).AND.&
            (.NOT.SublatticeTopologyMatches(iSublPhase, iBestSeed, iSeed))) then
            nBestTies = nBestTies + 1
        end if
    end do
    if ((nBestCount /= 2).AND.(nBestCount /= 4)) return
    if (nBestTies /= 1) return

    lOutsideSupported = .TRUE.
    do iOther = 1, nSublattice
        if (SublatticeTopologyMatches(iSublPhase, iBestSeed, iOther)) cycle
        if (nConstituentSublattice(iSublPhase,iOther) <= 1) cycle
        if (.NOT.SublatticeContainsVacancy(iSublPhase, iOther)) then
            lOutsideSupported = .FALSE.
            exit
        end if
    end do
    if (.NOT.lOutsideSupported) return

    nBest = nBestCount
    if ((nBest == 2).AND.((nSublattice == 2).OR.(nSublattice == 3))) then
        OrderDisorderProjectionTopologyClass = OD_PROJECTION_TOPOLOGY_B2
    else if ((nBest == 4).AND.((nSublattice == 4).OR.(nSublattice == 5))) then
        OrderDisorderProjectionTopologyClass = OD_PROJECTION_TOPOLOGY_L1
    end if

    return

contains

    logical function SublatticeTopologyMatches(iPhaseIn, iFirstSublattice, iSecondSublattice)
        integer, intent(in) :: iPhaseIn, iFirstSublattice, iSecondSublattice
        integer :: iConstituent, nConstituent

        SublatticeTopologyMatches = .FALSE.
        if (dStoichSublattice(iPhaseIn,iFirstSublattice) /= &
            dStoichSublattice(iPhaseIn,iSecondSublattice)) return
        nConstituent = nConstituentSublattice(iPhaseIn,iFirstSublattice)
        if (nConstituent /= nConstituentSublattice(iPhaseIn,iSecondSublattice)) return
        do iConstituent = 1, nConstituent
            if (TRIM(UpperConstituent(cConstituentNameSUB(&
                iPhaseIn,iFirstSublattice,iConstituent))) /= &
                TRIM(UpperConstituent(cConstituentNameSUB(&
                iPhaseIn,iSecondSublattice,iConstituent)))) return
        end do
        SublatticeTopologyMatches = .TRUE.

        return
    end function SublatticeTopologyMatches


    logical function SublatticeContainsVacancy(iPhaseIn, iSublatticeIn)
        integer, intent(in) :: iPhaseIn, iSublatticeIn
        integer :: iConstituent
        character(8) :: cName

        SublatticeContainsVacancy = .FALSE.
        do iConstituent = 1, nConstituentSublattice(iPhaseIn,iSublatticeIn)
            cName = UpperConstituent(cConstituentNameSUB(&
                iPhaseIn,iSublatticeIn,iConstituent))
            if ((TRIM(cName) == 'VA').OR.(TRIM(cName) == 'VACANCY')) then
                SublatticeContainsVacancy = .TRUE.
                return
            end if
        end do

        return
    end function SublatticeContainsVacancy


    character(8) function UpperConstituent(cName)
        character(*), intent(in) :: cName
        integer :: iChar, iCode, nChar

        UpperConstituent = ' '
        nChar = MIN(LEN_TRIM(cName), LEN(UpperConstituent))
        do iChar = 1, nChar
            iCode = IACHAR(cName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                UpperConstituent(iChar:iChar) = ACHAR(iCode - 32)
            else
                UpperConstituent(iChar:iChar) = cName(iChar:iChar)
            end if
        end do

        return
    end function UpperConstituent

end function OrderDisorderProjectionTopologyClass
