!
!
subroutine CompExcessGibbsEnergyRKMP(iSolnIndex)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergyRKMP.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of solution phase constituents in a RKMP
    !!           or RKMPM solution phase.
    !> \author  M.H.A. Piro
    !> \date    January 14, 2013
    !> \sa      CompExcessGibbsEnergy.f90
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/14/2013      M.H.A. Piro         Original code.
    !   05/14/2014      M.H.A. Piro         Fixed bug in computing binary term: cycle if dx = 0, which could
    !                                        cause problems when computing dxvmo = dx**(iExponent-1) when
    !                                        iExponent = 1.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the partial molar excess Gibbs energy of mixing
    !! (dPartialExcessGibbs) of all constituents in a non-ideal solution phase designated as 'RKMP'
    !! (Redlich-Kister-Muggiano-Polynomial) or 'RKMPM' (RKMP-Magnetic) phases.
    !!
    !! The molar excess Gibbs energy of mixing from a binary interaction term \f$ z \f$ is:
    !!
    !! \f$  g_{\lambda,z}^{ex} = x_1 x_2 \sum_{v=0}^{N_z} {^v} L_{1,2} (x_1 - x_2)^v  \f$
    !!
    !! where \f$ x_1 \f$ and \f$ x_2 \f$ are the mole fractions of constituents 1 and 2 in the binary parameter and
    !! \f$ {^v} L_{1,2} \f$ is the \f$ v \f$th order mixing parameter.
    !!
    !! The molar excess Gibbs energy of mixing from a ternary interaction term \f$ z \f$ is:
    !!
    !! \f$  g_{\lambda,z}^{ex} = x_1 x_2 x_3 \left( \frac{1-x_1 - x_2 - x_3}{3} + x_j \right) {^0}L_{1,2,3}  \f$
    !!
    !! where \f$ x_j \f$ is the \f$ j \f$th mixing term and \f$ j = 1,2,3 \f$.
    !!
    !! The molar excess Gibbs energy of mixing from a quaternary interaction term \f$ z \f$ is:
    !!
    !! \f$  g_{\lambda,z}^{ex} = x_1 x_2 x_3 x_4 {^0}L_{1,2,3,4}  \f$
    !!
    !! Finally, the total molar excess Gibbs energy of mixing for solution phase \f$ \lambda \f$ is
    !!
    !! \f$ g_{\lambda}^{ex} = \sum_{z=1}^Z g_{\lambda,z}^{ex} \f$
    !!
    !!
    !! For more information on the RKMP model, the Muggiano interpolation scheme and the molar excess Gibbs energy
    !! of mixing equations, the reader is referred to the following literature:
    !!
    !!      P. Chartrand and A. Pelton, "On the Choice of 'Geometric' Thermodynamic Models",
    !!      Journal of Phase Equilibria, 21, 2 (2000) 141-147.
    !!
!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iSolnIndex    Absolute index of a solution phase
    !
    ! nSpeciesPhase             An integec vector representing the last index of a species in a particular
    !                            solution phase
    ! nParamPhase               An integer vector representing the number of parameters for a particular phase.
    ! iParam                    An integer scalar representing the index number of a parameter.
    ! iRegularParam             An integer matrix representing the relative indices of species mixing for a
    !                            particular parmater.
    ! iFirstSpecies             An integer scalar representing the absolute index of the first species in a phase.
    ! iLastSpecies              An integer scalar representing the absolute index of the last species in a phase.
    ! dMolFraction              A double real vector representing hte mole fraction of all species in the system.
    ! dExcessGibbsParam         A double real vector representing the molar excess Gibbs energy of mixing for
    !                            each subsystem.
    ! dPartialExcessGibbs       Partial molar excess Gibbs energy of mixing of species.
    ! dMolFraction              Current estimated mole fraction.
    ! KD                        A double real scalar representing the Kronecker-Delta.  This is equal to unity
    !                            when a constituent in a ternary mixing term is equal to the higher order term.
    ! x1, x2, x3, x4            A double real scalar representing the mole fraction of the first, second, third
    !                            and fourth constituents in the parameter, respectively.
    ! dx                        A double real scalar representing the difference in x1 and x2.
    ! xprod                     A double real scalar representing the product of mole fractions of all
    !                            constituents associated with a parameter.
    ! dxvmo                     A double real scalar representing = (x1 - x2) ** (v-1).
    ! cSolnPhaseType            A character vector representing the solution phase type.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: i, j, iParam, iSolnIndex, iFirstSpecies, iLastSpecies, iExponent
    real(8) :: x1, x2, x3, x4, dx, xprod, xi, xj, KD, dxvmo, dx1,dx2,dx3
!
!
    ! Return control to the parent subroutine if there aren't any interaction parameters for this phase:
    if (nParamPhase(iSolnIndex) - nParamPhase(iSolnIndex-1) == 0) return
    
!
    ! Only proceed if the phase type is correct:
    IF_RKMP: if ((cSolnPhaseType(iSolnIndex) == 'RKMP').OR.(cSolnPhaseType(iSolnIndex) == 'RKMPM')) then
!
        ! Store the indices of the first and last species in this phase:
        iFirstSpecies = nSpeciesPhase(iSolnIndex-1) + 1
        iLastSpecies  = nSpeciesPhase(iSolnIndex)

        dPartialExcessGibbs(iFirstSpecies:iLastSpecies) = 0D0 !Added by SY
        dPartialEnthalpyXS(iFirstSpecies:iLastSpecies)  = 0D0 !Added by SY
        dPartialEntropyXS(iFirstSpecies:iLastSpecies)   = 0D0 !Added by SY
        dPartialHeatCapacityXS(iFirstSpecies:iLastSpecies) = 0D0 !Added by SY

!
        ! Loop through all interaction parameters in this phase:
        LOOP_Param: do iParam = nParamPhase(iSolnIndex-1)+1, nParamPhase(iSolnIndex)
            ! Compute temporary variables for sake of convenience:
            x1    = dMolFraction(iFirstSpecies + iRegularParam(iParam,2) - 1)
            x2    = dMolFraction(iFirstSpecies + iRegularParam(iParam,3) - 1)
            xprod = x1 * x2
            dx    = x1 - x2
!
            IF_ParamType: if (iRegularParam(iParam,1) == 2) then
!
                ! Binary parameter:
!
                ! Set to a small number when dx gets smaller
                if (dx == 0D0) dx=1D-10
!
                iExponent = iRegularParam(iParam,4)
                dxvmo     = dx**(iExponent-1)
                dxvmo     = DMIN1(dxvmo,1D30)
                dx1 = dxvmo * ((x2 - xprod)*dx + DFLOAT(iExponent)*xprod*(1D0 - dx))
                dx2 =dxvmo * ((x1 - xprod)*dx - DFLOAT(iExponent)*xprod*(1D0 + dx))
                dx3 =(dx**(iExponent)) * xprod * (1D0 + DFLOAT(iExponent))

!
                ! Loop through species in this phase:
                do i = iFirstSpecies, iLastSpecies
!
                    j = i - iFirstSpecies + 1
                    if (j == iRegularParam(iParam,2)) then
                        ! First species of parameter.
                        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) + dExcessGibbsParam(iParam)*dx1                        
                        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i) + dExcessHParam(iParam)*dx1
                        dPartialEntropyXS(i) = dPartialEntropyXS(i) + dExcessSParam(iParam)*dx1
                        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i) + dExcessCpParam(iParam)*dx1

                    elseif (j == iRegularParam(iParam,3)) then
                        ! Second species of parameter.
                        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) + dExcessGibbsParam(iParam) * dx2
                        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i) + dExcessHParam(iParam)* dx2
                        dPartialEntropyXS(i) = dPartialEntropyXS(i) + dExcessSParam(iParam)* dx2
                        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i) + dExcessCpParam(iParam)*dx2

                    else
                        ! This species does not belong to the parameter.
                        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) - dExcessGibbsParam(iParam) * dx3
                        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i) - dExcessHParam(iParam)* dx3
                        dPartialEntropyXS(i) = dPartialEntropyXS(i) - dExcessSParam(iParam) * dx3
                        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i) - dExcessCpParam(iParam)* dx3
                    end if
!
                end do
!
            elseif (iRegularParam(iParam,1) == 3) then
!
                ! Ternary parameter:
                x3    = dMolFraction(iFirstSpecies + iRegularParam(iParam,4) - 1)
                xj    = dMolFraction(iRegularParam(iParam,5) + iFirstSpecies - 1)
                xprod = xprod * x3
                
                dx2 = xprod * (x1 + x2 + x3 - 3D0*xj - 2D0/3D0)
!
                ! Loop through species in this phase:
                do i = iFirstSpecies, iLastSpecies
                    ! Relative species index in phase:
                    j = i - iFirstSpecies + 1!
                    ! Compute Kronecker-Delta term for ternary parameter:
                    KD = 0
                    if (iRegularParam(iParam,5) == j) KD = 1
!
                    
!
                    if ((j == iRegularParam(iParam,2)).OR.(j == iRegularParam(iParam,3)).OR. &
                        (j == iRegularParam(iParam,4))) then
!
                        ! This species contributes to the parameter.
                        xi = dMolFraction(i)
                        dx1 = xprod * (((1d0 / xi) - 3D0) * ((1D0 - x1 - x2 - x3)/3D0 + xj) + KD)
!
                        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) + dExcessGibbsParam(iParam) * dx1
                        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i) + dExcessHParam(iParam) * dx1
                        dPartialEntropyXS(i) = dPartialEntropyXS(i) + dExcessSParam(iParam) * dx1
                        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i) + dExcessCpParam(iParam) * dx1
!
                    else
                        ! This species does not contribute to the parameter.
                        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) + dExcessGibbsParam(iParam) * dx2
                        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i) + dExcessHParam(iParam) * dx2
                        dPartialEntropyXS(i) = dPartialEntropyXS(i) + dExcessSParam(iParam) * dx2
                        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i) + dExcessCpParam(iParam) * dx2
!
                    end if
                end do
!
            elseif (iRegularParam(iParam,1) == 4) then
!
                ! Quaternary parameter:
                x3    = dMolFraction(iFirstSpecies + iRegularParam(iParam,4) - 1)
                x4    = dMolFraction(iFirstSpecies + iRegularParam(iParam,5) - 1)
                xprod = xprod * x3 * x4
                
                dx2 = 3D0 * xprod 
!
                ! Loop through species in this phase:
                do i = iFirstSpecies, iLastSpecies
!
                    j = i - iFirstSpecies + 1
!
                    if ((j == iRegularParam(iParam,2)).OR.(j == iRegularParam(iParam,3)).OR. &
                        (j == iRegularParam(iParam,4)).OR.(j == iRegularParam(iParam,5))) then
!
                        ! This species contributes to the paramter:
                        xi = dMolFraction(i)
                        dx1 = xprod * ((1D0 / xi) - 3D0)

                        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) + dExcessGibbsParam(iParam) * dx1
                        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i) + dExcessHParam(iParam) * dx1
                        dPartialEntropyXS(i) = dPartialEntropyXS(i) + dExcessSParam(iParam) * dx1
                        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i) + dExcessCpParam(iParam) * dx1
                    else
                        ! This species does not contribute to the parameter:
                        dPartialExcessGibbs(i) = dPartialExcessGibbs(i) - dx2 * dExcessGibbsParam(iParam)
                        dPartialEnthalpyXS(i) = dPartialEnthalpyXS(i)  - dx2 * dExcessHParam(iParam)
                        dPartialEntropyXS(i) = dPartialEntropyXS(i)  - dx2 * dExcessSParam(iParam)
                        dPartialHeatCapacityXS(i) = dPartialHeatCapacityXS(i)  - dx2 * dExcessCpParam(iParam)
                    end if
!
                end do
            else
                ! The parameter index is not supported/recognized.  Report an error and exit.
                INFOThermo = 32
                exit LOOP_Param
!
            end if IF_ParamType
!
        end do LOOP_Param
!
    end if IF_RKMP
    
!
    return
!
end subroutine CompExcessGibbsEnergyRKMP
!
!
