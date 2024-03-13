!
!
subroutine CompGibbsMagneticSoln(iSolnPhaseIndex)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompGibbsMagneticSoln.f90
    !> \brief   Compute magnetic contributions to the Gibbs energy terms for a solution phase.
    !> \author  M.H.A. Piro
    !> \date    March 7, 2013
    !> \sa      CompExcessGibbsEnergy.f90
    !
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    03/07/2013    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the magnetic contribution to the standard molar
    !! Gibbs energy of a solution species.  This contribution is given by
    !! \f$ \Delta g_{mag} = RT ln(B_o + 1) g(\tau ) \f$, where \f$ B_o \f$ is the average magnetic moment per
    !! atom, \f$ \tau \f$  is the absolute temperature divided by the critical temperature (i.e., the Curie
    !! temperature for ferromagnetic materials or the Neel temperature for antiferromagnetic materials) and
    !! \f$ g \f$ is a function of \f$ \tau \f$, given by:
    !! \f$ g(\tau ) = 1 - \left( \frac{79\tau^-1}{140p} + \frac{474}{497}(\frac{1}{p} - 1)(\frac{\tau ^3}{6}
    !!   + \frac{\tau ^9}{135} \frac{\tau ^15}{600}   )  \right) /D, \tau \leq 1 \f$
    !! and \f$ g(\tau ) = - \left( \frac{\tau ^{-5}}{10} + \frac{\tau ^{-15}}{315} + \frac{\tau ^{-25}}{1500}
    !!   \right) /D, \tau > 1 \f$, where \f$ D = \frac{518}{1125} + \frac{11692}{15975} \left( \frac{1}{p} -1
    !!   \right) \f$.
    !!
    !
    ! References:
    ! ===========
    !
    !> \details The following references explain the magnetic contribution to the Gibbs energy term that is
    !! used in this subroutine:
    !!
    !!   M. Hillert and M. Jarl, "A Model for Alloying Effects in Ferromagnetic Alloys," CALPHAD, 2, 3
    !!   (1978) 227-238.
    !!
    !!   A.T. Dinsdale, "SGTE Data for Pure Elements," CALPHAD, 15, 4 (1991) 317-425.
    !!
    !!   H.L. Lukas, S.G. Fries and B. Sundman, "Computational Thermodynamics: The Calphad Method," Cambridge
    !!   University Press, New York (2007).
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iSolnPhaseIndex   An integer scalar representing the absolute solution phase index.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
!
    implicit none
!
    integer :: i, j, iFirst, iLast, iSolnPhaseIndex, iExponent, iParam
    integer :: k, m, n, s, c, sd, KD, l
    integer :: nSublattice, iChargedPhaseID
    integer :: iFirstParam, iSecondParam, iSubParam
    real(8) :: B, D, p, invpmone, tau, Tcritical, g, StructureFactor
    real(8) :: dTemp, dTempA, dTempB, dTempC, dTempD, dTempCoeff1, dTempCoeff2
    real(8) :: x1, x2, dx, xprod, dxvmo
    real(8) :: f,fp, fdp, ftp, qp, qdp, qtp
    real(8) :: dTempG, dTempH, dTempS, dTempCp
    real(8) :: dPreFactor, dPreFactorG, dPreFactorH, dPreFactorS, dPreFactorCp
    real(8) :: Param1, Param2, Param3
    logical :: lTn, lBn
!
!
    
!
    ! Store the first and last species indices:
    iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast  = nSpeciesPhase(iSolnPhaseIndex)

    ! Exit when only one endmember appears stable
    if(iLast-iFirst<=1) return

    ! Initialize variables:
    Tcritical = 0D0
    B         = 0D0
    dMagGibbsEnergy(iFirst:iLast)     = 0D0
    dMagEnthalpy(iFirst:iLast)        = 0D0
    dMagEntropy(iFirst:iLast)         = 0D0
    dMagHeatCapacity(iFirst:iLast)    = 0D0
    dTempG= 0D0
    dTempH= 0D0
    dTempS= 0D0
    dTempCp = 0D0
    

    ! Note: the structure factor and p factor are the same for all constituents in the phase.
    StructureFactor = dCoeffGibbsMagnetic(iFirst,3)
    p               = dCoeffGibbsMagnetic(iFirst,4)
    invpmone        = 1D0/p - 1D0
!
    call CompMagneticTemperatureMoment(iSolnPhaseIndex,Tcritical,B)
!
    ! ChemSage files store the critical temperature for antiferromagnetic materials
    ! (i.e., the Neel temperature) as a negative real value divided by the structure factor.
    lBn = .FALSE.
    lTn = .FALSE.
    if (B < 0D0) then
        lBn = .TRUE.
        B         = -B * StructureFactor
    end if
    if (Tcritical < 0D0) then
        lTn = .TRUE.
        Tcritical = -Tcritical * StructureFactor
    end if
!
    ! Only proceed if the critical temperature is not zero:
    IF_Proceed: if (Tcritical /= 0D0) then
!
        ! Compute model parameters:
        tau = dTemperature / Tcritical
        D   = (518D0/1125D0) + (11692D0/15975D0) * invpmone
!
        ! The magnetic model of Hillert and Jarl is empirical and depends on tau:
        IF_Tau: if (tau > 1D0) then
            dTempA = tau**(-5)                  ! tau^(-5)
            dTempB = dTempA**(3)                ! tau^(-15)
            dTempC = dTempA * dTempA * dTempB   ! tau^(-25)
            
            f   = -(dTempA/10D0 + dTempB/315D0 + dTempC/1500D0) / D
            fp  =  (1D0 / (D * Tcritical))*(dTempA/2D0 + dTempB/21D0 + dTempC/60D0)  ! tau * f'/Tc
            fdp = -(1D0 / (D * Tcritical))*(3D0*dTempA + 16D0*dTempB/21D0 + 13D0*dTempC/30D0)  ! tau**2 * f''/Tc
            ftp = +(1D0 / (D * Tcritical))*(21D0*dTempA + 272*dTempB/21D0 + 117D0*dTempC/10D0) ! tau**3 * f'''/Tc

            qp = fp
            qdp = fdp
            qtp = ftp
        else
            dTempA = tau**(3)                   ! tau^(3)
            dTempB = dTempA**(3)                ! tau^(9)
            dTempC = dTempA * dTempA * dTempB   ! tau^(15)

            f   = 1D0 - (79D0/(140D0*p*tau) + (474D0/497D0)*invpmone*(dTempA/6D0 + dTempB/135D0 + dTempC/600D0)) / D
            fp = (79D0/(140D0*p*dTemperature)&
                - (474D0/(497D0*Tcritical))*invpmone*(dTempA/2D0 + dTempB/15D0 + dTempC/40D0)) / D     ! tau * f'/Tc
            fdp = -(79D0/(70D0*p*dTemperature)&
                + (474D0/(497D0*Tcritical))*invpmone*(dTempA + 8D0*dTempB/15D0 + 7D0*dTempC/20D0)) / D  ! tau**2 * f''/Tc
            ftp = (237D0/(70D0*p*dTemperature)&
                +ftp- (474D0/(497D0*Tcritical))*invpmone*(dTempA + 56D0*dTempB/15D0 + 91D0*dTempC/20D0))/D ! tau**3 * f'''/Tc

            qp = -fp
            qdp = -fdp
            qtp = -ftp
        end if IF_Tau
!
        
        ! Loop through species in this phase and update the chemical potential:
        do i = iFirst, iLast
            dTempCoeff1 = dCoeffGibbsMagnetic(i,1)
            dTempCoeff2 = dCoeffGibbsMagnetic(i,2)
            if (lTn) then
                dTempCoeff1 = -dTempCoeff1 * StructureFactor
            end if
            if (lBn) then
                dTempCoeff2 = -dTempCoeff2 * StructureFactor
            end if
            Param1= (dTempCoeff2 - B) / (1D0 + B)
            Param2 = DLOG(1D0 + B)
            Param3 = -(Tcritical - dTempCoeff1)

            dTempG = f * Param1 + Param2 * (Param3* qp+ f)
            dTempH = -(Tcritical*Param1*fp+Param2*(Tcritical*fp+Param3*(qp+qdp)))
            dTempS = -Param1*(f+Tcritical*fp) - Param2*(f+Tcritical*fp+Param3*(2*qp+qdp))
            dTempCp = -Param1*Tcritical*(2*fp+fdp) - Param2*(Tcritical*(2*fp+fdp)+Param3*(2*qp+4*qdp+qtp))


            if (cSolnPhaseType(iSolnPhaseIndex) == 'RKMPM') then
                dMagGibbsEnergy(i) = dTempG
                dMagEnthalpy(i) = dTempH
                dMagEntropy(i) = dTempS
                dMagHeatCapacity(i) = dTempCp
            else if (cSolnPhaseType(iSolnPhaseIndex) == 'SUBLM') then
                m = i - iFirst + 1
                iChargedPhaseID = iPhaseSublattice(iSolnPhaseIndex)
                nSublattice     = nSublatticePhase(iChargedPhaseID)
                do j = iFirst, iLast
                    ! Relative species index:
                    n = j - iFirst + 1
!
                    ! Compute pre-factor term:
                    dPreFactor = 1D0 - DFLOAT(nSublattice)
!
                    ! Loop through sublattices:
                    do s = 1, nSublattice
                        ! Store constituent indices:
                        k = iConstituentSublattice(iChargedPhaseID,s,m)
                        l = iConstituentSublattice(iChargedPhaseID,s,n)
!
                        ! Effectively apply Kronecker-Delta term to pre-factor:
                        if (k == l)  dPreFactor = dPreFactor + 1D0 / dSiteFraction(iChargedPhaseID,s,k)
                    end do
!
                    ! Update the reference molar Gibbs energy:
                    dMagGibbsEnergy(j) = dMagGibbsEnergy(j) + dPreFactor * dMolFraction(i) * dTempG
                    dMagEnthalpy(j) = dMagEnthalpy(j) + dPreFactor * dMolFraction(i) * dTempH
                    dMagEntropy(j) = dMagEntropy(j) + dPreFactor * dMolFraction(i) * dTempS
                    dMagHeatCapacity(j) = dMagHeatCapacity(j) + dPreFactor * dMolFraction(i) * dTempCp
                end do
!
            end if
        end do
!
        LOOP_Param: do iParam = nMagParamPhase(iSolnPhaseIndex-1)+1, nMagParamPhase(iSolnPhaseIndex)
            dTempCoeff1 = dMagneticParam(iParam,1)
            dTempCoeff2 = dMagneticParam(iParam,2)
            if (lTn) then
                dTempCoeff1 = -dTempCoeff1 * StructureFactor
            end if
            if (lBn) then
                dTempCoeff2 = -dTempCoeff2 * StructureFactor
            end if
            Param1= dTempCoeff2 / (1D0 + B)
            Param2 = DLOG(1D0 + B)
            Param3 = dTempCoeff1
            dTempG = f * Param1 + Param2 * (Param3* qp)

            dTempH = -(Tcritical*Param1*fp+Param2*Param3*(qp+qdp))
            dTempS = -Param1*(f+Tcritical*fp) - Param2*Param3*(2*qp+qdp)
            dTempCp = -Param1*Tcritical*(2*fp+fdp) 
            dTempCp = dTempCp - Param2*Param3*(2*qp+4*qdp+qtp)

            if (cSolnPhaseType(iSolnPhaseIndex) == 'RKMPM') then
                ! Compute temporary variables for sake of convenience:
                if (iMagneticParam(iParam,1) == 2) then
                    ! Binary parameter:
                    x1    = dMolFraction(iFirst + iMagneticParam(iParam,2) - 1)
                    x2    = dMolFraction(iFirst + iMagneticParam(iParam,3) - 1)
                    iExponent = iMagneticParam(iParam,4)
                    xprod = x1 * x2
                    dx    = x1 - x2
                    dxvmo     = dx**(iExponent-1)
                    dxvmo     = DMIN1(dxvmo,1D30)
                    ! Cycle if dx = 0 to prevent calculating either an INF or a NAN:
                    if (dx == 0D0) cycle LOOP_Param
                    do i = iFirst, iLast
                        j = i - iFirst + 1
                        if (j == iMagneticParam(iParam,2)) then
                            ! First species of parameter.
                            dMagGibbsEnergy(i) = dMagGibbsEnergy(i) + dTempG * &
                                dxvmo * ((x2 - xprod)*dx + DFLOAT(iExponent)*xprod*(1D0 - dx))

                            dMagEnthalpy(i) = dMagEnthalpy(i) + dTempH * &
                                dxvmo * ((x2 - xprod)*dx + DFLOAT(iExponent)*xprod*(1D0 - dx))
                            dMagEntropy(i) = dMagEntropy(i) + dTempS * &
                                dxvmo * ((x2 - xprod)*dx + DFLOAT(iExponent)*xprod*(1D0 - dx))
                            dMagHeatCapacity(i) = dMagHeatCapacity(i) + dTempCp * &
                                dxvmo * ((x2 - xprod)*dx + DFLOAT(iExponent)*xprod*(1D0 - dx))
!
                        elseif (j == iMagneticParam(iParam,3)) then
                            ! Second species of parameter.
                            dMagGibbsEnergy(i) = dMagGibbsEnergy(i) + dTempG * &
                                dxvmo * ((x1 - xprod)*dx - DFLOAT(iExponent)*xprod*(1D0 + dx))

                            dMagEnthalpy(i) = dMagEnthalpy(i) + dTempH * &
                                dxvmo * ((x1 - xprod)*dx - DFLOAT(iExponent)*xprod*(1D0 + dx))
                            dMagEntropy(i) = dMagEntropy(i) + dTempS * &
                                dxvmo * ((x1 - xprod)*dx - DFLOAT(iExponent)*xprod*(1D0 + dx))
                            dMagHeatCapacity(i) = dMagHeatCapacity(i) + dTempCp * &
                                dxvmo * ((x1 - xprod)*dx - DFLOAT(iExponent)*xprod*(1D0 + dx))
                        else
                            ! This species does not belong to the parameter.
                            dMagGibbsEnergy(i) = dMagGibbsEnergy(i) - dTempG * &
                                (dx**(iExponent)) * xprod * (1D0 + DFLOAT(iExponent))
                            
                            dMagEnthalpy(i) = dMagEnthalpy(i) - dTempH * &
                                (dx**(iExponent)) * xprod * (1D0 + DFLOAT(iExponent))
                            dMagEntropy(i) = dMagEntropy(i) - dTempS * &
                                (dx**(iExponent)) * xprod * (1D0 + DFLOAT(iExponent))
                            dMagHeatCapacity(i) = dMagHeatCapacity(i) - dTempCp * &
                                (dx**(iExponent)) * xprod * (1D0 + DFLOAT(iExponent))
                        end if
!
                    end do
                else
                    ! The parameter index is not supported/recognized.  Report an error and exit.
                    INFOThermo = 43
                    exit LOOP_Param
                end if
            else if (cSolnPhaseType(iSolnPhaseIndex) == 'SUBLM') then
                iChargedPhaseID = iPhaseSublattice(iSolnPhaseIndex)
                nSublattice     = nSublatticePhase(iChargedPhaseID)
                ! Reinitialize temporary variable:
                dPreFactor = 1D0
!
                ! Store the number of constituents involved in this parameter:
                n = iMagneticParam(iParam,1)
                iExponent = iMagneticParam(iParam,n+2)
!
                ! Loop through constituents associated with this parameter:
                do k = 2, n + 1
!
                    ! Determine constituent and sublattice indices:
                    c = MOD(iMagneticParam(iParam,k), 10000)
                    s = iMagneticParam(iParam,k) - c
                    s = s / 10000
!
                    ! Compute prefactor term:
                    dPreFactor = dPreFactor * dSiteFraction(iChargedPhaseID,s,c)
!
                    ! Store the first and second site fractions:
                    ! This assumes that the constituents that are mixing are the first two listed:
                    if (k == 2) then
                        x1 = dSiteFraction(iChargedPhaseID,s,c)
                        iFirstParam = c
                        iSubParam   = s
                    elseif (k == 3) then
                        x2 = dSiteFraction(iChargedPhaseID,s,c)
                        iSecondParam = c
                    end if
                end do
!
                ! Multiply prefactor term by excess Gibbs energy parameter:
                dPreFactorG = dPreFactor *dTempG * (x1 - x2)**iExponent
                dPreFactorH = dPreFactor *dTempH * (x1 - x2)**iExponent
                dPreFactorS = dPreFactor *dTempS * (x1 - x2)**iExponent
                dPreFactorCp = dPreFactor *dTempCp * (x1 - x2)**iExponent
!
                LOOP_Param_Species: do i = iFirst, iLast
                    ! Reinitialize variables:
                    KD    = 0
                    m     = i - iFirst + 1
                    dTemp = -DFLOAT(nSublattice + iMagneticParam(iParam,n+2))
!
                    ! Loop through sublattices associated with this phase:
                    LOOP_Param_Sub: do s = 1, nSublattice
                        ! Store constituent index corresponding to component i on sublattice s:
                        c = iConstituentSublattice(iChargedPhaseID,s,m)
!
                        ! Assign Kronecker-Delta term if this component contains the constituent corresponding
                        ! to the mixing parameter:
                        if (s == iSubParam) then
                            if (c == iFirstParam) then
                                ! This is the first mixing constituent:
                                KD = 1
                            elseif (c == iSecondParam) then
                                ! This is the second mixing constituent:
                                KD = -1
                            else
                                ! The constituents don't match:
                                KD = 0
                            end if
                        end if
!
                        ! Loop through constituents involved in this mixing parameter:
                        LOOP_Param_Const: do j = 2, n + 1
                            k = MOD(iMagneticParam(iParam,j), 10000)  ! Constituent index corresponding to parameter.
                            sd = iMagneticParam(iParam,j) - k
                            sd = sd / 10000                       ! Sublattice index corresponding to parameter.
!
                            ! Cycle if they are on different sublattices:
                            if (sd /= s) cycle LOOP_Param_Const
!
                            ! Cycle if they are different constituents:
                            if (k /= c) cycle LOOP_Param_Const
!
                            ! Include contribution from this site fraction:
                            dTemp = dTemp + 1D0 / dSiteFraction(iChargedPhaseID,s,c)
!
                        end do LOOP_Param_Const ! j
                    end do LOOP_Param_Sub       ! s
!
                    ! Apply higher order terms (only if x1 and x2 are not the same):
                    if (x1 /= x2) then
                        dTemp = dTemp + DFLOAT(KD * iMagneticParam(iParam,n+2)) / (x1 - x2)
                    end if
!
                    ! Apply partial molar excess Gibbs energy of mixing:
                    dMagGibbsEnergy(i) = dMagGibbsEnergy(i) + dPreFactorG * dTemp 
                    dMagEnthalpy(i) = dMagEnthalpy(i) + dPreFactorH * dTemp
                    dMagEntropy(i) = dMagEntropy(i) + dPreFactorS * dTemp 
                    dMagHeatCapacity(i) = dMagHeatCapacity(i) + dPreFactorCp * dTemp 
!
                end do LOOP_Param_Species       ! i
            end if
        end do LOOP_Param
!
    end if IF_Proceed
!
    return
!
end subroutine CompGibbsMagneticSoln
!
!
