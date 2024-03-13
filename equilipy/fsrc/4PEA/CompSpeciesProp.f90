



subroutine CompSpeciesProp(i,l,G,H,S,Cp)

    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO
    implicit none
    integer, intent(in)                :: i, l
    real(8), intent(out)               :: G, H, S, Cp
    integer                            :: k, iCounterGibbsEqn
    real(8)                            :: dLogT, dLogP
    real(8), dimension(6)              :: dGibbsCoeff, dEnthalpyCoeff, dEntropyCoeff, dCpCoeff

    ! Initialize variables:
    iCounterGibbsEqn = 0
    G  = 0D0
    H  = 0D0
    S  = 0D0
    Cp = 0D0
!
    dLogT            = DLOG(dTemperature)              ! ln(T)
    dLogP            = DLOG(dPressure)                 ! ln(P)

    ! Compute Gibbs energy coefficients:
    dGibbsCoeff(1)   = 1D0                             ! A
    dGibbsCoeff(2)   = dTemperature                    ! B
    dGibbsCoeff(3)   = dTemperature*dLogT ! C
    dGibbsCoeff(4)   = dTemperature**2                 ! D
    dGibbsCoeff(5)   = dTemperature**3                 ! E
    dGibbsCoeff(6)   = 1D0 / dTemperature              ! F

    ! S = -dG/dT
    dEntropyCoeff(1)   = 0D0                           ! A
    dEntropyCoeff(2)   = -1                            ! B
    dEntropyCoeff(3)   = -(dLogT + 1)     ! C
    dEntropyCoeff(4)   = -2*dTemperature               ! D
    dEntropyCoeff(5)   = -3*dTemperature**2            ! E
    dEntropyCoeff(6)   = 1D0 / dTemperature**2         ! F

    ! H = G + TS
    dEnthalpyCoeff(1)   = 1D0                          ! A
    dEnthalpyCoeff(2)   = 0D0                          ! B
    dEnthalpyCoeff(3)   = -dTemperature                ! C
    dEnthalpyCoeff(4)   = -dTemperature**2             ! D
    dEnthalpyCoeff(5)   = -2*dTemperature**3           ! E
    dEnthalpyCoeff(6)   = 2D0 / dTemperature           ! F

    ! Cp = dH/dT = TdS/dT
    dCpCoeff(1)   = 0D0                                ! A
    dCpCoeff(2)   = 0D0                                ! B
    dCpCoeff(3)   = -1                                 ! C
    dCpCoeff(4)   = -2*dTemperature                    ! D
    dCpCoeff(5)   = -6*dTemperature**2                 ! E
    dCpCoeff(6)   = -2D0 / dTemperature**2             ! F


    do k = 2, 7
        G  = G  + dGibbsCoeffSpeciesTemp(k,l) * dGibbsCoeff(k-1)
        S  = S  + dGibbsCoeffSpeciesTemp(k,l) * dEntropyCoeff(k-1)
        H  = H  + dGibbsCoeffSpeciesTemp(k,l) * dEnthalpyCoeff(k-1)
        Cp = Cp + dGibbsCoeffSpeciesTemp(k,l) * dCpCoeff(k-1)
    end do

    ! Compute additional standard molar Gibbs energy terms:
    if (dGibbsCoeffSpeciesTemp(9,l) .EQ. 99) then
        G  = G  + dGibbsCoeffSpeciesTemp(8,l) * dLogT
        S  = S  - dGibbsCoeffSpeciesTemp(8,l) / dTemperature
        H  = H  + dGibbsCoeffSpeciesTemp(8,l) * (dLogT-1D0)
        Cp = Cp + dGibbsCoeffSpeciesTemp(8,l) / dTemperature
    else
        G  = G  + dGibbsCoeffSpeciesTemp(8,l) * dTemperature**dGibbsCoeffSpeciesTemp(9,l)
        if(dGibbsCoeffSpeciesTemp(9,l)==0D0) then
            S  = S
            H  = H  + dGibbsCoeffSpeciesTemp(8,l)
            Cp = Cp
        else
            S  = S  + dGibbsCoeffSpeciesTemp(8,l) * (-dGibbsCoeffSpeciesTemp(9,l))*dTemperature**(dGibbsCoeffSpeciesTemp(9,l)-1)
            H  = H  + dGibbsCoeffSpeciesTemp(8,l) * (1-dGibbsCoeffSpeciesTemp(9,l))*dTemperature**dGibbsCoeffSpeciesTemp(9,l)
            Cp = Cp + dGibbsCoeffSpeciesTemp(8,l) * &
            (1-dGibbsCoeffSpeciesTemp(9,l))*dGibbsCoeffSpeciesTemp(9,l)*dTemperature**(dGibbsCoeffSpeciesTemp(9,l)-1)
        end if
    end if

    if (dGibbsCoeffSpeciesTemp(11,l) .EQ. 99) then
        G  = G  + dGibbsCoeffSpeciesTemp(10,l) * dLogT
        S  = S  - dGibbsCoeffSpeciesTemp(10,l) / dTemperature
        H  = H  + dGibbsCoeffSpeciesTemp(10,l) * (dLogT-1D0)
        Cp = Cp + dGibbsCoeffSpeciesTemp(10,l) / dTemperature
    else
        G  = G  + dGibbsCoeffSpeciesTemp(10,l) * dTemperature**dGibbsCoeffSpeciesTemp(11,l)
        if(dGibbsCoeffSpeciesTemp(11,l)==0D0) then
            S  = S
            H  = H  + dGibbsCoeffSpeciesTemp(10,l)
            Cp = Cp 
        else
            S  = S  + dGibbsCoeffSpeciesTemp(10,l) * (-dGibbsCoeffSpeciesTemp(11,l))*dTemperature**(dGibbsCoeffSpeciesTemp(11,l)-1)
            H  = H  + dGibbsCoeffSpeciesTemp(10,l) * (1-dGibbsCoeffSpeciesTemp(11,l))*dTemperature**dGibbsCoeffSpeciesTemp(11,l)
            Cp = Cp + dGibbsCoeffSpeciesTemp(10,l) * &
            (1-dGibbsCoeffSpeciesTemp(11,l))*dGibbsCoeffSpeciesTemp(11,l)*dTemperature**(dGibbsCoeffSpeciesTemp(11,l)-1)
        end if
    end if

    if (dGibbsCoeffSpeciesTemp(13,l) .EQ. 99) then
        G  = G  + dGibbsCoeffSpeciesTemp(12,l) * dLogT
        S  = S  - dGibbsCoeffSpeciesTemp(12,l) / dTemperature
        H  = H  + dGibbsCoeffSpeciesTemp(12,l) * (dLogT-1D0)
        Cp = Cp + dGibbsCoeffSpeciesTemp(12,l) / dTemperature
    else
        G  = G  + dGibbsCoeffSpeciesTemp(12,l) * dTemperature**dGibbsCoeffSpeciesTemp(13,l)
        if(dGibbsCoeffSpeciesTemp(13,l)==0D0) then
            S  = S
            H  = H  + dGibbsCoeffSpeciesTemp(12,l)
            Cp = Cp 
        else
            S  = S  + dGibbsCoeffSpeciesTemp(12,l) * (-dGibbsCoeffSpeciesTemp(13,l))*dTemperature**(dGibbsCoeffSpeciesTemp(13,l)-1)
            H  = H  + dGibbsCoeffSpeciesTemp(12,l) * (1-dGibbsCoeffSpeciesTemp(13,l))*dTemperature**dGibbsCoeffSpeciesTemp(13,l)
            Cp = Cp + dGibbsCoeffSpeciesTemp(12,l) * &
            (1-dGibbsCoeffSpeciesTemp(13,l))*dGibbsCoeffSpeciesTemp(13,l)*dTemperature**(dGibbsCoeffSpeciesTemp(13,l)-1)
        end if
    end if

end subroutine CompSpeciesProp