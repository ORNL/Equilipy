!
!
subroutine CompGibbsMagneticSolnInit(i,j)
!
    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO
!
    implicit none
!
    integer :: i, j
    real(8) :: B, D, p, invpmone, tau, Tcritical, f, fp, fdp, StructureFactor
    real(8) :: dTempA, dTempB, dTempC
!
!
    ! Assign the critical temperature can either be the Curie temperature for ferromagnetic materials
    ! or the Neel temperature for antiferromagnetic materials:
    Tcritical       = dCoeffGibbsMagnetic(j,1)
    B               = dCoeffGibbsMagnetic(j,2)
    StructureFactor = dCoeffGibbsMagnetic(j,3)
    p               = dCoeffGibbsMagnetic(j,4)
    invpmone        = 1D0/p - 1D0
    dTempA          = 0D0
    dTempB          = 0D0
    dTempC          = 0D0
!
    ! ChemSage files store the critical (i.e., Neel) temperature for antiferromagnetic materials
    ! as a negative real value divided by the structure factor.  Correct Tcritical and B:
    if (Tcritical < 0D0) then
        Tcritical = -Tcritical * StructureFactor
        B         = -B * StructureFactor
    end if
!
    tau = dTemperature / Tcritical
    D   = 518D0/1125D0 + (11692D0/15975D0) * invpmone
!
    if (tau > 1D0) then
        dTempA = tau**(-5)                  ! tau^(-5)
        dTempB = dTempA**(3)                ! tau^(-15)
        dTempC = dTempA * dTempA * dTempB   ! tau^(-25)
        f   = -(dTempA/10D0 + dTempB/315D0 + dTempC/1500D0) / D
        fp  = (dTempA/2D0 + dTempB/21D0 + dTempC/60D0) / (D*tau)
        fdp = -(3*dTempA + 16*dTempB/21D0 + 13*dTempC/30D0) / (D*tau**2)
    else
        dTempA = tau**(3)                   ! tau^(3)
        dTempB = dTempA**(3)                ! tau^(9)
        dTempC = dTempA * dTempA * dTempB   ! tau^(15)
        f   = 1D0 - (79D0/(140D0*p*tau) + (474D0/497D0)*invpmone*(dTempA/6D0 + dTempB/135D0 + dTempC/600D0)) / D
        fp  = (79D0/(140D0*p*tau) - (474D0/497D0)*invpmone*(dTempA/2D0 + dTempB/15D0 + dTempC/40D0)) / (D*tau)
        fdp = -(79D0/(70D0*p*tau) + (474D0/497D0)*invpmone*(dTempA + 8*dTempB/15D0 + 7*dTempC/20D0)) / (D*tau**2)
    end if
!
    ! Add the magnetic contribution to the Gibbs energy term (unitless)
    dMagGibbsEnergy(j)  = DLOG(B + 1D0) * f
    dMagEnthalpy(j)     = -DLOG(B + 1D0) * tau*fp
    dMagEntropy(j)      = -DLOG(B + 1D0) * (f + tau*fp)
    dMagHeatCapacity(j) = -DLOG(B + 1D0) * tau * (2*fp + tau*fdp)
!
    return
!
end subroutine CompGibbsMagneticSolnInit
!
!
