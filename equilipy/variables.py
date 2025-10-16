

# This function is based on list-directed IO

# Global variables for ParseCS73
# Variables to pass on to equilifort (required for minimization)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# 1. Integers
global nElementsCS, nSpeciesCS,nSolnPhasesSysCS,nCountSublatticeCS,\
    nMaxSpeciesPhaseCS,nParamCS,nMagParamCS

# parameters: Make sure that these parametersa are consistnent in ModuleParseCS
global nMaxSublatticeCS,nSolnPhasesSysMax, nGibbsCoeff, nMaxGibbsEqs, nParamMax
nMaxSublatticeCS  = int(5)
nSolnPhasesSysMax = int(500)
nGibbsCoeff       = int(13)
nMaxGibbsEqs      = int(6)
nParamMax         = int(4)

# 2 Arrays (1D tensor)
global iPhaseCS,iParticlesPerMoleCS,iPhaseSublatticeCS,nGibbsEqSpecies,nSublatticePhaseCS,\
    nParamPhaseCS,nSpeciesPhaseCS,nMagParamPhaseCS,dAtomicMass,cRegularParamCS,cElementNameCS,\
    cSolnPhaseTypeCS,cSolnPhaseNameCS,cSpeciesNameCS

# 3 Matrix (2D tensor)
global nPairsSROCS,nConstituentSublatticeCS, nSublatticeElementsCS, dGibbsMagneticCS,\
    dStoichSublatticeCS,dZetaSpeciesCS,dGibbsCoeffSpeciesTemp,dStoichSpeciesCS,\
    iMagneticParamCS,dMagneticParamCS,iRegularParamCS,dRegularParamCS
global cPairNameCS

# 4. 3D tensor 
global iConstituentSublatticeCS, iPairIDCS, iChemicalGroupCS, dSublatticeChargeCS,\
    dStoichPairsCS, dConstituentCoefficientsCS,dCoordinationNumberCS, cConstituentNameSUBCS
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Variables for internal use>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

# 1 Integers
global nPureSpeciesCS,nSolnPhaseCS,indx,INFO,DataBase

# 2 Arrays (1D tensor)
global iCounterGibbsEqn, iMagParamPassCS,iParamPassCS, iElementSys, iElementSysIndex,\
    iSys2DBSoln, iSys2DBSpecies, iSys2DBComp,iPhaseSys,iSolnPS, iSys2DBOld, \
    iSys2DBSolnDefault, iSys2DBSpeciesDefault, iSys2DBCompDefault, iElementDBIndex
# Questions
# what do iparampasscs and imagparampasscs do? These are defined by nParamCS,nMagParamCS

# 3 Matrix (2D tensor) 
global dTempVec, dStoichConstituentCS, dMolarVolume, dConditionSys

# 4. String/Charactors for internal use
global cPhaseNames,cSolnPhaseTypeSupport,PhaseNameSys,cEndmemberNameSys,\
    cEndmemberNameCS, cPhaseNameSys,cSpeciesNameSys, cComponentNameSys, cPeriodicTable,\
    cPhaseNameSysDefault

global FunctionsHSCp, FactSage8Plus, CompoundGibbs

global dSpeciesFractionOut, dPhaseAmountOut, dTOut, dPOut
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

cPeriodicTable=dict({
    'e': [0, 0.0],
    'H': [1, 1.008],
    'He': [2, 4.002602],
    'Li': [3, 6.94],
    'Be': [4, 9.0121831],
    'B': [5, 10.81],
    'C': [6, 12.011],
    'N': [7, 14.007],
    'O': [8, 15.999],
    'F': [9, 18.99840316],
    'Ne': [10, 20.1797],
    'Na': [11, 22.98976928],
    'Mg': [12, 24.305],
    'Al': [13, 26.9815384],
    'Si': [14, 28.085],
    'P': [15, 30.973762],
    'S': [16, 32.06],
    'Cl': [17, 35.45],
    'Ar': [18, 39.948],
    'K': [19, 39.0983],
    'Ca': [20, 40.078],
    'Sc': [21, 44.955908],
    'Ti': [22, 47.867],
    'V': [23, 50.9415],
    'Cr': [24, 51.9961],
    'Mn': [25, 54.938043],
    'Fe': [26, 55.845],
    'Co': [27, 58.933194],
    'Ni': [28, 58.6934],
    'Cu': [29, 63.546],
    'Zn': [30, 65.38],
    'Ga': [31, 69.723],
    'Ge': [32, 72.63],
    'As': [33, 74.921595],
    'Se': [34, 78.971],
    'Br': [35, 79.904],
    'Kr': [36, 83.798],
    'Rb': [37, 85.4678],
    'Sr': [38, 87.62],
    'Y': [39, 88.90584],
    'Zr': [40, 91.224],
    'Nb': [41, 92.90637],
    'Mo': [42, 95.95],
    'Tc': [43, 97.0],
    'Ru': [44, 101.07],
    'Rh': [45, 102.90549],
    'Pd': [46, 106.42],
    'Ag': [47, 107.8682],
    'Cd': [48, 112.414],
    'In': [49, 114.818],
    'Sn': [50, 118.71],
    'Sb': [51, 121.76],
    'Te': [52, 127.6],
    'I': [53, 126.90447],
    'Xe': [54, 131.293],
    'Cs': [55, 132.905452],
    'Ba': [56, 137.327],
    'La': [57, 138.90547],
    'Ce': [58, 140.116],
    'Pr': [59, 140.90766],
    'Nd': [60, 144.242],
    'Pm': [61, 145.0],
    'Sm': [62, 150.36],
    'Eu': [63, 151.964],
    'Gd': [64, 157.25],
    'Tb': [65, 158.925354],
    'Dy': [66, 162.5],
    'Ho': [67, 164.930328],
    'Er': [68, 167.259],
    'Tm': [69, 168.934218],
    'Yb': [70, 173.045],
    'Lu': [71, 174.9668],
    'Hf': [72, 178.486],
    'Ta': [73, 180.94788],
    'W': [74, 183.84],
    'Re': [75, 186.207],
    'Os': [76, 190.23],
    'Ir': [77, 192.217],
    'Pt': [78, 195.084],
    'Au': [79, 196.96657],
    'Hg': [80, 200.592],
    'Tl': [81, 204.38],
    'Pb': [82, 207.2],
    'Bi': [83, 208.9804],
    'Po': [84, 209.0],
    'At': [85, 210.0],
    'Rn': [86, 222.0],
    'Fr': [87, 223.0],
    'Ra': [88, 226.0],
    'Ac': [89, 227.0],
    'Th': [90, 232.0377],
    'Pa': [91, 231.03588],
    'U': [92, 238.02891],
    'Np': [93, 237.0],
    'Pu': [94, 244.0],
    'Am': [95, 243.0],
    'Cm': [96, 247.0],
    'Bk': [97, 247.0],
    'Cf': [98, 251.0],
    'Es': [99, 252.0],
    'Fm': [100, 257.0],
    'Md': [101, 258.0],
    'No': [102, 259.0],
    'Lr': [103, 262.0],
    'Rf': [104, 267.0],
    'Db': [105, 270.0],
    'Sg': [106, 269.0],
    'Bh': [107, 270.0],
    'Hs': [108, 270.0],
    'Mt': [109, 278.0],
    'Ds': [110, 281.0],
    'Rg': [111, 281.0],
    'Cn': [112, 285.0],
    'Nh': [113, 286.0],
    'Fl': [114, 289.0],
    'Mc': [115, 289.0],
    'Lv': [116, 293.0],
    'Ts': [117, 293.0],
    'Og': [118, 294.0]
})

def to_dict():
    db = dict({})
    db['nElementsCS'] = nElementsCS
    db['nSpeciesCS'] = nSpeciesCS
    db['nSolnPhasesSysCS'] = nSolnPhasesSysCS
    db['nCountSublatticeCS'] = nCountSublatticeCS
    db['nMaxSpeciesPhaseCS'] = nMaxSpeciesPhaseCS
    db['nParamCS'] = nParamCS
    db['nMagParamCS'] = nMagParamCS
    db['nMaxSublatticeCS'] = nMaxSublatticeCS
    db['nSolnPhasesSysMax'] = nSolnPhasesSysMax
    db['nGibbsCoeff'] = nGibbsCoeff
    db['nMaxGibbsEqs'] = nMaxGibbsEqs
    db['nParamMax'] = nParamMax
    db['iPhaseCS'] = iPhaseCS
    db['iParticlesPerMoleCS'] = iParticlesPerMoleCS
    db['iPhaseSublatticeCS'] = iPhaseSublatticeCS
    db['nGibbsEqSpecies'] = nGibbsEqSpecies
    db['nSublatticePhaseCS'] = nSublatticePhaseCS
    db['nParamPhaseCS'] = nParamPhaseCS
    db['nSpeciesPhaseCS'] = nSpeciesPhaseCS
    db['nMagParamPhaseCS'] = nMagParamPhaseCS
    db['dAtomicMass'] = dAtomicMass
    db['cRegularParamCS'] = cRegularParamCS
    db['cElementNameCS'] = cElementNameCS
    db['cSolnPhaseTypeCS'] = cSolnPhaseTypeCS
    db['cSolnPhaseNameCS'] = cSolnPhaseNameCS
    db['cSpeciesNameCS'] = cSpeciesNameCS
    db['nPairsSROCS'] = nPairsSROCS
    db['nConstituentSublatticeCS'] = nConstituentSublatticeCS
    db['nSublatticeElementsCS'] = nSublatticeElementsCS
    db['dGibbsMagneticCS'] = dGibbsMagneticCS
    db['dStoichSublatticeCS'] = dStoichSublatticeCS
    db['dZetaSpeciesCS'] = dZetaSpeciesCS
    db['dGibbsCoeffSpeciesTemp'] = dGibbsCoeffSpeciesTemp
    db['dStoichSpeciesCS'] = dStoichSpeciesCS
    db['iMagneticParamCS'] = iMagneticParamCS
    db['dMagneticParamCS'] = dMagneticParamCS
    db['iRegularParamCS'] = iRegularParamCS
    db['dRegularParamCS'] = dRegularParamCS
    db['cPairNameCS'] = cPairNameCS
    db['iConstituentSublatticeCS'] = iConstituentSublatticeCS
    db['iPairIDCS'] = iPairIDCS
    db['iChemicalGroupCS'] = iChemicalGroupCS
    db['dSublatticeChargeCS'] = dSublatticeChargeCS
    db['dStoichPairsCS'] = dStoichPairsCS
    db['dConstituentCoefficientsCS'] = dConstituentCoefficientsCS
    db['dCoordinationNumberCS'] = dCoordinationNumberCS
    db['cConstituentNameSUBCS'] = cConstituentNameSUBCS
    db['cPhaseNames'] = cPhaseNames
    db['nPureSpeciesCS'] = nPureSpeciesCS
    db['nSolnPhaseCS'] = nSolnPhaseCS
    # db['indx'] = indx
    db['INFO'] = INFO
    db['DataBase'] = DataBase
    db['iCounterGibbsEqn'] = iCounterGibbsEqn
    db['cPeriodicTable'] = cPeriodicTable
    return db
