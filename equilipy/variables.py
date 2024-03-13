

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
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

cPeriodicTable=dict({
    'Ac': [89, 227.0],
    'Ag': [47, 107.8682],
    'Al': [13, 26.9815384],
    'Am': [95, 243.0],
    'Ar': [18, 39.948],
    'As': [33, 74.921595],
    'At': [85, 210.0],
    'Au': [79, 196.96657],
    'B': [5, 10.81],
    'Ba': [56, 137.327],
    'Be': [4, 9.0121831],
    'Bh': [107, 270.0],
    'Bi': [83, 208.9804],
    'Bk': [97, 247.0],
    'Br': [35, 79.904],
    'C': [6, 12.011],
    'Ca': [20, 40.078],
    'Cd': [48, 112.414],
    'Ce': [58, 140.116],
    'Cf': [98, 251.0],
    'Cl': [17, 35.45],
    'Cm': [96, 247.0],
    'Cn': [112, 285.0],
    'Co': [27, 58.933194],
    'Cr': [24, 51.9961],
    'Cs': [55, 132.905452],
    'Cu': [29, 63.546],
    'Db': [105, 270.0],
    'Ds': [110, 281.0],
    'Dy': [66, 162.5],
    'Er': [68, 167.259],
    'Es': [99, 252.0],
    'Eu': [63, 151.964],
    'F': [9, 18.99840316],
    'Fe': [26, 55.845],
    'Fl': [114, 289.0],
    'Fm': [100, 257.0],
    'Fr': [87, 223.0],
    'Ga': [31, 69.723],
    'Gd': [64, 157.25],
    'Ge': [32, 72.63],
    'H': [1, 1.008],
    'He': [2, 4.002602],
    'Hf': [72, 178.486],
    'Hg': [80, 200.592],
    'Ho': [67, 164.930328],
    'Hs': [108, 270.0],
    'I': [53, 126.90447],
    'In': [49, 114.818],
    'Ir': [77, 192.217],
    'K': [19, 39.0983],
    'Kr': [36, 83.798],
    'La': [57, 138.90547],
    'Li': [3, 6.94],
    'Lr': [103, 262.0],
    'Lu': [71, 174.9668],
    'Lv': [116, 293.0],
    'Mc': [115, 289.0],
    'Md': [101, 258.0],
    'Mg': [12, 24.305],
    'Mn': [25, 54.938043],
    'Mo': [42, 95.95],
    'Mt': [109, 278.0],
    'N': [7, 14.007],
    'Na': [11, 22.98976928],
    'Nb': [41, 92.90637],
    'Nd': [60, 144.242],
    'Ne': [10, 20.1797],
    'Nh': [113, 286.0],
    'Ni': [28, 58.6934],
    'No': [102, 259.0],
    'Np': [93, 237.0],
    'O': [8, 15.999],
    'Og': [118, 294.0],
    'Os': [76, 190.23],
    'P': [15, 30.973762],
    'Pa': [91, 231.03588],
    'Pb': [82, 207.2],
    'Pd': [46, 106.42],
    'Pm': [61, 145.0],
    'Po': [84, 209.0],
    'Pr': [59, 140.90766],
    'Pt': [78, 195.084],
    'Pu': [94, 244.0],
    'Ra': [88, 226.0],
    'Rb': [37, 85.4678],
    'Re': [75, 186.207],
    'Rf': [104, 267.0],
    'Rg': [111, 281.0],
    'Rh': [45, 102.90549],
    'Rn': [86, 222.0],
    'Ru': [44, 101.07],
    'S': [16, 32.06],
    'Sb': [51, 121.76],
    'Sc': [21, 44.955908],
    'Se': [34, 78.971],
    'Sg': [106, 269.0],
    'Si': [14, 28.085],
    'Sm': [62, 150.36],
    'Sn': [50, 118.71],
    'Sr': [38, 87.62],
    'Ta': [73, 180.94788],
    'Tb': [65, 158.925354],
    'Tc': [43, 97.0], 
    'Te': [52, 127.6],
    'Th': [90, 232.0377],
    'Ti': [22, 47.867],
    'Tl': [81, 204.38],
    'Tm': [69, 168.934218],
    'Ts': [117, 293.0],
    'U': [92, 238.02891],
    'V': [23, 50.9415],
    'W': [74, 183.84], 
    'Xe': [54, 131.293],
    'Y': [39, 88.90584],
    'Yb': [70, 173.045],
    'Zn': [30, 65.38],
    'Zr': [40, 91.224]
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
    db['indx'] = indx
    db['INFO'] = INFO
    db['DataBase'] = DataBase
    db['iCounterGibbsEqn'] = iCounterGibbsEqn
    db['cPeriodicTable'] = cPeriodicTable
    return db
