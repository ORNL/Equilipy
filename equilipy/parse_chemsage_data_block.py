"""Parse ChemSage data block sections."""

from __future__ import annotations

import numpy as np

import equilipy.variables as var

from .parse_chemsage_data_block_gibbs import ParseCSDataBlockGibbs
from .parse_chemsage_data_block_qkto import ParseCSDataBlockQKTO
from .parse_chemsage_data_block_rkmp import ParseCSDataBlockRKMP
from .parse_chemsage_data_block_subg import ParseCSDataBlockSUBG
from .parse_chemsage_data_block_subl import ParseCSDataBlockSUBL

_SUBLATTICE_MODELS = {"SUBL", "SUBLM", "SUBOM", "SUBG", "SUBQ"}
_MAGNETIC_HEADER_MODELS = {"RKMPM", "QKTOM", "SUBLM", "SUBOM"}
_PAIR_FRACTION_MODELS = {"SUBG", "SUBQ"}


def ParseCSDataBlock():
    """Parse the solution and pure-species data block."""
    # Intialize variables:
    var.nCountSublatticeCS = int(0)
    var.iCounterGibbsEqn = int(0)
    var.nParamCS = int(0)
    var.nMagParamCS = int(0)
    var.dTempVec = np.zeros(15, dtype=float)
    var.dGibbsMagneticCS = np.zeros((var.nSpeciesCS, 4), dtype=float)

    # Loop through all solution phases
    for i in range(1, var.nSolnPhasesSysCS + 1):
        # Read phase name
        var.cSolnPhaseNameCS[i - 1] = "{:<25}".format(var.DataBase.pop(0))

        # Read model name
        var.cSolnPhaseTypeCS[i - 1] = "{:<8}".format(var.DataBase.pop(0))
        phase_type = var.cSolnPhaseTypeCS[i - 1].strip().upper()

        if phase_type == "SUBOM":
            var.cDisorderedPhaseNameCS[i - 1] = str(var.DataBase.pop(0)).strip()

        # Check if the model is supported
        if var.cSolnPhaseTypeCS[i - 1] in var.cSolnPhaseTypeSupport:
            var.INFO = 0

        else:
            var.INFO = 17

        # Dummy parameter
        if i <= 1:
            var.indx = int(0) + 1
        else:
            var.indx = var.nSpeciesPhaseCS[i - 2] + 1

        # Count sublattice phases
        if phase_type in _SUBLATTICE_MODELS:
            var.nCountSublatticeCS = var.nCountSublatticeCS + 1
            var.iPhaseSublatticeCS[i - 1] = var.nCountSublatticeCS

        # Read magnetic terms if necessary
        if phase_type in _MAGNETIC_HEADER_MODELS:
            # Magnetic ordering is considered
            var.dTempVec[:2] = np.asarray(var.DataBase[:2], dtype=float)
            del var.DataBase[:2]
            for j in range(var.indx, var.nSpeciesPhaseCS[i - 1] + 1):
                var.dGibbsMagneticCS[j - 1, 2:4] = var.dTempVec[:2]
        elif phase_type == "SUBG":
            # Quadruplet model
            # Read zeta (FNN/SNN ratio)
            # In SUBG zeta is the same for every FNN pair
            # This method for keeping track of zeta is not very safe,
            # relies on having a predictable FNN pair order in the dat file
            dummy = float(var.DataBase.pop(0))
            var.dZetaSpeciesCS[var.nCountSublatticeCS - 1, : var.nMaxSpeciesPhaseCS] = (
                dummy * np.ones(var.nMaxSpeciesPhaseCS)
            )
            # Read integers for the number of species and the number of pairs.
            var.nPairsSROCS[var.nCountSublatticeCS - 1, :2] = np.asarray(
                var.DataBase[:2], dtype=int
            )
            del var.DataBase[:2]
        elif phase_type == "SUBQ":
            # MQM1 model
            # SUBQ phase data files seem not to have the magnetic term,
            # so skip this part.
            # Read integers for the number of species and the number of pairs.
            var.nPairsSROCS[var.nCountSublatticeCS - 1, :2] = np.asarray(
                var.DataBase[:2], dtype=int
            )
            del var.DataBase[:2]

        # Loop through species in solution phase:
        for j in range(var.indx, var.nSpeciesPhaseCS[i - 1] + 1):
            # SUBG and SUBQ phases contain a certain number of species,
            # necessarily fewer than the number of pair fractions. The # of
            # species indicated in the
            # header file actually represents the number of pairs. Therefore, there are
            # fewer species listed than what has been allocated.
            if phase_type in _PAIR_FRACTION_MODELS:
                if j >= var.indx + var.nPairsSROCS[var.nCountSublatticeCS - 1, 0]:
                    # Assign ID for all pairs
                    var.iPhaseCS[j - 1 : var.nSpeciesPhaseCS[i - 1] + 1] = i
                    break
            # Store the magnetic ordering terms for each solution:
            k = var.indx

            var.dGibbsMagneticCS[j - 1, 2:4] = var.dGibbsMagneticCS[k - 1, 2:4]

            # Store the phase index corresponding to the current species:
            var.iPhaseCS[j - 1] = i

            if phase_type in _PAIR_FRACTION_MODELS:
                # Parse the Gibbs energy equations (entries 3-5).
                ParseCSDataBlockGibbs(i, j)

                # Get pair stoichiometry in terms of constituents
                var.dConstituentCoefficientsCS[
                    var.nCountSublatticeCS - 1, j - var.indx, :5
                ] = np.asarray(var.DataBase[:5], dtype=float)
                del var.DataBase[:5]
                if phase_type == "SUBQ":
                    var.dZetaSpeciesCS[var.nCountSublatticeCS - 1, j - var.indx] = (
                        float(var.DataBase.pop(0))
                    )
            else:
                # Parse the Gibbs energy equations (entries 3-5).
                ParseCSDataBlockGibbs(i, j)
            if "QKTO" in var.cSolnPhaseTypeCS[i - 1]:
                var.dTempVec = np.asarray(var.DataBase[:2], dtype=float)
                del var.DataBase[:2]
                if "QKTOM" in var.cSolnPhaseTypeCS[i - 1]:
                    if sum(var.dTempVec - [1, 2]) > 1e-10:
                        var.INFO = 1600 + i
                        print("Error: QKTOM")

        # Check the type of solution phase to interpret mixing parameters:

        # Ideal mixing (IDMX)
        if var.cSolnPhaseTypeCS[i - 1] == var.cSolnPhaseTypeSupport[0]:
            pass
        # Quasichemical Kohler-Toop model
        elif "QKTO" in var.cSolnPhaseTypeCS[i - 1]:
            ParseCSDataBlockQKTO(i)
        # Redlich-Kister-Muggiano-Polynomial model
        elif "RKMP" in var.cSolnPhaseTypeCS[i - 1]:
            ParseCSDataBlockRKMP(i)
        # Compound Energy Formalism (sublattice) model
        elif phase_type in {"SUBL", "SUBLM", "SUBOM"}:
            ParseCSDataBlockSUBL(i)
        # Quadruplet quasichemical model:
        elif "SUBG" in var.cSolnPhaseTypeCS[i - 1]:
            ParseCSDataBlockSUBG(i)
        # Modified quasichemical model:
        elif "SUBQ" in var.cSolnPhaseTypeCS[i - 1]:
            ParseCSDataBlockSUBG(i)
        else:
            # The solution phase type is not supported. Report an error.
            var.INFO = 17
            return
        # Record the index of the mixing parameter for this phase:
        var.nParamPhaseCS[i - 1] = var.nParamCS
        var.nMagParamPhaseCS[i - 1] = var.nMagParamCS

    # END of Loop_SolnPhases
    phase_indices_by_name: dict[str, list[int]] = {}
    for index, name in enumerate(var.cSolnPhaseNameCS, start=1):
        if name.strip():
            phase_indices_by_name.setdefault(name.strip().upper(), []).append(index)
    occurrence_by_pair: dict[tuple[str, str], int] = {}
    for index, phase_name in enumerate(var.cDisorderedPhaseNameCS):
        if phase_name.strip():
            ordered_name = var.cSolnPhaseNameCS[index].strip().upper()
            disordered_name = phase_name.strip().upper()
            pair = (ordered_name, disordered_name)
            occurrence = occurrence_by_pair.get(pair, 0)
            disordered_indices = phase_indices_by_name.get(disordered_name, [])
            if disordered_indices:
                var.iDisorderedPhaseCS[index] = disordered_indices[
                    min(occurrence, len(disordered_indices) - 1)
                ]
            occurrence_by_pair[pair] = occurrence + 1

    # Begin parsing pure condensed phases:
    for j in range(
        var.nSpeciesPhaseCS[var.nSolnPhasesSysCS - 1] + 1, var.nSpeciesCS + 1
    ):
        # The phase index of a pure separate phase is set to zero:
        var.iPhaseCS[j - 1] = 0
        ParseCSDataBlockGibbs(var.iPhaseCS[j - 1], j)

    return None
    # Nest functions
