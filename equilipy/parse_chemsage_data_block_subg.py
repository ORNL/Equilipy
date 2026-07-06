"""Parse SUBG and SUBQ solution model data blocks."""

from __future__ import annotations

import numpy as np
import regex as re

import equilipy.variables as var


def ParseCSDataBlockSUBG(i):
    """Parse SUBG or SUBQ parameters for one solution phase."""
    # Modify from here

    def ParseConstituentName(k, constituent_name):
        # Check if there is charge number

        if "+" in constituent_name or "-" in constituent_name:
            # split by square brackets
            cNameL = list(filter(None, re.split(r"[\[\]]", constituent_name)))
            cNameList = list(filter(None, re.split(r"(\d+)", cNameL[0])))
        else:
            cNameList = list(filter(None, re.split(r"(\d+)", constituent_name)))

        for e in range(len(cNameList)):
            if cNameList[e].strip() == "Va":
                return
            iNumber = 1
            try:
                iNumber = int(cNameList[e])
                iElement = var.cElementNameCS.index(
                    "{:<3}".format(cNameList[e - 1].strip())
                )
            except ValueError:
                iNumber = 1
                iElement = var.cElementNameCS.index(
                    "{:<3}".format(cNameList[e].strip())
                )
            var.dStoichConstituentCS[k - 1, iElement] = float(iNumber)
        return None

    nCSCS = var.nCountSublatticeCS
    # This line contains N integers (where N is the number of sublattices)
    # where each integer represents the number of constituents on the respective
    # sublattice. There are always two sublattices for SUBG phases.
    var.nSublatticeElementsCS[nCSCS - 1, :2] = np.asarray(var.DataBase[:2], dtype=int)
    del var.DataBase[:2]
    var.nConstituentSublatticeCS[nCSCS - 1, :2] = var.nSublatticeElementsCS[
        nCSCS - 1, :2
    ]
    var.nSublatticePhaseCS[nCSCS - 1] = 2
    nTotalConst = (
        var.nConstituentSublatticeCS[nCSCS - 1, 0]
        + var.nConstituentSublatticeCS[nCSCS - 1, 1]
    )
    var.dStoichConstituentCS = np.zeros((nTotalConst, var.nElementsCS))

    nPairs = (
        var.nSublatticeElementsCS[nCSCS - 1, 0]
        * var.nSublatticeElementsCS[nCSCS - 1, 1]
    )
    # Read in names of constituents on first sublattice:
    var.cConstituentNameSUBCS[
        nCSCS - 1, 0, : var.nSublatticeElementsCS[nCSCS - 1, 0], :
    ] = np.asarray(
        [
            "{:<8}".format(x)
            for x in var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 0]]
        ],
        np.dtype("c"),
    )
    del var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 0]]

    # Set stoichiometry matrix. The order must follow var.cElementNameCS
    for k in range(1, var.nSublatticeElementsCS[nCSCS - 1, 0] + 1):
        # This block works similar to ParseConstituentName(k,1)
        cName = (
            var.cConstituentNameSUBCS[nCSCS - 1, 0, k - 1, :].tobytes().decode().strip()
        )
        ParseConstituentName(k, cName)

    # Read in names of constituents on second sublattice: (ignore for now):
    var.cConstituentNameSUBCS[
        nCSCS - 1, 1, : var.nSublatticeElementsCS[nCSCS - 1, 1], :
    ] = np.asarray(
        [
            "{:<8}".format(x)
            for x in var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 1]]
        ],
        np.dtype("c"),
    )
    del var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 1]]

    # Match elements on 2nd sublattice with elements in dat file order
    for k in range(1, var.nSublatticeElementsCS[nCSCS - 1, 1] + 1):
        # This block works similar to ParseConstituentName(k,1)
        cName = (
            var.cConstituentNameSUBCS[nCSCS - 1, 1, k - 1, :].tobytes().decode().strip()
        )
        k = k + var.nSublatticeElementsCS[nCSCS - 1, 0]
        ParseConstituentName(k, cName)

    # Read in the charge of each constituent on the first sublattice.
    var.dSublatticeChargeCS[nCSCS - 1, 0, : var.nSublatticeElementsCS[nCSCS - 1, 0]] = (
        np.asarray(var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 0]], dtype=float)
    )
    del var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 0]]

    # Chemical groups on sublattice 1:
    var.iChemicalGroupCS[nCSCS - 1, 0, : var.nSublatticeElementsCS[nCSCS - 1, 0]] = (
        np.asarray(var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 0]], dtype=int)
    )
    del var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 0]]

    # Read in the charge of each constituent on the second sublattice.
    var.dSublatticeChargeCS[nCSCS - 1, 1, : var.nSublatticeElementsCS[nCSCS - 1, 1]] = (
        np.asarray(var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 1]], dtype=float)
    )
    del var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 1]]

    # Chemical groups on sublattice 2:
    var.iChemicalGroupCS[nCSCS - 1, 1, : var.nSublatticeElementsCS[nCSCS - 1, 1]] = (
        np.asarray(var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 1]], dtype=int)
    )
    del var.DataBase[: var.nSublatticeElementsCS[nCSCS - 1, 1]]

    # This entry appears to represent IDs matching constituents on the first
    # sublattice to species:
    var.iConstituentSublatticeCS[nCSCS - 1, 0, :nPairs] = np.asarray(
        var.DataBase[:nPairs], dtype=int
    )
    del var.DataBase[:nPairs]

    # This entry appears to represent IDs matching constituents on the second
    # sublattice to species:
    var.iConstituentSublatticeCS[nCSCS - 1, 1, :nPairs] = np.asarray(
        var.DataBase[:nPairs], dtype=int
    )
    del var.DataBase[:nPairs]

    # Set up default pair IDs and coordination numbers
    for y in range(1, var.nSublatticeElementsCS[nCSCS - 1, 1] + 1):
        # LOOP_sroPairsOuter
        for x in range(1, var.nSublatticeElementsCS[nCSCS - 1, 1] + 1):
            if x == y:
                p = (x - 1) * int(
                    var.nSublatticeElementsCS[nCSCS - 1, 0]
                    * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                    / 2
                )
            elif x > y:
                continue
            else:
                p = (
                    var.nSublatticeElementsCS[nCSCS - 1, 1]
                    + (x - 1)
                    + int(((y - 2) * (y - 1) / 2))
                ) * int(
                    var.nSublatticeElementsCS[nCSCS - 1, 0]
                    * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                    / 2
                )
            for k in range(1, var.nSublatticeElementsCS[nCSCS - 1, 0] + 1):
                for j in range(1, var.nSublatticeElementsCS[nCSCS - 1, 0] + 1):
                    if j == k:
                        pair_index = j
                    elif j > k:
                        continue
                    else:
                        pair_index = (
                            var.nSublatticeElementsCS[nCSCS - 1, 0]
                            + j
                            + int((k - 2) * (k - 1) / 2)
                        )

                    # var.iPairIDCS is defined using Python numeric sequence:
                    # e.g. (1,1,4,4)  -> [0,0,3,3]
                    var.iPairIDCS[nCSCS - 1, pair_index + p - 1, 0] = j
                    var.iPairIDCS[nCSCS - 1, pair_index + p - 1, 1] = k
                    var.iPairIDCS[nCSCS - 1, pair_index + p - 1, 2] = (
                        x + var.nSublatticeElementsCS[nCSCS - 1, 0]
                    )
                    var.iPairIDCS[nCSCS - 1, pair_index + p - 1, 3] = (
                        y + var.nSublatticeElementsCS[nCSCS - 1, 0]
                    )

    # Parse the co-ordination numbers corresponding to all pairs in the phase.
    # Since these lines correspond to pairs, there will always be the same
    # number of integers and reals on a line, but the number of lines
    # corresponds to the number of pairs.
    # The SUBG model considers quadruplets, which is why there are four sets.
    # Note that a quadruplet must satisfy the following constraint:
    # q(i)/Z(i) + q(j)/Z(j) =  q(x)/Z(x) + q(y)/Z(y)
    lPairSet = np.zeros(var.nSpeciesPhaseCS[i - 1] - var.indx + 1, dtype=bool)

    for _ in range(1, var.nPairsSROCS[nCSCS - 1, 1] + 1):
        j = int(var.DataBase.pop(0))
        k = int(var.DataBase.pop(0))
        x = int(var.DataBase.pop(0))
        y = int(var.DataBase.pop(0))
        dTempVec = np.asarray(var.DataBase[:4], dtype=float)
        del var.DataBase[:4]
        x = x - int(var.nSublatticeElementsCS[nCSCS - 1, 0])
        y = y - int(var.nSublatticeElementsCS[nCSCS - 1, 0])
        if x == y:
            p = (x - 1) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
        elif x > y:
            continue
        else:
            p = (
                var.nSublatticeElementsCS[nCSCS - 1, 1]
                + (x - 1)
                + int(((y - 2) * (y - 1) / 2))
            ) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
        #  p = (var.nSublatticeElementsCS[nCSCS-1,1] + (x - 1)
        #       + int((y-2)*(y-1)/2)) * int(
        #           var.nSublatticeElementsCS[nCSCS-1,0]
        #           * (var.nSublatticeElementsCS[nCSCS-1,0] + 1) / 2
        #       )

        if j == k:
            pair_index = j
        elif j > k:
            continue
        else:
            pair_index = (
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                + j
                + int(((k - 2) * (k - 1) / 2))
            )

        var.dCoordinationNumberCS[nCSCS - 1, pair_index + p - 1, 0] = dTempVec[0]
        var.dCoordinationNumberCS[nCSCS - 1, pair_index + p - 1, 1] = dTempVec[1]
        var.dCoordinationNumberCS[nCSCS - 1, pair_index + p - 1, 2] = dTempVec[2]
        var.dCoordinationNumberCS[nCSCS - 1, pair_index + p - 1, 3] = dTempVec[3]
        lPairSet[pair_index + p - 1] = True

    # Increase pairs counter to include default pairs
    var.nPairsSROCS[nCSCS - 1, 1] = var.nSpeciesPhaseCS[i - 1] - var.indx + 1

    # This loop sets default coordination numbers for quadruplets not explicitly
    # listed in data file.
    for k in range(1, var.nPairsSROCS[nCSCS - 1, 1] + 1):
        # If coordinations already set, skip rest
        if lPairSet[k - 1]:
            continue

        # Constituent indices:
        a = var.iPairIDCS[nCSCS - 1, k - 1, 0]
        b = var.iPairIDCS[nCSCS - 1, k - 1, 1]
        x = var.iPairIDCS[nCSCS - 1, k - 1, 2] - var.nSublatticeElementsCS[nCSCS - 1, 0]
        y = var.iPairIDCS[nCSCS - 1, k - 1, 3] - var.nSublatticeElementsCS[nCSCS - 1, 0]

        # Constituent charges
        qa = var.dSublatticeChargeCS[nCSCS - 1, 0, a - 1]
        qb = var.dSublatticeChargeCS[nCSCS - 1, 0, b - 1]
        qx = var.dSublatticeChargeCS[nCSCS - 1, 1, x - 1]
        qy = var.dSublatticeChargeCS[nCSCS - 1, 1, y - 1]

        if a != b and x == y:
            p = (x - 1) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
            za = var.dCoordinationNumberCS[nCSCS - 1, p + a - 1, 0]
            zb = var.dCoordinationNumberCS[nCSCS - 1, p + b - 1, 0]

            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 0] = za
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 1] = zb
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 2] = (qx + qy) / (
                (qa / za) + (qb / zb)
            )
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 3] = (qx + qy) / (
                (qa / za) + (qb / zb)
            )
        elif a == b and x != y:
            p = (x - 1) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
            zx = var.dCoordinationNumberCS[nCSCS - 1, p + a - 1, 2]
            p = (y - 1) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
            zy = var.dCoordinationNumberCS[nCSCS - 1, p + a - 1, 2]

            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 0] = (qa + qb) / (
                (qx / zx) + (qy / zy)
            )
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 1] = (qa + qb) / (
                (qx / zx) + (qy / zy)
            )
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 2] = zx
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 3] = zy
        elif a != b and x != y:
            # Indices for AA/XY and BB/XY
            p = (
                var.nSublatticeElementsCS[nCSCS - 1, 1]
                + (x - 1)
                + int(((y - 2) * (y - 1) / 2))
            ) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
            iaaxy = a + p
            ibbxy = b + p

            # Indices for AB/XX and AB/YY
            ab_pair_index = (
                var.nSublatticeElementsCS[nCSCS - 1, 0] + a + int((b - 2) * (b - 1) / 2)
            )
            p = (x - 1) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
            iabxx = ab_pair_index + p
            p = (y - 1) * int(
                var.nSublatticeElementsCS[nCSCS - 1, 0]
                * (var.nSublatticeElementsCS[nCSCS - 1, 0] + 1)
                / 2
            )
            iabyy = ab_pair_index + p

            # Coordinations of specific species for the above quadruplets
            za = var.dCoordinationNumberCS[nCSCS - 1, iaaxy - 1, 0]
            zb = var.dCoordinationNumberCS[nCSCS - 1, ibbxy - 1, 0]
            zx = var.dCoordinationNumberCS[nCSCS - 1, iabxx - 1, 2]
            zy = var.dCoordinationNumberCS[nCSCS - 1, iabyy - 1, 2]

            # Equation 24 from part iv paper
            dF = (1 / 8) * ((qa / za) + (qb / zb) + (qx / zx) + (qy / zy))

            # Equation 23 from part iv paper
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 0] = 1 / (
                (
                    (zx / (qx * var.dCoordinationNumberCS[nCSCS - 1, iabxx - 1, 0]))
                    + (zy / (qy * var.dCoordinationNumberCS[nCSCS - 1, iabyy - 1, 0]))
                )
                * dF
            )
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 1] = 1 / (
                (
                    (zx / (qx * var.dCoordinationNumberCS[nCSCS - 1, iabxx - 1, 1]))
                    + (zy / (qy * var.dCoordinationNumberCS[nCSCS - 1, iabyy - 1, 1]))
                )
                * dF
            )
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 2] = 1 / (
                (
                    (zx / (qx * var.dCoordinationNumberCS[nCSCS - 1, iabxx - 1, 2]))
                    + (zy / (qy * var.dCoordinationNumberCS[nCSCS - 1, iabyy - 1, 2]))
                )
                * dF
            )
            var.dCoordinationNumberCS[nCSCS - 1, k - 1, 3] = 1 / (
                (
                    (zx / (qx * var.dCoordinationNumberCS[nCSCS - 1, iabxx - 1, 3]))
                    + (zy / (qy * var.dCoordinationNumberCS[nCSCS - 1, iabyy - 1, 3]))
                )
                * dF
            )

    # var.cPairNameCS[nCSCS-1,:var.nPairsSROCS[nCSCS-1,0],:] =\
    #   np.asarray(
    #       [
    #           "{:<30}".format(x)
    #           for x in var.cSpeciesNameCS[
    #               (var.indx-1):(var.indx+var.nPairsSROCS[nCSCS-1,0]-1)
    #           ]
    #       ],
    #       np.dtype("c"),
    #   )
    var.dStoichPairsCS[
        nCSCS - 1, : var.nPairsSROCS[nCSCS - 1, 1], : var.nElementsCS
    ] = var.dStoichSpeciesCS[
        (var.indx - 1) : var.nSpeciesPhaseCS[i - 1], : var.nElementsCS
    ]
    var.dStoichSpeciesCS[
        (var.indx - 1) : var.nSpeciesPhaseCS[i - 1], : var.nElementsCS
    ] = np.zeros((var.nSpeciesPhaseCS[i - 1] - var.indx + 1, var.nElementsCS))

    # Loop through all pairs to calculate stoichiometry entries for quadruplets:
    for j in range(1, var.nPairsSROCS[nCSCS - 1, 1] + 1):
        a = var.iPairIDCS[nCSCS - 1, j - 1, 0]
        b = var.iPairIDCS[nCSCS - 1, j - 1, 1]
        x = var.iPairIDCS[nCSCS - 1, j - 1, 2]
        y = var.iPairIDCS[nCSCS - 1, j - 1, 3]

        species_index = j + var.indx - 1
        # Use constituent stoichiometries to construct the species matrix.
        # To be clear, the species are the quadruplets.
        for k in range(1, var.nElementsCS + 1):
            var.dStoichSpeciesCS[species_index - 1, k - 1] = (
                var.dStoichSpeciesCS[species_index - 1, k - 1]
                + var.dStoichConstituentCS[a - 1, k - 1]
                / var.dCoordinationNumberCS[nCSCS - 1, j - 1, 0]
            )
            var.dStoichSpeciesCS[species_index - 1, k - 1] = (
                var.dStoichSpeciesCS[species_index - 1, k - 1]
                + var.dStoichConstituentCS[b - 1, k - 1]
                / var.dCoordinationNumberCS[nCSCS - 1, j - 1, 1]
            )
            var.dStoichSpeciesCS[species_index - 1, k - 1] = (
                var.dStoichSpeciesCS[species_index - 1, k - 1]
                + var.dStoichConstituentCS[x - 1, k - 1]
                / var.dCoordinationNumberCS[nCSCS - 1, j - 1, 2]
            )
            var.dStoichSpeciesCS[species_index - 1, k - 1] = (
                var.dStoichSpeciesCS[species_index - 1, k - 1]
                + var.dStoichConstituentCS[y - 1, k - 1]
                / var.dCoordinationNumberCS[nCSCS - 1, j - 1, 3]
            )

        # Create quadruplet names
        var.cSpeciesNameCS[species_index - 1] = "{:<30}".format(
            var.cConstituentNameSUBCS[nCSCS - 1, 0, a - 1, :].tobytes().decode().strip()
            + "-"
            + var.cConstituentNameSUBCS[nCSCS - 1, 0, b - 1].tobytes().decode().strip()
            + "-"
            + var.cConstituentNameSUBCS[
                nCSCS - 1, 1, x - var.nSublatticeElementsCS[nCSCS - 1, 0] - 1
            ]
            .tobytes()
            .decode()
            .strip()
            + "-"
            + var.cConstituentNameSUBCS[
                nCSCS - 1, 1, y - var.nSublatticeElementsCS[nCSCS - 1, 0] - 1
            ]
            .tobytes()
            .decode()
            .strip()
        )
        var.cEndmemberNameCS[species_index - 1] = var.cSpeciesNameCS[species_index - 1]
    # Copy previously-read end member info into appropriate variables before it
    # gets overwritten by
    # quadruplet data calculated below.
    var.cPairNameCS[nCSCS - 1, : var.nPairsSROCS[nCSCS - 1, 1], :] = np.asarray(
        [
            "{:<30}".format(x)
            for x in var.cSpeciesNameCS[
                (var.indx - 1) : (var.indx + var.nPairsSROCS[nCSCS - 1, 1] - 1)
            ]
        ],
        np.dtype("c"),
    )

    # Loop through excess mixing parameters:
    j = 0
    while True:
        j += 1
        # Read in number of constituents involved in parameter:
        var.iRegularParamCS[var.nParamCS, 0] = int(var.DataBase.pop(0))

        # The end of the parameter listing is marked by "0"
        # or a negative number indicating the number of extra parameter lines.
        # These lines indicate interpolation schemes. This parameters overwrite
        # the chemical group interpolation scheme.
        if var.iRegularParamCS[var.nParamCS, 0] <= 0:
            for _ in range(-var.iRegularParamCS[var.nParamCS, 0]):
                # This overwrites the interpolation scheme for the first sublattice.
                OverwriteInterpolationScheme(nCSCS, var.nParamCS)
            break

        # Check if the parameter is binary or ternary (3:Binary, 4:Ternary)
        if var.iRegularParamCS[var.nParamCS, 0] in [3, 4]:
            # Mixing terms
            var.cRegularParamCS[var.nParamCS] = str(var.DataBase.pop(0))[0]
            var.iRegularParamCS[var.nParamCS, 1:9] = np.asarray(
                var.DataBase[:8], dtype=int
            )
            del var.DataBase[:8]

            if var.cRegularParamCS[var.nParamCS] not in ["G", "Q", "R", "B", "H"]:
                var.INFO = 10000 + 1000 * j + i
                return print("Error: Excess Gibbs energy model is not supported")

            # According to Patrice Chartrand, he has no idea what these two
            # lines mean. Ignore.
            del var.DataBase[:12]

            # Read in the excess gibbs energy of mixing terms.
            var.iRegularParamCS[var.nParamCS, 9:11] = np.asarray(
                var.DataBase[:2], dtype=int
            )
            del var.DataBase[:2]
            var.dRegularParamCS[var.nParamCS, :6] = np.asarray(
                var.DataBase[:6], dtype=float
            )
            del var.DataBase[:6]

            # Count the number of parameters
            var.nParamCS += 1
        else:
            var.INFO = 10000 + 1000 * j + i
            return print("Error: Excess Gibbs energy parameters are not supported")

    return None


#
def OverwriteInterpolationScheme(n_cscs, n_param_cs):
    """
    Parse interpolation scheme and update var.iChemicalGroupCS.

    Format:
    [const1, const2, const3, (Tx, c1, c2), (K/Tx, c1, c2),
    (Tx, c1, c2), sublat2_const, ...]

    - First 3 integers: constituent IDs on first sublattice for this ternary parameter
    - Then groups of (interpolation_type, const1, const2) for each binary pair
    - 'Tn' = Toop with apex at ternary position n; that constituent is
      asymmetric
    - 'K' = Kohler → constituents are symmetric (must be in same group)
    - Last integer: constituent on second sublattice

    Example: ['2', '1', '3T2', '2', '1K', '2', '3T3', '1', '3', '4']
    Parsed as: [2, 1, 3] + [(T2, 2, 1), (K, 2, 3), (T3, 1, 3)] + [4]

    T2 means position 2 in ternary list (constituent '1') is asymmetric in pair (2,1)
    K between (2,3) means 2 and 3 must be in the same group (overrides T3)
    → Only Constituent 1 gets a different chemical group
    """
    import regex as re

    import equilipy.variables as var

    # Read the 10 entries for interpolation scheme
    entries = var.DataBase[:10]
    del var.DataBase[:10]

    # Parse entries to extract constituent IDs and interpolation markers
    constituentList = []  # First sublattice constituents for this ternary
    asymmetricConstituents = set()  # Constituents marked as asymmetric by Toop
    kohlerPairs = []  # Pairs of constituents that must be in same group (from Kohler)

    i = 0
    while i < len(entries):
        entry = str(entries[i]).strip()
        if not entry or entry == "0":
            i += 1
            continue

        # Check for interpolation marker (T or K) embedded in the string
        match = re.match(r"^(\d+)([TK])(\d+)?$", entry)

        if match:
            constId = int(match.group(1))
            interpType = match.group(2)  # 'T' or 'K'
            interpArg = int(match.group(3)) if match.group(3) else None

            # This constituent is part of a pair definition
            if len(constituentList) < 3:
                constituentList.append(constId)

            # Read the next two entries as the pair constituents
            i += 1
            pairConst1 = None
            pairConst2 = None
            if i < len(entries):
                try:
                    pairConst1 = int(entries[i])
                    i += 1
                except (TypeError, ValueError):
                    pass
            if i < len(entries):
                try:
                    # Check if this has another marker
                    nextMatch = re.match(r"^(\d+)([TK])?", str(entries[i]))
                    if nextMatch:
                        pairConst2 = int(nextMatch.group(1))
                        if not nextMatch.group(2):  # No marker, consume it
                            i += 1
                except (TypeError, ValueError):
                    pass

            if interpType == "T" and interpArg is not None:
                # Toop interpolation: the apex position is asymmetric.
                if interpArg <= len(constituentList):
                    apexConstId = constituentList[interpArg - 1]
                    asymmetricConstituents.add(apexConstId)
            elif interpType == "K":
                # Kohler: the two constituents in this pair must be in the same group
                if pairConst1 is not None and pairConst2 is not None:
                    kohlerPairs.append((pairConst1, pairConst2))
        else:
            # Pure integer - constituent ID
            try:
                constId = int(entry)
                if len(constituentList) < 3:
                    constituentList.append(constId)
                i += 1
            except ValueError:
                i += 1

    # Build union-find for Kohler constraints (constituents that must be in same group)
    parent = {}

    def find(x):
        if x not in parent:
            parent[x] = x
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Apply Kohler constraints
    for c1, c2 in kohlerPairs:
        union(c1, c2)

    # Determine which constituents are truly asymmetric
    # A constituent is asymmetric only if it is NOT in the same group as any
    # other constituent
    # that was NOT marked asymmetric

    # Find which asymmetric constituents are isolated, meaning not linked by
    # Kohler to non-asymmetric ones.
    trueAsymmetric = set()
    for constId in asymmetricConstituents:
        # Check whether this constituent is linked to any non-asymmetric
        # constituent via Kohler.
        linkedToSymmetric = False
        for other in constituentList:
            if other != constId and other not in asymmetricConstituents:
                if find(constId) == find(other):
                    linkedToSymmetric = True
                    break
        if not linkedToSymmetric:
            # Also check whether another asymmetric constituent is in the same
            # Kohler group.
            # If so, they should share their asymmetry (stay in same group)
            sameGroupAsymmetric = [
                c
                for c in asymmetricConstituents
                if c != constId and find(c) == find(constId)
            ]
            if not sameGroupAsymmetric:
                trueAsymmetric.add(constId)

    # Reset chemical groups for constituents involved in this scheme to 1.
    # This ensures that if a constituent was previously asymmetric (group 2+)
    # but is now symmetric, it gets reset.
    for constId in constituentList:
        var.iChemicalGroupCS[n_cscs - 1, 0, constId - 1] = 1

    # Update var.iChemicalGroupCS for truly asymmetric constituents
    nextGroup = 2
    for constId in sorted(list(trueAsymmetric)):  # Sort for deterministic order
        # constId is 1-indexed, array is 0-indexed
        var.iChemicalGroupCS[n_cscs - 1, 0, constId - 1] = nextGroup
        nextGroup += 1

    return None
