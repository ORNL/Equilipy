"""Runtime adapters from DatabaseIR to Equilipy calculation dictionaries."""

from __future__ import annotations

import re
from collections.abc import Collection
from dataclasses import dataclass, field, replace
from hashlib import sha1
from itertools import permutations, product

import numpy as np

from equilipy import variables as var
from equilipy.exceptions import DatabaseParsingError

from . import tdb as _tdb
from .model import DatabaseIR, FunctionDefinition, GibbsRange, Parameter, Phase
from .tdb_canonical import DISORDERED_PHASE_CANONICAL_NAMES

_PSEUDO_ELEMENTS = {"VA", "/-", "E", "E-", "ELECTRON_GAS"}
_NUMBER_RE = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?"
_T_LN_T_SENTINEL = 98.0
_LN_T_SENTINEL = 99.0
_RUNTIME_NAME_WIDTH = 30


@dataclass(frozen=True)
class _RuntimeCompound:
    name: str
    stoichiometry: dict[str, float]
    ranges: list[GibbsRange]
    phase_name: str = ""
    target: str = ""
    magnetic: tuple[float, float, float, float] | None = None


@dataclass(frozen=True)
class _RuntimeSolutionSpecies:
    name: str
    target: str
    stoichiometry: dict[str, float]
    ranges: list[GibbsRange]
    magnetic: tuple[float, float, float, float] | None = None
    sublattice_constituent_indices: tuple[int, ...] = ()


@dataclass(frozen=True)
class _RuntimeInteraction:
    arity: int
    species_indices: tuple[int, ...]
    order: int
    coefficients: np.ndarray
    ternary_weighted: bool = False


@dataclass(frozen=True)
class _RuntimeSolution:
    name: str
    phase_type: str
    species: list[_RuntimeSolutionSpecies]
    interactions: list[_RuntimeInteraction]
    magnetic_interactions: list[_RuntimeInteraction] = field(default_factory=list)
    sublattice_site_ratios: tuple[float, ...] = ()
    sublattice_constituents: tuple[tuple[str, ...], ...] = ()


@dataclass(frozen=True)
class _CompiledTargetPattern:
    """Compiled CEF target pattern for fast wildcard endmember resolution."""

    specificity: tuple[int, int]
    allowed_sections: tuple[frozenset[str], ...]
    parameter: Parameter


@dataclass
class _CefPhaseParameterResolver:
    """Per-phase index for CEF endmember and magnetic parameter lookup."""

    phase: Phase
    sublattice_constituents: tuple[tuple[str, ...], ...]
    exact_parameters: dict[str, dict[str, Parameter]]
    wildcard_g_parameters: list[_CompiledTargetPattern]
    symmetry_groups: list[list[int]]
    wildcard_g_by_target: dict[tuple[str, ...], Parameter] = field(
        default_factory=dict
    )
    symmetric_variants_cache: dict[tuple[str, ...], list[tuple[str, ...]]] = field(
        default_factory=dict
    )

    @classmethod
    def from_parameters(
        cls,
        phase: Phase,
        parameters: list[Parameter],
        sublattice_constituents: tuple[tuple[str, ...], ...],
    ) -> "_CefPhaseParameterResolver":
        """Compile exact, wildcard, and symmetry-aware parameter lookup data."""
        exact_parameters: dict[str, dict[str, Parameter]] = {}
        wildcard_g_parameters: list[_CompiledTargetPattern] = []
        for parameter in parameters:
            if parameter.order != 0 or not parameter.target:
                continue
            parameter_type = parameter.parameter_type.upper()
            target = parameter.target[0]
            normalized_target = _normalize_target(target)
            exact_by_target = exact_parameters.setdefault(parameter_type, {})
            exact_by_target.setdefault(normalized_target, parameter)
            if parameter_type != "G" or "," in target:
                continue
            pattern = _target_pattern_sections(target, sublattice_constituents)
            if pattern is None:
                continue
            wildcard_g_parameters.append(
                _CompiledTargetPattern(
                    specificity=_target_pattern_specificity(pattern),
                    allowed_sections=tuple(frozenset(section) for section in pattern),
                    parameter=parameter,
                )
            )
        wildcard_g_parameters.sort(
            key=lambda item: item.specificity,
            reverse=True,
        )
        resolver = cls(
            phase=phase,
            sublattice_constituents=sublattice_constituents,
            exact_parameters=exact_parameters,
            wildcard_g_parameters=wildcard_g_parameters,
            symmetry_groups=_runtime_equivalent_sublattice_groups(phase),
        )
        resolver.wildcard_g_by_target = resolver._build_wildcard_g_coverage()
        return resolver

    def gibbs_parameter(self, target: str) -> Parameter | None:
        """Return the G parameter covering ``target`` under TDB precedence."""
        parameter = self._exact_or_symmetric_parameter("G", target)
        if parameter is not None:
            return parameter
        return self.wildcard_g_by_target.get(self._target_sections(target))

    def _build_wildcard_g_coverage(self) -> dict[tuple[str, ...], Parameter]:
        """Expand wildcard G records to concrete normalized endmember targets."""
        coverage: dict[tuple[str, ...], Parameter] = {}
        for pattern in self.wildcard_g_parameters:
            section_options = [
                tuple(sorted(allowed)) for allowed in pattern.allowed_sections
            ]
            for covered_sections in product(*section_options):
                for variant in self._symmetric_variants(tuple(covered_sections)):
                    coverage.setdefault(variant, pattern.parameter)
        return coverage

    def magnetic_data(
        self,
        target: str,
        functions_by_name: dict[str, FunctionDefinition],
        phase_magnetic_options: dict[str, tuple[float, float]],
    ) -> tuple[float, float, float, float] | None:
        """Return runtime magnetic tuple for a CEF endmember target."""
        option = phase_magnetic_options.get(self.phase.name.upper())
        if option is None:
            return None
        tc_parameter = self._exact_or_symmetric_parameter("TC", target)
        moment_parameter = (
            self._exact_or_symmetric_parameter("BMAGN", target)
            or self._exact_or_symmetric_parameter("BMAG", target)
            or self._exact_or_symmetric_parameter("BM", target)
        )
        if tc_parameter is None and moment_parameter is None:
            return (0.0, 0.0, option[0], option[1])
        tc_value = (
            _constant_parameter_value(tc_parameter, functions_by_name)
            if tc_parameter is not None
            else 0.0
        )
        moment_value = (
            _constant_parameter_value(moment_parameter, functions_by_name)
            if moment_parameter is not None
            else 0.0
        )
        return (
            _runtime_magnetic_scalar_value(tc_value, option),
            _runtime_magnetic_scalar_value(moment_value, option),
            option[0],
            option[1],
        )

    def _exact_or_symmetric_parameter(
        self,
        parameter_type: str,
        target: str,
    ) -> Parameter | None:
        exact_by_target = self.exact_parameters.get(parameter_type.upper(), {})
        normalized_target = _normalize_target(target)
        exact = exact_by_target.get(normalized_target)
        if exact is not None:
            return exact
        for variant in self._symmetric_variants(self._target_sections(target)):
            variant_target = ":".join(variant)
            if variant_target == normalized_target:
                continue
            exact = exact_by_target.get(variant_target)
            if exact is not None:
                return exact
        return None

    def _target_sections(self, target: str) -> tuple[str, ...]:
        return tuple(_normalize_target(section) for section in target.split(":"))

    def _symmetric_variants(
        self,
        sections: tuple[str, ...],
    ) -> list[tuple[str, ...]]:
        cached = self.symmetric_variants_cache.get(sections)
        if cached is not None:
            return cached
        variants = _symmetric_section_variants(sections, self.symmetry_groups)
        self.symmetric_variants_cache[sections] = variants
        return variants


def to_chemsage_compounds(
    database: DatabaseIR,
    *,
    include_stoichiometric_compounds: bool = True,
    include_endmember_compounds: bool = False,
    endmember_compound_phases: Collection[str] | str | None = None,
) -> dict[str, object]:
    """Convert TDB records into a read_dat-style calculation dictionary."""
    element_symbols = _runtime_element_symbols(database)
    if not element_symbols:
        raise DatabaseParsingError("TDB does not define any runtime elements.")

    functions_by_name = {
        function.name.upper(): function for function in database.functions
    }
    disordered_phase_by_ordered = _phase_disordered_parts(database)
    phase_magnetic_options = _phase_magnetic_options(database)
    solutions = _runtime_solutions(
        database,
        element_symbols,
        functions_by_name,
        phase_magnetic_options,
        disordered_phase_by_ordered,
    )
    compounds = _runtime_compounds(
        database,
        element_symbols,
        include_stoichiometric_compounds=include_stoichiometric_compounds,
        include_endmember_compounds=include_endmember_compounds,
        endmember_compound_phases=endmember_compound_phases,
        functions_by_name=functions_by_name,
        phase_magnetic_options=phase_magnetic_options,
    )
    if not compounds and not solutions:
        raise DatabaseParsingError(
            "TDB does not define any executable PHASE records."
        )

    n_solution_species = sum(len(solution.species) for solution in solutions)
    n_species = n_solution_species + len(compounds)
    n_elements = len(element_symbols)
    n_soln_phases = len(solutions)
    soln_dim = max(1, n_soln_phases)
    regular_parameter_capacity = max(
        1000,
        sum(len(solution.interactions) for solution in solutions),
    )
    magnetic_parameter_capacity = max(
        1000,
        sum(len(solution.magnetic_interactions) for solution in solutions),
    )
    n_max_species_phase = max(
        [
            1,
            *(len(solution.species) for solution in solutions),
            *(
                len(constituents)
                for solution in solutions
                for constituents in solution.sublattice_constituents
            ),
        ]
    )
    n_max_gibbs_eqs = max(
        var.nMaxGibbsEqs,
        max(
            [
                *(
                    len(species.ranges)
                    for solution in solutions
                    for species in solution.species
                ),
                *(len(compound.ranges) for compound in compounds),
            ]
        ),
    )

    d_gibbs = np.zeros(
        (var.nGibbsCoeff, n_species * n_max_gibbs_eqs),
        dtype=float,
    )
    d_magnetic = np.zeros((n_species, 4), dtype=float)
    n_gibbs_eq = np.zeros(n_species, dtype=int)
    d_stoich = np.zeros((n_species, n_elements), dtype=float)
    i_phase = np.zeros(n_species, dtype=int)
    n_species_phase = np.zeros(soln_dim, dtype=int)
    n_soln_phase = np.zeros(soln_dim, dtype=int)
    n_param_phase = np.zeros(soln_dim, dtype=int)
    n_mag_param_phase = np.zeros(soln_dim, dtype=int)
    i_disordered_phase = np.zeros(soln_dim, dtype=int)
    i_regular_param = np.zeros(
        (regular_parameter_capacity, var.nParamMax * 2 + 3),
        dtype=int,
    )
    d_regular_param = np.zeros((regular_parameter_capacity, 6), dtype=float)
    i_magnetic_param = np.zeros(
        (magnetic_parameter_capacity, var.nParamMax * 2 + 3),
        dtype=int,
    )
    d_magnetic_param = np.zeros((magnetic_parameter_capacity, 2), dtype=float)
    i_phase_sublattice = np.zeros(soln_dim, dtype=int)
    n_sublattice_phase = np.zeros(soln_dim, dtype=int)
    n_constituent_sublattice = np.zeros(
        (soln_dim, var.nMaxSublatticeCS),
        dtype=int,
    )
    n_sublattice_elements = np.zeros(
        (soln_dim, var.nMaxSublatticeCS),
        dtype=int,
    )
    d_stoich_sublattice = np.zeros(
        (soln_dim, var.nMaxSublatticeCS),
        dtype=float,
    )
    i_constituent_sublattice = np.zeros(
        (soln_dim, var.nMaxSublatticeCS, n_max_species_phase),
        dtype=int,
    )
    c_constituent_name_subl = np.full(
        (soln_dim, var.nMaxSublatticeCS, n_max_species_phase, 8),
        b" ",
        dtype="c",
    )

    species_names: list[str] = []
    gibbs_column = 0
    species_index = 0
    param_index = 0
    mag_param_index = 0
    sublattice_record_count = 0
    for phase_index, solution in enumerate(solutions, start=1):
        n_soln_phase[phase_index - 1] = len(solution.species)
        n_species_phase[phase_index - 1] = species_index + len(solution.species)
        if solution.sublattice_constituents:
            sublattice_record_count += 1
            sublattice_index = sublattice_record_count - 1
            i_phase_sublattice[phase_index - 1] = sublattice_record_count
            n_sublattice_phase[sublattice_index] = len(
                solution.sublattice_constituents
            )
            d_stoich_sublattice[
                sublattice_index,
                : len(solution.sublattice_site_ratios),
            ] = solution.sublattice_site_ratios
            for set_index, constituents in enumerate(
                solution.sublattice_constituents
            ):
                n_constituent_sublattice[sublattice_index, set_index] = len(
                    constituents
                )
                n_sublattice_elements[sublattice_index, set_index] = len(
                    constituents
                )
                for constituent_index, constituent in enumerate(constituents):
                    c_constituent_name_subl[
                        sublattice_index,
                        set_index,
                        constituent_index,
                        :,
                    ] = np.asarray(f"{constituent:<8}"[:8], dtype="c")
        for solution_species in solution.species:
            i_phase[species_index] = phase_index
            species_names.append(_fixed_width(solution_species.name, 30))
            n_gibbs_eq[species_index] = len(solution_species.ranges)
            if solution_species.sublattice_constituent_indices:
                for set_index, constituent_index in enumerate(
                    solution_species.sublattice_constituent_indices
                ):
                    i_constituent_sublattice[
                        i_phase_sublattice[phase_index - 1] - 1,
                        set_index,
                        species_index - n_species_phase[phase_index - 2]
                        if phase_index > 1
                        else species_index,
                    ] = constituent_index
            if solution_species.magnetic is not None:
                d_magnetic[species_index, :] = solution_species.magnetic
            for element_index, symbol in enumerate(element_symbols):
                d_stoich[species_index, element_index] = (
                    solution_species.stoichiometry.get(symbol, 0.0)
                )
            for gibbs_range in solution_species.ranges:
                d_gibbs[:, gibbs_column] = _gibbs_coefficients(gibbs_range.Gibbs)
                d_gibbs[0, gibbs_column] = gibbs_range.T_max
                gibbs_column += 1
            species_index += 1
        for interaction in solution.interactions:
            if param_index >= i_regular_param.shape[0]:
                raise DatabaseParsingError(
                    "Runtime storage was allocated with too few solution "
                    "interaction parameter rows."
                )
            i_regular_param[param_index, 0] = interaction.arity
            i_regular_param[param_index, 1 : 1 + interaction.arity] = (
                interaction.species_indices
            )
            if solution.sublattice_constituents:
                i_regular_param[param_index, interaction.arity + 1] = (
                    interaction.order
                )
                if interaction.ternary_weighted:
                    i_regular_param[param_index, interaction.arity + 2] = 1
            elif interaction.arity == 2:
                i_regular_param[param_index, 3] = interaction.order
            elif interaction.arity == 3:
                order_index = min(interaction.order, interaction.arity - 1)
                i_regular_param[param_index, 4] = interaction.species_indices[
                    order_index
                ]
            elif interaction.arity == 4:
                i_regular_param[param_index, 5] = interaction.order
            d_regular_param[param_index, :] = interaction.coefficients
            param_index += 1
        n_param_phase[phase_index - 1] = param_index
        for interaction in solution.magnetic_interactions:
            if mag_param_index >= i_magnetic_param.shape[0]:
                raise DatabaseParsingError(
                    "Runtime storage was allocated with too few magnetic "
                    "solution interaction parameter rows."
                )
            i_magnetic_param[mag_param_index, 0] = interaction.arity
            i_magnetic_param[mag_param_index, 1 : 1 + interaction.arity] = (
                interaction.species_indices
            )
            if solution.sublattice_constituents:
                i_magnetic_param[mag_param_index, interaction.arity + 1] = (
                    interaction.order
                )
            elif interaction.arity == 2:
                i_magnetic_param[mag_param_index, 3] = interaction.order
            d_magnetic_param[mag_param_index, :] = interaction.coefficients[:2]
            mag_param_index += 1
        n_mag_param_phase[phase_index - 1] = mag_param_index

    for compound in compounds:
        species_names.append(_fixed_width(compound.name, 30))
        n_gibbs_eq[species_index] = len(compound.ranges)
        if compound.magnetic is not None:
            d_magnetic[species_index, :] = compound.magnetic
        for element_index, symbol in enumerate(element_symbols):
            d_stoich[species_index, element_index] = compound.stoichiometry.get(
                symbol,
                0.0,
            )
        for gibbs_range in compound.ranges:
            d_gibbs[:, gibbs_column] = _gibbs_coefficients(gibbs_range.Gibbs)
            d_gibbs[0, gibbs_column] = gibbs_range.T_max
            gibbs_column += 1
        species_index += 1

    element_by_upper = {
        element.symbol.upper(): element for element in database.elements
    }
    atomic_masses = np.asarray(
        [
            var.cPeriodicTable.get(symbol, [0, None])[1]
            if var.cPeriodicTable.get(symbol, [0, None])[1] is not None
            else (
                element_by_upper.get(symbol.upper()).atomic_mass
                if element_by_upper.get(symbol.upper()) is not None
                else 0.0
            )
            for symbol in element_symbols
        ],
        dtype=float,
    )
    solution_names = [_fixed_width(solution.name, 25) for solution in solutions]
    solution_types = [_fixed_width(solution.phase_type, 8) for solution in solutions]
    solution_index_by_name = {
        solution.name.upper(): index
        for index, solution in enumerate(solutions, start=1)
    }
    for ordered_name, disordered_name in disordered_phase_by_ordered.items():
        ordered_index = solution_index_by_name.get(ordered_name)
        disordered_index = solution_index_by_name.get(disordered_name)
        if ordered_index is not None and disordered_index is not None:
            i_disordered_phase[ordered_index - 1] = disordered_index

    return {
        "nElementsCS": n_elements,
        "nSpeciesCS": n_species,
        "nSolnPhasesSysCS": n_soln_phases,
        "nCountSublatticeCS": sublattice_record_count,
        "nMaxSpeciesPhaseCS": n_max_species_phase,
        "nParamCS": param_index,
        "nMagParamCS": mag_param_index,
        "nMaxSublatticeCS": var.nMaxSublatticeCS,
        "nSolnPhasesSysMax": var.nSolnPhasesSysMax,
        "nGibbsCoeff": var.nGibbsCoeff,
        "nMaxGibbsEqs": n_max_gibbs_eqs,
        "nParamMax": var.nParamMax,
        "iPhaseCS": i_phase,
        "iParticlesPerMoleCS": np.ones(n_species, dtype=int),
        "iPhaseSublatticeCS": i_phase_sublattice,
        "iDisorderedPhaseCS": i_disordered_phase,
        "nGibbsEqSpecies": n_gibbs_eq,
        "nSublatticePhaseCS": n_sublattice_phase,
        "nParamPhaseCS": n_param_phase,
        "nSpeciesPhaseCS": n_species_phase,
        "nMagParamPhaseCS": n_mag_param_phase,
        "dAtomicMass": atomic_masses,
        "cRegularParamCS": [_fixed_width("", 1)] * regular_parameter_capacity,
        "cElementNameCS": [_fixed_width(symbol, 3) for symbol in element_symbols],
        "cSolnPhaseTypeCS": [
            *(solution_types),
            *([_fixed_width("", 8)] * (soln_dim - n_soln_phases)),
        ],
        "cSolnPhaseNameCS": [
            *(solution_names),
            *([_fixed_width("", 25)] * (soln_dim - n_soln_phases)),
        ],
        "cSpeciesNameCS": species_names,
        "nPairsSROCS": np.zeros((soln_dim, 2), dtype=int),
        "nConstituentSublatticeCS": n_constituent_sublattice,
        "nSublatticeElementsCS": n_sublattice_elements,
        "dGibbsMagneticCS": d_magnetic,
        "dStoichSublatticeCS": d_stoich_sublattice,
        "dZetaSpeciesCS": np.zeros((soln_dim, n_max_species_phase), dtype=float),
        "dGibbsCoeffSpeciesTemp": d_gibbs,
        "dStoichSpeciesCS": d_stoich,
        "iMagneticParamCS": i_magnetic_param,
        "dMagneticParamCS": d_magnetic_param,
        "iRegularParamCS": i_regular_param,
        "dRegularParamCS": d_regular_param,
        "cPairNameCS": np.full((soln_dim, n_max_species_phase, 30), b" ", dtype="c"),
        "iConstituentSublatticeCS": i_constituent_sublattice,
        "iPairIDCS": np.zeros((soln_dim, n_max_species_phase, 4), dtype=int),
        "iChemicalGroupCS": np.zeros(
            (soln_dim, var.nMaxSublatticeCS, n_max_species_phase),
            dtype=int,
        ),
        "dSublatticeChargeCS": np.zeros(
            (soln_dim, var.nMaxSublatticeCS, n_max_species_phase),
            dtype=float,
        ),
        "dStoichPairsCS": np.zeros(
            (soln_dim, n_max_species_phase, n_elements),
            dtype=float,
        ),
        "dConstituentCoefficientsCS": np.zeros(
            (soln_dim, n_max_species_phase, 5),
            dtype=float,
        ),
        "dCoordinationNumberCS": np.zeros(
            (soln_dim, n_max_species_phase, 4),
            dtype=float,
        ),
        "cConstituentNameSUBCS": c_constituent_name_subl,
        "cPhaseNames": [
            *(solution.name for solution in solutions),
            *(compound.name for compound in compounds),
        ],
        "cOrderDisorderHelperPhaseNames": _duplicate_order_disorder_helper_names(
            solutions,
            disordered_phase_by_ordered,
        ),
        "cEndmemberNameCS": species_names.copy(),
        "nPureSpeciesCS": len(compounds),
        "nSolnPhaseCS": n_soln_phase,
        "indx": 0,
        "INFO": 0,
        "DataBase": [],
        "iCounterGibbsEqn": gibbs_column,
        "cPeriodicTable": var.cPeriodicTable,
    }


def _duplicate_order_disorder_helper_names(
    solutions: list[_RuntimeSolution],
    disordered_phase_by_ordered: dict[str, str],
) -> list[str]:
    """Return DIS_PART helper aliases that duplicate a canonical phase.

    TDB order/disorder phases may use helper names such as ``A2_BCC`` for the
    disordered reference of an ordered phase while also defining the physical
    phase ``BCC_A2``.  The helper must remain in the runtime arrays so the
    ordered-phase correction uses its own reference parameters, but it should
    not be advertised as an independent user-facing phase.
    """
    solution_names = {solution.name.upper(): solution.name for solution in solutions}
    helper_names: list[str] = []
    for disordered_name in disordered_phase_by_ordered.values():
        helper_key = disordered_name.upper()
        canonical_name = DISORDERED_PHASE_CANONICAL_NAMES.get(helper_key)
        if canonical_name is None:
            continue
        if helper_key not in solution_names or canonical_name not in solution_names:
            continue
        helper_names.append(solution_names[helper_key])
    return sorted(set(helper_names))


def _runtime_element_symbols(database: DatabaseIR) -> list[str]:
    symbols: list[str] = []
    for element in database.elements:
        if element.symbol.upper() in _PSEUDO_ELEMENTS:
            continue
        symbols.append(element.symbol)
    return symbols


def _runtime_solutions(
    database: DatabaseIR,
    element_symbols: list[str],
    functions_by_name: dict[str, FunctionDefinition],
    phase_magnetic_options: dict[str, tuple[float, float]],
    disordered_phase_by_ordered: dict[str, str],
) -> list[_RuntimeSolution]:
    """Return TDB solution phases as ChemSage runtime records."""
    parameters_by_phase: dict[str, list[Parameter]] = {}
    for parameter in database.parameters:
        parameters_by_phase.setdefault(parameter.phase_name.upper(), []).append(
            parameter
        )

    solutions: list[_RuntimeSolution] = []
    for phase in database.phases:
        parameters = parameters_by_phase.get(phase.name.upper(), [])
        if _phase_is_runtime_rkmp_solution(phase):
            species = _runtime_solution_species(
                phase,
                parameters,
                database,
                element_symbols,
                functions_by_name,
                phase_magnetic_options,
            )
            interactions = _runtime_solution_interactions(
                phase,
                parameters,
                {
                    item.target.upper(): index
                    for index, item in enumerate(species, start=1)
                },
                functions_by_name,
            )
            has_magnetic_species = any(item.magnetic is not None for item in species)
            if has_magnetic_species:
                phase_type = "RKMPM"
            elif interactions:
                phase_type = "RKMP"
            else:
                phase_type = "IDMX"
            solutions.append(
                _RuntimeSolution(
                    name=phase.name,
                    phase_type=phase_type,
                    species=species,
                    interactions=interactions,
                )
            )
        elif _phase_is_runtime_cef_solution(phase):
            solutions.append(
                _runtime_sublattice_solution(
                    phase,
                    parameters,
                    database,
                    element_symbols,
                    functions_by_name,
                    phase_magnetic_options,
                    disordered_phase_by_ordered,
                )
            )
    return solutions


def _phase_is_runtime_rkmp_solution(phase: Phase) -> bool:
    """Return whether a phase can be represented by the simple RKMP adapter."""
    if phase.model.upper() != "SOLUTION":
        return False
    if len(phase.constituents) != 1:
        return False
    species = [
        str(name).strip()
        for name in phase.constituents[0].species
        if str(name).strip() and str(name).strip() != "*"
    ]
    if len(species) <= 1:
        return False
    return not any(name.upper() in _PSEUDO_ELEMENTS for name in species)


def _phase_is_runtime_cef_solution(phase: Phase) -> bool:
    """Return whether a phase can be represented by ChemSage SUBL/SUBLM."""
    if phase.model.upper() != "CEF":
        return False
    if len(phase.constituents) < 2:
        return False
    if len(phase.constituents) > var.nMaxSublatticeCS:
        raise DatabaseParsingError(
            f"CEF phase '{phase.name}' has {len(phase.constituents)} "
            f"sublattices; runtime supports at most {var.nMaxSublatticeCS}."
        )
    return True


def _runtime_sublattice_solution(
    phase: Phase,
    parameters: list[Parameter],
    database: DatabaseIR,
    element_symbols: list[str],
    functions_by_name: dict[str, FunctionDefinition],
    phase_magnetic_options: dict[str, tuple[float, float]],
    disordered_phase_by_ordered: dict[str, str],
) -> _RuntimeSolution:
    """Build one ChemSage SUBL/SUBLM phase from a CEF TDB phase."""
    constituent_sets = sorted(phase.constituents, key=lambda item: item.sublattice)
    raw_sublattice_constituents = tuple(
        tuple(
            name
            for name in constituent_set.species
            if str(name).strip() and str(name).strip() != "*"
        )
        for constituent_set in constituent_sets
    )
    if any(not constituents for constituents in raw_sublattice_constituents):
        raise DatabaseParsingError(
            f"CEF phase '{phase.name}' has an empty sublattice constitution."
        )

    site_ratios = tuple(float(item.site_ratio) for item in constituent_sets)
    sublattice_constituents = raw_sublattice_constituents
    constituent_indices_by_sublattice = [
        {
            _normalize_target(name): index
            for index, name in enumerate(constituents, start=1)
        }
        for constituents in sublattice_constituents
    ]
    expected_targets = _cef_expected_runtime_targets(
        sublattice_constituents,
        phase,
        database,
        element_symbols,
    )
    parameter_resolver = _CefPhaseParameterResolver.from_parameters(
        phase,
        parameters,
        sublattice_constituents,
    )
    species: list[_RuntimeSolutionSpecies] = []
    for target in expected_targets:
        target_sections = target.split(":")
        parameter = parameter_resolver.gibbs_parameter(target)
        if parameter is None:
            parameter = _implicit_zero_endmember_parameter(phase, target)
        stoichiometry = _endmember_stoichiometry(
            phase,
            target,
            database,
            element_symbols,
        )
        species.append(
            _RuntimeSolutionSpecies(
                name=target,
                target=target,
                stoichiometry=stoichiometry,
                ranges=_expanded_parameter_ranges(parameter, functions_by_name),
                magnetic=parameter_resolver.magnetic_data(
                    target,
                    functions_by_name,
                    phase_magnetic_options,
                ),
                sublattice_constituent_indices=tuple(
                    constituent_indices_by_sublattice[index][_normalize_target(name)]
                    for index, name in enumerate(target_sections)
                ),
            ),
        )

    interactions = _runtime_sublattice_solution_interactions(
        phase,
        parameters,
        sublattice_constituents,
        functions_by_name,
    )
    has_magnetic_species = any(item.magnetic is not None for item in species)
    magnetic_interactions = _runtime_sublattice_magnetic_interactions(
        phase,
        parameters,
        sublattice_constituents,
        functions_by_name,
        phase_magnetic_options,
    )
    has_magnetic_terms = has_magnetic_species or bool(magnetic_interactions)
    if phase.name.upper() in disordered_phase_by_ordered:
        phase_type = "SUBOM"
    elif has_magnetic_terms:
        phase_type = "SUBLM"
    else:
        phase_type = "SUBL"
    return _RuntimeSolution(
        name=phase.name,
        phase_type=phase_type,
        species=species,
        interactions=interactions,
        magnetic_interactions=magnetic_interactions,
        sublattice_site_ratios=site_ratios,
        sublattice_constituents=sublattice_constituents,
    )


def _implicit_zero_endmember_parameter(phase: Phase, target: str) -> Parameter:
    """Return a runtime-only zero G record for an omitted CEF endmember.

    Commercial TDB readers commonly treat omitted CEF endmember G terms as zero
    compatibility records.  Equilipy should be able to execute those databases
    without materializing thousands of ``ZERO#`` parameters in the IR or writer
    output, so the fallback lives only in the runtime adapter.
    """
    return Parameter(
        phase_name=phase.name,
        parameter_type="G",
        target=[target],
        order=0,
        expression="298.15 0; 6000 N",
        source=phase.source,
        metadata={
            "generated_by": "tdb_runtime",
            "reason": "implicit_missing_endmember_zero",
        },
    )


def _cef_gibbs_parameter_for_target(
    parameters: list[Parameter],
    target: str,
    sublattice_constituents: tuple[tuple[str, ...], ...],
    phase: Phase | None = None,
) -> Parameter | None:
    exact = _target_parameter(parameters, "G", target)
    if exact is not None:
        return exact
    target_sections = [_normalize_target(section) for section in target.split(":")]
    symmetry_groups = _runtime_equivalent_sublattice_groups(phase)
    for variant in _symmetric_section_variants(target_sections, symmetry_groups):
        variant_target = ":".join(variant)
        if variant_target == target:
            continue
        exact = _target_parameter(parameters, "G", variant_target)
        if exact is not None:
            return exact
    wildcard_matches: list[tuple[tuple[int, int], Parameter]] = []
    for parameter in parameters:
        if parameter.parameter_type.upper() != "G":
            continue
        if parameter.order != 0 or not parameter.target:
            continue
        if "," in parameter.target[0]:
            continue
        pattern = _target_pattern_sections(
            parameter.target[0],
            sublattice_constituents,
        )
        if pattern is None:
            continue
        if _target_pattern_covers_symmetric(pattern, target_sections, symmetry_groups):
            wildcard_matches.append((_target_pattern_specificity(pattern), parameter))
    if wildcard_matches:
        wildcard_matches.sort(key=lambda item: item[0], reverse=True)
        return wildcard_matches[0][1]
    return None


def _single_constituent_target_sections(
    target: str,
    sublattice_constituents: tuple[tuple[str, ...], ...],
) -> list[str] | None:
    sections = [section.strip() for section in target.split(":")]
    if len(sections) != len(sublattice_constituents):
        return None
    canonical_sections: list[str] = []
    for section, constituents in zip(sections, sublattice_constituents, strict=True):
        tokens = [
            _normalize_target(token)
            for token in section.split(",")
            if token.strip()
        ]
        if len(tokens) != 1 or tokens[0] == "*":
            return None
        by_upper = {_normalize_target(name): name for name in constituents}
        canonical = by_upper.get(tokens[0])
        if canonical is None:
            return None
        canonical_sections.append(canonical)
    return canonical_sections


def _target_pattern_sections(
    target: str,
    sublattice_constituents: tuple[tuple[str, ...], ...],
) -> list[list[str]] | None:
    sections = [section.strip() for section in target.split(":")]
    if len(sections) != len(sublattice_constituents):
        return None
    patterns: list[list[str]] = []
    for section, constituents in zip(sections, sublattice_constituents, strict=True):
        tokens = [
            _normalize_target(token)
            for token in section.split(",")
            if token.strip()
        ]
        if not tokens:
            return None
        allowed = {_normalize_target(name) for name in constituents}
        if "*" in tokens:
            patterns.append(sorted(allowed))
            continue
        if not set(tokens).issubset(allowed):
            return None
        patterns.append(tokens)
    return patterns


def _target_pattern_covers(
    pattern: list[list[str]],
    target_sections: list[str],
) -> bool:
    if len(pattern) != len(target_sections):
        return False
    return all(
        target in allowed
        for allowed, target in zip(pattern, target_sections, strict=True)
    )


def _target_pattern_covers_symmetric(
    pattern: list[list[str]],
    target_sections: list[str],
    symmetry_groups: list[list[int]],
) -> bool:
    return any(
        _target_pattern_covers(pattern, list(variant))
        for variant in _symmetric_section_variants(target_sections, symmetry_groups)
    )


def _runtime_equivalent_sublattice_groups(phase: Phase | None) -> list[list[int]]:
    """Return zero-based ``:F``-equivalent sublattice groups for runtime matching."""
    return _tdb._phase_equivalent_sublattice_groups(phase)


def _symmetric_section_variants(
    sections: list[str] | tuple[str, ...],
    symmetry_groups: list[list[int]],
) -> list[tuple[str, ...]]:
    """Return unique section permutations allowed by ``PHASE ...:F`` symmetry."""
    variants = [tuple(sections)]
    for group in symmetry_groups:
        if any(index >= len(sections) for index in group):
            continue
        next_variants: list[tuple[str, ...]] = []
        seen_for_group: set[tuple[str, ...]] = set()
        for variant in variants:
            for permutation in permutations(
                [variant[index] for index in group],
                len(group),
            ):
                candidate = list(variant)
                for index, value in zip(group, permutation, strict=True):
                    candidate[index] = value
                key = tuple(candidate)
                if key in seen_for_group:
                    continue
                seen_for_group.add(key)
                next_variants.append(key)
        variants = next_variants

    unique: list[tuple[str, ...]] = []
    seen: set[tuple[str, ...]] = set()
    for variant in variants:
        if variant in seen:
            continue
        seen.add(variant)
        unique.append(variant)
    return unique


def _target_pattern_specificity(pattern: list[list[str]]) -> tuple[int, int]:
    exact_sections = sum(1 for allowed in pattern if len(allowed) == 1)
    total_width = sum(len(allowed) for allowed in pattern)
    return exact_sections, -total_width


def _cef_expected_runtime_targets(
    sublattice_constituents: tuple[tuple[str, ...], ...],
    _phase: Phase,
    _database: DatabaseIR,
    _element_symbols: list[str],
) -> list[str]:
    targets: list[str] = []
    for sections in product(*sublattice_constituents):
        targets.append(":".join(sections))
    return targets


def _runtime_sublattice_solution_interactions(
    phase: Phase,
    parameters: list[Parameter],
    sublattice_constituents: tuple[tuple[str, ...], ...],
    functions_by_name: dict[str, FunctionDefinition],
) -> list[_RuntimeInteraction]:
    interactions: list[_RuntimeInteraction] = []
    max_width = var.nParamMax * 2 + 3
    offset = 0
    constituent_lookup_by_sublattice: list[dict[str, int]] = []
    for constituents in sublattice_constituents:
        constituent_lookup_by_sublattice.append(
            {
                _normalize_target(name): offset + index
                for index, name in enumerate(constituents, start=1)
            }
        )
        offset += len(constituents)

    symmetry_groups = _runtime_equivalent_sublattice_groups(phase)
    for parameter in parameters:
        if parameter.parameter_type.upper() not in {"G", "L"} or not parameter.target:
            continue
        sections = parameter.target[0].split(":")
        if len(sections) != len(sublattice_constituents):
            continue
        section_option_variants = []
        for variant_sections in _symmetric_section_variants(sections, symmetry_groups):
            section_options = _sublattice_interaction_section_options(
                phase,
                parameter.target[0],
                list(variant_sections),
                constituent_lookup_by_sublattice,
            )
            if section_options is not None:
                section_option_variants.append(section_options)
        if not section_option_variants:
            continue
        base_coefficients: np.ndarray | None = None
        emitted_indices: set[tuple[int, ...]] = set()
        for section_options in section_option_variants:
            for option_sections in product(*section_options):
                mix_lengths = [
                    len(tokens) for tokens in option_sections if len(tokens) > 1
                ]
                if not mix_lengths:
                    continue
                if mix_lengths not in ([2], [3], [2, 2]):
                    raise DatabaseParsingError(
                        f"CEF phase '{phase.name}' has unsupported interaction "
                        f"target '{parameter.target[0]}'. Runtime SUBL supports "
                        "one binary set, one ternary set, or two binary sets."
                    )
                indices = tuple(index for tokens in option_sections for index in tokens)
                if indices in emitted_indices:
                    continue
                emitted_indices.add(indices)
                if len(indices) + 2 > max_width:
                    raise DatabaseParsingError(
                        f"CEF phase '{phase.name}' has interaction target "
                        f"'{parameter.target[0]}' with {len(indices)} constituents; "
                        f"runtime storage supports at most {max_width - 2}."
                    )
                if base_coefficients is None:
                    base_coefficients = _regular_parameter_coefficients(
                        parameter,
                        functions_by_name,
                    )
                # ``:F`` marks equivalent sublattices; each emitted variant
                # keeps the original TDB parameter value.
                interactions.append(
                    _RuntimeInteraction(
                        arity=len(indices),
                        species_indices=indices,
                        order=parameter.order,
                        coefficients=base_coefficients,
                    )
                )
    return _mark_complete_ternary_sublattice_interactions(
        interactions,
        sublattice_constituents,
    )


def _mark_complete_ternary_sublattice_interactions(
    interactions: list[_RuntimeInteraction],
    sublattice_constituents: tuple[tuple[str, ...], ...],
) -> list[_RuntimeInteraction]:
    """Mark TC-style complete ternary ``;0``/``;1``/``;2`` CEF sets.

    A lone TDB ternary ``;0`` record is an unweighted ``yA*yB*yC*L0`` term.
    When the same target is supplied with all three orders, TC evaluates the
    complete set with Muggianu redistribution in higher-component sublattices:
    ``yA*yB*yC*sum_i((y_i + (1 - yA - yB - yC)/3)*L_i)``.  The runtime stores
    this as the original three rows plus a marker bit so the Fortran evaluator
    can distinguish complete sets from lone ``;0`` records.
    """
    complete_keys = _complete_ternary_sublattice_interaction_keys(
        interactions,
        sublattice_constituents,
    )
    if not complete_keys:
        return interactions
    return [
        replace(interaction, ternary_weighted=True)
        if (
            interaction.order in (0, 1, 2)
            and tuple(interaction.species_indices) in complete_keys
        )
        else interaction
        for interaction in interactions
    ]


def _complete_ternary_sublattice_interaction_keys(
    interactions: list[_RuntimeInteraction],
    sublattice_constituents: tuple[tuple[str, ...], ...],
) -> set[tuple[int, ...]]:
    orders_by_key: dict[tuple[int, ...], set[int]] = {}
    for interaction in interactions:
        if not _is_single_section_ternary_interaction(
            interaction,
            sublattice_constituents,
        ):
            continue
        orders_by_key.setdefault(tuple(interaction.species_indices), set()).add(
            interaction.order,
        )
    return {key for key, orders in orders_by_key.items() if orders == {0, 1, 2}}


def _is_single_section_ternary_interaction(
    interaction: _RuntimeInteraction,
    sublattice_constituents: tuple[tuple[str, ...], ...],
) -> bool:
    if interaction.order not in (0, 1, 2):
        return False
    by_sublattice: dict[int, int] = {}
    for encoded in interaction.species_indices:
        remaining = int(encoded)
        matched_sublattice_index = 0
        for candidate_sublattice_index, constituents in enumerate(
            sublattice_constituents,
            start=1,
        ):
            width = len(constituents)
            if 1 <= remaining <= width:
                matched_sublattice_index = candidate_sublattice_index
                break
            remaining -= width
        else:
            return False
        by_sublattice[matched_sublattice_index] = (
            by_sublattice.get(matched_sublattice_index, 0) + 1
        )
    counts = sorted(by_sublattice.values())
    return counts.count(3) == 1 and all(count in (1, 3) for count in counts)


def _runtime_sublattice_magnetic_interactions(
    phase: Phase,
    parameters: list[Parameter],
    sublattice_constituents: tuple[tuple[str, ...], ...],
    functions_by_name: dict[str, FunctionDefinition],
    phase_magnetic_options: dict[str, tuple[float, float]],
) -> list[_RuntimeInteraction]:
    """Build paired TC/BM CEF magnetic interaction terms.

    TDB files commonly provide only ``TC`` interaction parameters.  The
    missing magnetic-moment partner is a zero contribution, not an instruction
    to drop the Curie/Neel-temperature interaction.
    """
    moment_by_key: dict[tuple[str, int], Parameter] = {}
    for parameter in parameters:
        if not parameter.target:
            continue
        parameter_type = parameter.parameter_type.upper()
        if parameter_type not in {"BMAGN", "BMAG", "BM"}:
            continue
        key = (_normalize_target(parameter.target[0]), parameter.order)
        if key not in moment_by_key or _magnetic_moment_priority(
            parameter_type,
        ) > _magnetic_moment_priority(moment_by_key[key].parameter_type):
            moment_by_key[key] = parameter

    interactions: list[_RuntimeInteraction] = []
    tc_keys: set[tuple[str, int]] = set()
    for tc_parameter in parameters:
        if tc_parameter.parameter_type.upper() != "TC" or not tc_parameter.target:
            continue
        key = (_normalize_target(tc_parameter.target[0]), tc_parameter.order)
        tc_keys.add(key)
        moment_parameter = moment_by_key.get(key)
        interactions.extend(
            _runtime_sublattice_paired_interactions(
                phase,
                tc_parameter,
                moment_parameter,
                sublattice_constituents,
                functions_by_name,
                phase_magnetic_options,
            )
        )
    for key, moment_parameter in moment_by_key.items():
        if key in tc_keys:
            continue
        interactions.extend(
            _runtime_sublattice_paired_interactions(
                phase,
                None,
                moment_parameter,
                sublattice_constituents,
                functions_by_name,
                phase_magnetic_options,
            )
        )
    return interactions


def _runtime_sublattice_paired_interactions(
    phase: Phase,
    first_parameter: Parameter | None,
    second_parameter: Parameter | None,
    sublattice_constituents: tuple[tuple[str, ...], ...],
    functions_by_name: dict[str, FunctionDefinition],
    phase_magnetic_options: dict[str, tuple[float, float]],
) -> list[_RuntimeInteraction]:
    interactions: list[_RuntimeInteraction] = []
    max_width = var.nParamMax * 2 + 3
    offset = 0
    constituent_lookup_by_sublattice: list[dict[str, int]] = []
    for constituents in sublattice_constituents:
        constituent_lookup_by_sublattice.append(
            {
                _normalize_target(name): offset + index
                for index, name in enumerate(constituents, start=1)
            }
        )
        offset += len(constituents)

    template_parameter = (
        first_parameter if first_parameter is not None else second_parameter
    )
    if template_parameter is None:
        return interactions
    sections = template_parameter.target[0].split(":")
    if len(sections) != len(sublattice_constituents):
        return interactions
    symmetry_groups = _runtime_equivalent_sublattice_groups(phase)
    section_option_variants = []
    for variant_sections in _symmetric_section_variants(sections, symmetry_groups):
        section_options = _sublattice_interaction_section_options(
            phase,
            template_parameter.target[0],
            list(variant_sections),
            constituent_lookup_by_sublattice,
        )
        if section_options is not None:
            section_option_variants.append(section_options)
    if not section_option_variants:
        return interactions

    base_coefficients = np.zeros(6, dtype=float)
    if first_parameter is not None:
        base_coefficients[0] = _constant_parameter_value(
            first_parameter,
            functions_by_name,
        )
    if second_parameter is not None:
        base_coefficients[1] = _constant_parameter_value(
            second_parameter,
            functions_by_name,
        )
    emitted_indices: set[tuple[int, ...]] = set()
    for section_options in section_option_variants:
        for option_sections in product(*section_options):
            mix_lengths = [len(tokens) for tokens in option_sections if len(tokens) > 1]
            if not mix_lengths:
                continue
            if mix_lengths not in ([2], [3], [2, 2]):
                raise DatabaseParsingError(
                    f"CEF phase '{phase.name}' has unsupported magnetic interaction "
                    f"target '{template_parameter.target[0]}'. Runtime SUBLM supports "
                    "one binary set, one ternary set, or two binary sets."
                )
            ordered_option_sections = _magnetic_option_sections_for_fortran(
                option_sections,
            )
            indices = tuple(
                index for tokens in ordered_option_sections for index in tokens
            )
            if indices in emitted_indices:
                continue
            emitted_indices.add(indices)
            if len(indices) + 2 > max_width:
                raise DatabaseParsingError(
                    f"CEF phase '{phase.name}' has magnetic interaction target "
                    f"'{template_parameter.target[0]}' with "
                    f"{len(indices)} constituents; "
                    f"runtime storage supports at most {max_width - 2}."
                )
            # ``:F`` marks equivalent sublattices; each emitted variant keeps
            # the original TC/BM parameter values.
            interactions.append(
                _RuntimeInteraction(
                    arity=len(indices),
                    species_indices=indices,
                    order=template_parameter.order,
                    coefficients=base_coefficients,
                )
            )
    return interactions


def _magnetic_option_sections_for_fortran(
    option_sections: tuple[tuple[int, ...], ...],
) -> tuple[tuple[int, ...], ...]:
    """Place the mixed sublattice first for SUBLM/SUBOM magnetic terms.

    ``CompGibbsMagneticSoln`` treats the first two encoded constituents as the
    binary magnetic pair and all later constituents as prefactors.  TDB targets
    may put a fixed sublattice before the mixed one, for example
    ``Al:Al,Fe:Va``.  Reorder only the emitted runtime indices; coefficient
    scaling is still computed from the original sublattice positions.
    """
    return tuple(
        [
            *[tokens for tokens in option_sections if len(tokens) > 1],
            *[tokens for tokens in option_sections if len(tokens) <= 1],
        ],
    )


def _sublattice_interaction_section_options(
    phase: Phase,
    target: str,
    sections: list[str],
    constituent_lookup_by_sublattice: list[dict[str, int]],
) -> list[list[tuple[int, ...]]] | None:
    section_options: list[list[tuple[int, ...]]] = []
    for section, lookup in zip(
        sections,
        constituent_lookup_by_sublattice,
        strict=True,
    ):
        tokens = [
            _normalize_target(token)
            for token in section.split(",")
            if token.strip()
        ]
        if not tokens:
            return None
        if "*" in tokens:
            section_options.append([(index,) for index in lookup.values()])
            continue
        try:
            section_options.append([tuple(lookup[token] for token in tokens)])
        except KeyError:
            return None
    return section_options


def _runtime_solution_species(
    phase: Phase,
    parameters: list[Parameter],
    database: DatabaseIR,
    element_symbols: list[str],
    functions_by_name: dict[str, FunctionDefinition],
    phase_magnetic_options: dict[str, tuple[float, float]],
) -> list[_RuntimeSolutionSpecies]:
    gibbs_parameters = _gibbs_parameters(parameters)
    species: list[_RuntimeSolutionSpecies] = []
    for name in phase.constituents[0].species:
        if name.upper() in _PSEUDO_ELEMENTS or name == "*":
            continue
        parameter = _target_parameter(gibbs_parameters, "G", name)
        if parameter is None:
            raise DatabaseParsingError(
                f"Solution phase '{phase.name}' is missing a G endmember "
                f"parameter for '{name}'."
            )
        species.append(
            _RuntimeSolutionSpecies(
                name=name,
                target=name,
                stoichiometry=_endmember_stoichiometry(
                    phase,
                    name,
                    database,
                    element_symbols,
                ),
                ranges=_expanded_parameter_ranges(parameter, functions_by_name),
                magnetic=_magnetic_data_for_target(
                    phase,
                    name,
                    parameters,
                    functions_by_name,
                    phase_magnetic_options,
                ),
            )
        )
    return species


def _runtime_solution_interactions(
    phase: Phase,
    parameters: list[Parameter],
    species_indices_by_name: dict[str, int],
    functions_by_name: dict[str, FunctionDefinition],
) -> list[_RuntimeInteraction]:
    interactions: list[_RuntimeInteraction] = []
    for parameter in parameters:
        if parameter.parameter_type.upper() not in {"G", "L"} or not parameter.target:
            continue
        sections = parameter.target[0].split(":")
        if len(sections) != 1:
            continue
        tokens = [_normalize_target(token) for token in sections[0].split(",")]
        tokens = [token for token in tokens if token]
        arity = len(tokens)
        if arity < 2:
            continue
        if arity > 4:
            raise DatabaseParsingError(
                f"Solution phase '{phase.name}' has unsupported L parameter "
                f"target '{parameter.target[0]}'. Runtime RKMP supports binary, "
                "ternary, and quaternary interactions."
            )
        try:
            indices = tuple(species_indices_by_name[token] for token in tokens)
        except KeyError as exc:
            raise DatabaseParsingError(
                f"Solution phase '{phase.name}' has L parameter target "
                f"'{parameter.target[0]}' for a species not present in the "
                "phase constitution."
            ) from exc
        interactions.append(
            _RuntimeInteraction(
                arity=arity,
                species_indices=indices,
                order=parameter.order,
                coefficients=_regular_parameter_coefficients(
                    parameter,
                    functions_by_name,
                ),
            )
        )
    return interactions


def _regular_parameter_coefficients(
    parameter: Parameter,
    functions_by_name: dict[str, FunctionDefinition],
) -> np.ndarray:
    ranges = _expanded_parameter_ranges(parameter, functions_by_name)
    if len(ranges) != 1:
        raise DatabaseParsingError(
            f"Solution parameter {parameter.parameter_type}({parameter.phase_name}, "
            f"{parameter.target[0] if parameter.target else ''};{parameter.order}) "
            "expands to multiple temperature ranges. The current RKMP runtime "
            "adapter supports one expression over the full range."
        )
    coefficients = _gibbs_coefficients(ranges[0].Gibbs)
    if np.count_nonzero(np.abs(coefficients[7:]) > 1e-12):
        raise DatabaseParsingError(
            f"Solution parameter {parameter.parameter_type}({parameter.phase_name}) "
            "contains custom Gibbs power terms that cannot be represented in the "
            "RKMP interaction slot."
        )
    return coefficients[1:7]


def _runtime_compounds(
    database: DatabaseIR,
    element_symbols: list[str],
    *,
    include_stoichiometric_compounds: bool,
    include_endmember_compounds: bool,
    endmember_compound_phases: Collection[str] | str | None,
    functions_by_name: dict[str, FunctionDefinition],
    phase_magnetic_options: dict[str, tuple[float, float]],
) -> list[_RuntimeCompound]:
    parameters_by_phase: dict[str, list[Parameter]] = {}
    for parameter in database.parameters:
        if parameter.order != 0:
            continue
        parameters_by_phase.setdefault(parameter.phase_name.upper(), []).append(
            parameter
        )

    compounds: list[_RuntimeCompound] = []
    if include_stoichiometric_compounds:
        for phase in database.phases:
            if phase.model.upper() != "COMPOUND":
                continue
            target = _compound_target(phase)
            parameter = _compound_gibbs_parameter(
                phase,
                target,
                _gibbs_parameters(parameters_by_phase.get(phase.name.upper(), [])),
            )
            ranges = _expanded_parameter_ranges(parameter, functions_by_name)
            stoichiometry = _compound_stoichiometry(phase, database, element_symbols)
            magnetic = _magnetic_data_for_target(
                phase,
                target,
                parameters_by_phase.get(phase.name.upper(), []),
                functions_by_name,
                phase_magnetic_options,
            )
            compounds.append(
                _RuntimeCompound(
                    name=phase.name,
                    stoichiometry=stoichiometry,
                    ranges=ranges,
                    phase_name=phase.name,
                    target=target,
                    magnetic=magnetic,
                )
            )

    selected_endmember_phases = _normalize_phase_filter(endmember_compound_phases)
    if include_endmember_compounds or selected_endmember_phases is not None:
        compounds.extend(
            _solution_endmember_compounds(
                database,
                element_symbols,
                functions_by_name,
                used_names={compound.name.upper() for compound in compounds},
                selected_phase_names=selected_endmember_phases,
                phase_magnetic_options=phase_magnetic_options,
            )
        )

    return compounds


def _solution_endmember_compounds(
    database: DatabaseIR,
    element_symbols: list[str],
    functions_by_name: dict[str, FunctionDefinition],
    used_names: set[str],
    selected_phase_names: set[str] | None,
    phase_magnetic_options: dict[str, tuple[float, float]],
) -> list[_RuntimeCompound]:
    phases_by_name = {phase.name.upper(): phase for phase in database.phases}
    compounds: list[_RuntimeCompound] = []
    for parameter in database.parameters:
        if parameter.parameter_type.upper() != "G" or parameter.order != 0:
            continue
        phase = phases_by_name.get(parameter.phase_name.upper())
        if phase is None or phase.model.upper() == "COMPOUND":
            continue
        if (
            selected_phase_names is not None
            and phase.name.upper() not in selected_phase_names
        ):
            continue
        if not parameter.target:
            continue
        target = parameter.target[0]
        if _target_is_not_endmember(target):
            continue
        stoichiometry = _endmember_stoichiometry(
            phase,
            target,
            database,
            element_symbols,
        )
        if not stoichiometry:
            continue
        name = _unique_runtime_name(
            f"{phase.name}_{_compact_target(target)}_ENDMBR",
            used_names,
        )
        used_names.add(name.upper())
        ranges = _expanded_parameter_ranges(parameter, functions_by_name)
        magnetic = _magnetic_data_for_target(
            phase,
            target,
            [
                candidate
                for candidate in database.parameters
                if candidate.phase_name.upper() == phase.name.upper()
                and candidate.order == 0
            ],
            functions_by_name,
            phase_magnetic_options,
        )
        compounds.append(
            _RuntimeCompound(
                name=name,
                stoichiometry=stoichiometry,
                ranges=ranges,
                phase_name=phase.name,
                target=target,
                magnetic=magnetic,
            )
        )
    return compounds


def _normalize_phase_filter(names: Collection[str] | str | None) -> set[str] | None:
    if names is None:
        return None
    if isinstance(names, str):
        names = [names]
    return {str(name).strip().upper() for name in names if str(name).strip()}


def _gibbs_parameters(parameters: list[Parameter]) -> list[Parameter]:
    return [
        parameter for parameter in parameters if parameter.parameter_type.upper() == "G"
    ]


def _phase_magnetic_options(database: DatabaseIR) -> dict[str, tuple[float, float]]:
    code_options: dict[str, tuple[float, float]] = {}
    direct_options: dict[str, tuple[float, float]] = {}
    disorder_parts = _phase_disordered_parts(database)

    for command in database.tdb_commands:
        if not command.active:
            continue
        parsed = command.parsed
        action = parsed.get("action", "").upper()
        if action == "MAGNETIC":
            option = _magnetic_option_from_command(command)
            if option is None:
                continue
            code = parsed.get("code", "").strip()
            target_phase = parsed.get("target_phase", "").strip()
            if code and target_phase == "@":
                code_options[code.upper()] = option
            elif target_phase:
                direct_options[target_phase.upper()] = option

    phases_by_name = {phase.name.upper(): phase for phase in database.phases}
    options: dict[str, tuple[float, float]] = {}
    for phase in database.phases:
        option = _phase_magnetic_option(
            phase,
            phases_by_name,
            code_options,
            direct_options,
            disorder_parts,
            seen=set(),
        )
        if option is not None:
            options[phase.name.upper()] = option
    return options


def _phase_disordered_parts(database: DatabaseIR) -> dict[str, str]:
    disorder_parts: dict[str, str] = {}
    phase_names = {phase.name.upper() for phase in database.phases}
    for command in database.tdb_commands:
        if not command.active:
            continue
        parsed = command.parsed
        if parsed.get("action", "").upper() != "DIS_PART":
            continue
        ordered = parsed.get("ordered_phase", "").strip()
        disordered = parsed.get("disordered_phase", "").strip()
        if ordered and disordered:
            disordered_key = disordered.upper()
            canonical = DISORDERED_PHASE_CANONICAL_NAMES.get(disordered_key)
            if disordered_key not in phase_names and canonical in phase_names:
                disordered_key = canonical
            disorder_parts[ordered.upper()] = disordered_key
    return disorder_parts


def _magnetic_option_from_command(
    command,
) -> tuple[float, float] | None:
    try:
        magnetic_model = float(command.parsed.get("magnetic_model", ""))
        structure_factor = abs(magnetic_model)
        if magnetic_model < 0.0 and structure_factor > 1.0:
            structure_factor = 1.0 / structure_factor
        return (
            structure_factor,
            float(command.parsed.get("magnetic_factor", "")),
        )
    except ValueError:
        return None


def _phase_magnetic_option(
    phase: Phase,
    phases_by_name: dict[str, Phase],
    code_options: dict[str, tuple[float, float]],
    direct_options: dict[str, tuple[float, float]],
    disorder_parts: dict[str, str],
    *,
    seen: set[str],
) -> tuple[float, float] | None:
    phase_key = phase.name.upper()
    if phase_key in seen:
        return None
    seen.add(phase_key)
    if phase_key in direct_options:
        return direct_options[phase_key]

    for code in _phase_type_codes(phase):
        if code.upper() in code_options:
            return code_options[code.upper()]

    disordered_name = disorder_parts.get(phase_key)
    if disordered_name:
        disordered_phase = phases_by_name.get(disordered_name)
        if disordered_phase is not None:
            return _phase_magnetic_option(
                disordered_phase,
                phases_by_name,
                code_options,
                direct_options,
                disorder_parts,
                seen=seen,
            )
    return None


def _phase_type_codes(phase: Phase) -> list[str]:
    match = re.search(
        r"(?i)(?:^|\s)PHASE\s+\S+\s+(\S+)",
        str(phase.source.command or ""),
    )
    token = match.group(1).strip() if match is not None else ""
    if not token.startswith("%"):
        return []
    return list(token[1:])


def _magnetic_data_for_target(
    phase: Phase,
    target: str,
    parameters: list[Parameter],
    functions_by_name: dict[str, FunctionDefinition],
    phase_magnetic_options: dict[str, tuple[float, float]],
) -> tuple[float, float, float, float] | None:
    option = phase_magnetic_options.get(phase.name.upper())
    if option is None:
        return None
    tc_parameter = _target_parameter_symmetric(parameters, "TC", target, phase)
    moment_parameter = (
        _target_parameter_symmetric(parameters, "BMAGN", target, phase)
        or _target_parameter_symmetric(parameters, "BMAG", target, phase)
        or _target_parameter_symmetric(parameters, "BM", target, phase)
    )
    if tc_parameter is None and moment_parameter is None:
        return (0.0, 0.0, option[0], option[1])
    tc_value = (
        _constant_parameter_value(tc_parameter, functions_by_name)
        if tc_parameter is not None
        else 0.0
    )
    moment_value = (
        _constant_parameter_value(moment_parameter, functions_by_name)
        if moment_parameter is not None
        else 0.0
    )
    return (
        _runtime_magnetic_scalar_value(tc_value, option),
        _runtime_magnetic_scalar_value(moment_value, option),
        option[0],
        option[1],
    )


def _runtime_magnetic_scalar_value(
    value: float,
    magnetic_option: tuple[float, float] | None,
) -> float:
    """Return the raw TDB magnetic scalar for runtime storage.

    TC-style TDB magnetic parameters are aggregated first.  If the aggregate
    Curie/Neel temperature or magnetic moment is negative, the magnetic model
    factor from ``TYPE_DEF ... MAGNETIC`` is then applied.  The existing
    Equilipy magnetic routine already performs the aggregate sign correction,
    so the parser must keep individual ``TC``/``BMAG`` values unscaled and
    encode any antiferromagnetic model factor in the phase magnetic option.
    """
    return value


def _target_parameter_symmetric(
    parameters: list[Parameter],
    parameter_type: str,
    target: str,
    phase: Phase | None,
) -> Parameter | None:
    exact = _target_parameter(parameters, parameter_type, target)
    if exact is not None:
        return exact
    target_sections = [_normalize_target(section) for section in target.split(":")]
    for variant in _symmetric_section_variants(
        target_sections,
        _runtime_equivalent_sublattice_groups(phase),
    ):
        variant_target = ":".join(variant)
        if variant_target == target:
            continue
        exact = _target_parameter(parameters, parameter_type, variant_target)
        if exact is not None:
            return exact
    return None


def _target_parameter(
    parameters: list[Parameter],
    parameter_type: str,
    target: str,
) -> Parameter | None:
    normalized_target = _normalize_target(target)
    normalized_type = parameter_type.upper()
    for parameter in parameters:
        if parameter.parameter_type.upper() != normalized_type:
            continue
        if parameter.order != 0 or not parameter.target:
            continue
        if _normalize_target(parameter.target[0]) == normalized_target:
            return parameter
    return None


def _magnetic_moment_priority(parameter_type: str) -> int:
    return {"BM": 0, "BMAG": 1, "BMAGN": 2}.get(parameter_type.upper(), -1)


def _constant_parameter_value(
    parameter: Parameter,
    functions_by_name: dict[str, FunctionDefinition],
) -> float:
    ranges = _expanded_parameter_ranges(parameter, functions_by_name)
    values = [_constant_gibbs_expression_value(item.Gibbs) for item in ranges]
    first = values[0]
    if any(abs(value - first) > 1e-9 for value in values[1:]):
        raise DatabaseParsingError(
            f"Magnetic parameter {parameter.parameter_type}({parameter.phase_name}) "
            "has temperature-range dependent values, which cannot be represented "
            "in the compound runtime magnetic slot."
        )
    return first


def _constant_gibbs_expression_value(expression: str) -> float:
    coefficients = _gibbs_coefficients(expression)
    if np.count_nonzero(np.abs(coefficients[2:]) > 1e-12):
        raise DatabaseParsingError(
            f"Magnetic scalar expression '{expression}' contains temperature "
            "terms and cannot be represented in the compound runtime magnetic slot."
        )
    return float(coefficients[1])


def _compound_target(phase: Phase) -> str:
    target_parts: list[str] = []
    for constituent_set in sorted(phase.constituents, key=lambda item: item.sublattice):
        if len(constituent_set.species) != 1:
            raise DatabaseParsingError(
                f"Phase '{phase.name}' is marked as a compound, but sublattice "
                f"{constituent_set.sublattice} has {len(constituent_set.species)} "
                "constituents."
            )
        target_parts.append(constituent_set.species[0])
    if not target_parts:
        raise DatabaseParsingError(
            f"Compound phase '{phase.name}' has no constituents."
        )
    return ":".join(target_parts)


def _compound_gibbs_parameter(
    phase: Phase,
    target: str,
    parameters: list[Parameter],
) -> Parameter:
    if not parameters:
        raise DatabaseParsingError(
            f"Compound phase '{phase.name}' is missing a G parameter."
        )
    normalized_target = _normalize_target(target)
    for parameter in parameters:
        if (
            parameter.target
            and _normalize_target(parameter.target[0]) == normalized_target
        ):
            return parameter
    if len(parameters) == 1:
        return parameters[0]
    targets = ", ".join(
        parameter.target[0] if parameter.target else "<empty>"
        for parameter in parameters
    )
    raise DatabaseParsingError(
        f"Compound phase '{phase.name}' has no G parameter for target "
        f"'{target}'. Available targets: {targets}."
    )


def _expanded_parameter_ranges(
    parameter: Parameter,
    functions_by_name: dict[str, FunctionDefinition],
) -> list[GibbsRange]:
    expression = parameter.expression.replace("#", "")
    source_ranges = _tdb._parse_gibbs_ranges(expression)
    if not source_ranges:
        raise DatabaseParsingError(
            f"Parameter {parameter.parameter_type}({parameter.phase_name}) has no "
            "parseable temperature range."
        )
    dependency_names = _tdb._referenced_function_names(expression, functions_by_name)
    dependency_bounds = _tdb._dependency_temperature_bounds(
        [(name, 1.0) for name in dependency_names],
        functions_by_name,
    )

    ranges: list[GibbsRange] = []
    expansion_cache: dict[tuple[str, float, float], str] = {}
    for source_range in source_ranges:
        bounds = _tdb._merged_bounds(
            [
                source_range.T_min,
                *[
                    bound
                    for bound in dependency_bounds
                    if source_range.T_min < bound < source_range.T_max
                ],
                source_range.T_max,
            ]
        )
        for lower, upper in zip(bounds, bounds[1:], strict=False):
            expanded = _tdb._expand_composite_expression(
                source_range.Gibbs,
                lower,
                upper,
                functions_by_name,
                seen=set(),
                cache=expansion_cache,
            )
            expanded = _tdb._simplify_gibbs_expression(expanded.replace("#", ""))
            ranges.append(GibbsRange(lower, upper, expanded, source_range.status))
    return ranges


def _compound_stoichiometry(
    phase: Phase,
    database: DatabaseIR,
    element_symbols: list[str],
) -> dict[str, float]:
    species_by_name = {species.name.upper(): species for species in database.species}
    element_by_upper = {symbol.upper(): symbol for symbol in element_symbols}
    stoichiometry = {symbol: 0.0 for symbol in element_symbols}
    for constituent_set in phase.constituents:
        site_ratio = float(constituent_set.site_ratio)
        for name in constituent_set.species:
            upper_name = name.upper()
            if upper_name == "VA":
                continue
            species = species_by_name.get(upper_name)
            if species is None and upper_name in element_by_upper:
                stoichiometry[element_by_upper[upper_name]] += site_ratio
                continue
            if species is None:
                raise DatabaseParsingError(
                    f"Compound phase '{phase.name}' references unknown species "
                    f"'{name}'."
                )
            for symbol, amount in species.composition.items():
                if symbol.upper() == "VA":
                    continue
                canonical = element_by_upper.get(symbol.upper())
                if canonical is None:
                    continue
                stoichiometry[canonical] += site_ratio * float(amount)
    return {symbol: amount for symbol, amount in stoichiometry.items() if amount}


def _endmember_stoichiometry(
    phase: Phase,
    target: str,
    database: DatabaseIR,
    element_symbols: list[str],
) -> dict[str, float]:
    species_by_name = {species.name.upper(): species for species in database.species}
    element_by_upper = {symbol.upper(): symbol for symbol in element_symbols}
    constituent_sets = sorted(phase.constituents, key=lambda item: item.sublattice)
    target_parts = [part.strip() for part in target.split(":")]
    if len(target_parts) != len(constituent_sets):
        raise DatabaseParsingError(
            f"Endmember target '{target}' for phase '{phase.name}' has "
            f"{len(target_parts)} sublattice part(s), but the phase has "
            f"{len(constituent_sets)} sublattice(s)."
        )

    stoichiometry = {symbol: 0.0 for symbol in element_symbols}
    for constituent_set, name in zip(constituent_sets, target_parts, strict=True):
        upper_name = name.upper()
        if upper_name == "VA":
            continue
        site_ratio = float(constituent_set.site_ratio)
        species = species_by_name.get(upper_name)
        if species is None and upper_name in element_by_upper:
            stoichiometry[element_by_upper[upper_name]] += site_ratio
            continue
        if species is None:
            raise DatabaseParsingError(
                f"Endmember target '{target}' for phase '{phase.name}' references "
                f"unknown species '{name}'."
            )
        for symbol, amount in species.composition.items():
            if symbol.upper() == "VA":
                continue
            canonical = element_by_upper.get(symbol.upper())
            if canonical is None:
                continue
            stoichiometry[canonical] += site_ratio * float(amount)
    return {symbol: amount for symbol, amount in stoichiometry.items() if amount}


def _target_is_not_endmember(target: str) -> bool:
    return any(token in target for token in {",", "*"})


def _compact_target(target: str) -> str:
    parts = []
    for part in target.split(":"):
        token = re.sub(r"[^A-Za-z0-9]+", "", part.strip())
        parts.append(token or "X")
    return "".join(parts)


def _unique_runtime_name(base: str, used_names: set[str]) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_]+", "_", base).strip("_") or "ENDMBR"
    candidate = sanitized[:_RUNTIME_NAME_WIDTH]
    if candidate.upper() not in used_names:
        return candidate
    digest = sha1(sanitized.encode("utf-8")).hexdigest()[:6].upper()
    stem = sanitized[: _RUNTIME_NAME_WIDTH - len(digest) - 1].rstrip("_")
    candidate = f"{stem}_{digest}"
    if candidate.upper() not in used_names:
        return candidate
    index = 2
    while True:
        suffix = f"_{index}"
        stem = sanitized[: _RUNTIME_NAME_WIDTH - len(digest) - len(suffix) - 1]
        candidate = f"{stem}_{digest}{suffix}"
        if candidate.upper() not in used_names:
            return candidate
        index += 1


def _gibbs_coefficients(expression: str) -> np.ndarray:
    coefficients = np.zeros(var.nGibbsCoeff, dtype=float)
    custom_terms: list[tuple[float, float]] = []
    for term in _tdb._split_additive_terms(_normalize_expression(expression)):
        parsed = _parse_gibbs_term(term)
        if parsed is None:
            raise DatabaseParsingError(
                f"Cannot convert Gibbs term '{term}' from expression '{expression}' "
                "to runtime coefficients."
            )
        coefficient, power = parsed
        if abs(power) < 1e-12:
            coefficients[1] += coefficient
        elif abs(power - 1.0) < 1e-12:
            coefficients[2] += coefficient
        elif abs(power - _T_LN_T_SENTINEL) < 1e-12:
            coefficients[3] += coefficient
        elif abs(power - _LN_T_SENTINEL) < 1e-12:
            _append_custom_gibbs_term(custom_terms, coefficient, power)
        elif abs(power - 2.0) < 1e-12:
            coefficients[4] += coefficient
        elif abs(power - 3.0) < 1e-12:
            coefficients[5] += coefficient
        elif abs(power + 1.0) < 1e-12:
            coefficients[6] += coefficient
        else:
            _append_custom_gibbs_term(custom_terms, coefficient, power)

    if len(custom_terms) > 3:
        raise DatabaseParsingError(
            "Runtime Gibbs storage supports at most three custom power terms; "
            f"found {len(custom_terms)} in '{expression}'."
        )
    for index, (coefficient, power) in enumerate(custom_terms):
        row = 7 + 2 * index
        coefficients[row] = coefficient
        coefficients[row + 1] = power
    return coefficients


def _append_custom_gibbs_term(
    custom_terms: list[tuple[float, float]],
    coefficient: float,
    power: float,
) -> None:
    # Keep any nonzero coefficient: tiny coefficients can be significant when
    # multiplied by high powers of T.
    for index, (existing_coefficient, existing_power) in enumerate(custom_terms):
        if abs(existing_power - power) < 1e-12:
            combined = existing_coefficient + coefficient
            if combined == 0.0:
                del custom_terms[index]
            else:
                custom_terms[index] = (combined, existing_power)
            return
    if coefficient != 0.0:
        custom_terms.append((coefficient, power))


def _parse_gibbs_term(term: str) -> tuple[float, float] | None:
    sign = -1.0 if term.startswith("-") else 1.0
    body = term[1:].strip() if term[:1] in "+-" else term.strip()
    body = body.replace(" ", "")
    if not body:
        return None
    constant = _tdb._constant_with_r_value(body)
    if constant is not None:
        return sign * constant, 0.0

    for pattern in (
        r"(?:(?P<coeff>.+)\*)?T\*LN\(T\)$",
        r"(?:(?P<coeff>.+)\*)?LN\(T\)\*T$",
    ):
        match = re.fullmatch(pattern, body, flags=re.IGNORECASE)
        if match is not None:
            coefficient = _coefficient(match.group("coeff"))
            return sign * coefficient, _T_LN_T_SENTINEL

    match = re.fullmatch(r"(?:(?P<coeff>.+)\*)?LN\(T\)$", body, flags=re.IGNORECASE)
    if match is not None:
        coefficient = _coefficient(match.group("coeff"))
        return sign * coefficient, _LN_T_SENTINEL

    match = re.fullmatch(
        rf"(?:(?P<coeff>.+)\*)?T\*\*\(?(?P<power>[+\-]?{_NUMBER_RE})\)?$",
        body,
        flags=re.IGNORECASE,
    )
    if match is not None:
        coefficient = _coefficient(match.group("coeff"))
        power = float(match.group("power"))
        return sign * coefficient, power

    if body == "T":
        return sign, 1.0
    if body.endswith("*T"):
        coefficient = _tdb._constant_with_r_value(body[:-2])
        if coefficient is not None:
            return sign * coefficient, 1.0
    if body.startswith("T*"):
        coefficient = _tdb._constant_with_r_value(body[2:])
        if coefficient is not None:
            return sign * coefficient, 1.0
    return None


def _coefficient(text: str | None) -> float:
    if text is None or not text.strip():
        return 1.0
    value = _tdb._constant_with_r_value(text)
    if value is None:
        raise DatabaseParsingError(f"Unsupported Gibbs coefficient '{text}'.")
    return value


def _normalize_expression(expression: str) -> str:
    text = expression.replace("#", "")
    text = re.sub(r"(?<=\d)[Dd](?=[+\-]?\d)", "E", text)
    text = re.sub(r"(?i)\bLOG\(", "LN(", text)
    return text


def _normalize_target(target: str) -> str:
    return target.replace(" ", "").upper()


def _fixed_width(value: object, width: int) -> str:
    return f"{str(value):<{width}}"[:width]
