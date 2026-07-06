"""Equilib result public model names and phase result records."""

from __future__ import annotations

from collections import defaultdict
from typing import Any, Dict, List, Optional, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from typing_extensions import TypedDict

import equilipy.equilifort as fort
import equilipy.variables as var
from equilipy.database_ir.tdb_canonical import DISORDERED_PHASE_CANONICAL_NAMES
from equilipy.exceptions import PostProcessError

from .capture import FortranCaptureState
from .common import (
    clean_quantity_mapping,
    component_names_from_amount_maps,
    normalize_phase_amount,
    unique_preserve_order,
)
from .context import ResultContext
from .phase import CompositionFractions, PhaseCollection, PhaseSpeciesProperty
from .serialization import (
    EQUILIB_BUNDLE_KINDS,
    context_from_payload,
    context_to_payload,
    single_point_from_payload,
    single_point_to_payload,
    validate_result_bundle,
)
from .table import ResultTable


class PhaseResult(BaseModel):
    """
    Pydantic model for a single phase, solution or compound.

    Holds the thermochemical state of one phase with lowercase public field
    names and clean composition accessors.
    """

    id: int
    name: str
    amount_n: float = Field(
        default=0.0,
        ge=0.0,
    )
    amount_w: float = Field(
        default=0.0,
        ge=0.0,
    )
    amount_n_basis: str = "phase_moles"
    stability: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
    )
    parent_model_id: Optional[int] = None
    parent_model_name: Optional[str] = None
    composition_set_id: int = Field(default=1, ge=0)
    display_label: Optional[str] = None
    ordering_degree: Optional[float] = None
    endmembers_x: Dict[str, float] = Field(default_factory=dict)
    endmembers_w: Dict[str, float] = Field(default_factory=dict)
    elements_x: Dict[str, float] = Field(default_factory=dict)
    elements_w: Dict[str, float] = Field(default_factory=dict)
    G: Optional[float] = None
    H: Optional[float] = None
    S: Optional[float] = None
    Cp: Optional[float] = None
    partial_gibbs: Dict[str, float] = Field(default_factory=dict)
    standard_gibbs_energy: Dict[str, float] = Field(default_factory=dict)
    activity: Dict[str, float] = Field(default_factory=dict)
    partial_enthalpy: Dict[str, float] = Field(default_factory=dict)
    partial_entropy: Dict[str, float] = Field(default_factory=dict)
    partial_heat_capacity: Dict[str, float] = Field(default_factory=dict)

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        extra="forbid",
        populate_by_name=True,
    )

    @field_validator("amount_n", "amount_w", mode="before")
    @classmethod
    def _normalize_tiny_phase_amount(cls, value: Any) -> float:
        """Clamp phase amount roundoff before non-negative validation."""
        return normalize_phase_amount(value)

    @property
    def endmembers(self) -> CompositionFractions:
        """Return endmember mole and mass fractions."""
        return CompositionFractions(
            x_i=dict(self.endmembers_x),
            w_i=dict(self.endmembers_w),
        )

    @property
    def elements(self) -> CompositionFractions:
        """Return element mole and mass fractions."""
        return CompositionFractions(
            x_i=dict(self.elements_x),
            w_i=dict(self.elements_w),
        )


class BatchPhaseResult:
    """Phase-like view that aggregates one phase across batch rows."""

    def __init__(self, name: str, phases: list[PhaseResult | None]) -> None:
        self.name = name
        self.id = _collect_phase_scalar(phases, "id", np.nan)
        self.amount_n = _collect_phase_scalar(phases, "amount_n", 0.0)
        self.amount_w = _collect_phase_scalar(phases, "amount_w", 0.0)
        self.amount_n_basis = _collect_phase_scalar(
            phases,
            "amount_n_basis",
            "",
        )
        self.stability = _collect_phase_scalar(phases, "stability", 0.0)
        self.parent_model_id = _collect_phase_scalar(phases, "parent_model_id", np.nan)
        self.parent_model_name = _collect_phase_scalar(phases, "parent_model_name", "")
        self.composition_set_id = _collect_phase_scalar(phases, "composition_set_id", 0)
        self.display_label = _collect_phase_scalar(phases, "display_label", "")
        self.ordering_degree = _collect_phase_scalar(phases, "ordering_degree", np.nan)
        self.endmembers_x = _collect_phase_mapping(phases, "endmembers_x")
        self.endmembers_w = _collect_phase_mapping(phases, "endmembers_w")
        self.elements_x = _collect_phase_mapping(phases, "elements_x")
        self.elements_w = _collect_phase_mapping(phases, "elements_w")
        self.G = _collect_phase_scalar(phases, "G", np.nan)
        self.H = _collect_phase_scalar(phases, "H", np.nan)
        self.S = _collect_phase_scalar(phases, "S", np.nan)
        self.Cp = _collect_phase_scalar(phases, "Cp", np.nan)
        self.partial_gibbs = _collect_phase_mapping(
            phases,
            "partial_gibbs",
        )
        self.standard_gibbs_energy = _collect_phase_mapping(
            phases,
            "standard_gibbs_energy",
        )
        self.activity = _collect_phase_mapping(phases, "activity")
        self.partial_enthalpy = _collect_phase_mapping(phases, "partial_enthalpy")
        self.partial_entropy = _collect_phase_mapping(phases, "partial_entropy")
        self.partial_heat_capacity = _collect_phase_mapping(
            phases,
            "partial_heat_capacity",
        )

    @property
    def endmembers(self) -> CompositionFractions:
        """Return endmember mole and mass fractions across batch rows."""
        return CompositionFractions(
            x_i=dict(self.endmembers_x),
            w_i=dict(self.endmembers_w),
        )

    @property
    def elements(self) -> CompositionFractions:
        """Return element mole and mass fractions across batch rows."""
        return CompositionFractions(
            x_i=dict(self.elements_x),
            w_i=dict(self.elements_w),
        )


def _collect_phase_scalar(
    phases: list[PhaseResult | None],
    attribute_name: str,
    missing_value: float,
) -> list[Any]:
    """Collect one scalar phase attribute across batch rows."""
    return [
        getattr(phase, attribute_name, missing_value)
        if phase is not None
        else missing_value
        for phase in phases
    ]


def _collect_phase_mapping(
    phases: list[PhaseResult | None],
    attribute_name: str,
) -> dict[str, list[Any]]:
    """Collect one species mapping across batch rows."""
    keys: list[str] = []
    seen: set[str] = set()
    for phase in phases:
        if phase is None:
            continue
        values = getattr(phase, attribute_name, {})
        for key in values:
            if key not in seen:
                seen.add(key)
                keys.append(key)

    return {
        key: [
            getattr(phase, attribute_name, {}).get(key, np.nan)
            if phase is not None
            else np.nan
            for phase in phases
        ]
        for key in keys
    }


def _batch_phase_result_map(
    points: list["EquilibPoint"],
) -> dict[str, BatchPhaseResult]:
    """Return phase-like batch views keyed by phase name."""
    phase_names: list[str] = []
    for point in points:
        phase_names.extend(point.phase_map)
    phase_names = unique_preserve_order(phase_names)

    return {
        phase_name: BatchPhaseResult(
            phase_name,
            [point.phase_map.get(phase_name) for point in points],
        )
        for phase_name in phase_names
    }


def count_unstable_compounds(i):
    """Return the number of skipped unstable compounds before species index i."""
    n = 0
    iphase = fort.modulethermo.iphase
    iLastSoln = fort.modulethermo.nspeciesphase[len(var.iSys2DBSoln)]
    for k in range(iLastSoln, i):
        if iphase[k] < 0:
            n = n + 1
    return n


def get_assemblage_name(assemblage_ids):
    """
    Return phase names based on assemblage IDs.

    Variables
    =========

    Input
    assemblage_ids: list of phase assemblage id defined after selecting phases
        based on input elements

    Output
    AssemblageNames
    """
    AssemblageNames = list([])
    for i in assemblage_ids:
        if i == 0:
            # Empty phase: place holder
            AssemblageNames.append("{:<1}".format(""))
        elif i < 0:
            # Solution PhaseResult
            # Note that the solution phase id changes with phase selection
            AssemblageNames.append(
                _public_phase_name(var.cPhaseNameSys[int(-(i + 1))])
            )
        else:
            # Compound PhaseResult

            # Note that the compund id doesn't change with phase selection
            # Correction is made by counting the number of ignored compound phases
            nSoln = len(var.iSys2DBSoln)
            iLastSoln = fort.modulethermo.nspeciesphase[nSoln]
            id = -iLastSoln + i + nSoln - 1
            n = count_unstable_compounds(i)
            id = int(id - n)
            AssemblageNames.append(_public_phase_name(var.cPhaseNameSys[id]))
    return AssemblageNames


def _public_phase_name(phase_name: str) -> str:
    """Return the public name for a runtime phase alias."""
    name = str(phase_name).strip()
    return DISORDERED_PHASE_CANONICAL_NAMES.get(name.upper(), name)


def _prefer_phase_result(candidate: PhaseResult, current: PhaseResult | None) -> bool:
    """Return whether candidate should represent a public phase alias."""
    if current is None:
        return True
    if candidate.stability > current.stability:
        return True
    if (
        candidate.stability == current.stability
        and candidate.amount_n > current.amount_n
    ):
        return True
    return False


def _stable_phase_id_by_name(
    capture: FortranCaptureState,
    phase_name: str,
) -> Optional[int]:
    """Return the Fortran assemblage ID for a stable phase name."""
    for phase_id, stable_name in zip(
        capture.assemblage_ids,
        get_assemblage_name(capture.assemblage_ids),
        strict=False,
    ):
        if stable_name.strip() == phase_name.strip():
            return int(phase_id)
    return None


def _phase_property_dict(
    species_names: list[str],
    values: np.ndarray,
) -> Dict[str, float]:
    """Return a clean species-property mapping for one phase."""
    return {
        name: float(value)
        for name, value in zip(species_names, values, strict=False)
    }


def _phase_weighted_property(
    fractions: list[float],
    values: np.ndarray,
) -> float:
    """Return a phase molar property from species fractions and partial values."""
    if len(fractions) == 0:
        return np.nan
    fraction_array = np.asarray(fractions, dtype=float)
    value_array = np.asarray(values, dtype=float)
    if len(fraction_array) != len(value_array):
        return np.nan
    return float(np.dot(fraction_array, value_array))


def _species_mass_fractions(
    species_indices: np.ndarray,
    species_fractions: list[float] | np.ndarray,
    elements: list[str],
) -> list[float]:
    """Return mass fractions for a species mixture on the active element basis."""
    fractions = np.asarray(species_fractions, dtype=float)
    if len(species_indices) != len(fractions) or len(fractions) == 0:
        return list(fractions)

    try:
        stoich = np.asarray(var.dStoichSpeciesCS, dtype=float)
        active_elements = np.asarray(var.iElementDBIndex, dtype=int)
        masses = np.asarray(
            [var.cPeriodicTable[element.strip()][1] for element in elements],
            dtype=float,
        )
    except (AttributeError, KeyError, TypeError, ValueError):
        return list(fractions)

    species_indices = np.asarray(species_indices, dtype=int)
    if np.any(species_indices < 0) or np.any(species_indices >= len(stoich)):
        return list(fractions)
    if active_elements.size != len(masses) or np.any(
        active_elements >= stoich.shape[1]
    ):
        return list(fractions)

    species_masses = stoich[species_indices][:, active_elements] @ masses
    weights = fractions * species_masses
    total = float(np.sum(weights))
    if total <= 0.0 or not np.isfinite(total):
        return list(fractions)
    return list(weights / total)


def _element_mass_fractions(
    elements_x: dict[str, float],
) -> dict[str, float]:
    """Convert element mole fractions to mass fractions."""
    weights: dict[str, float] = {}
    total = 0.0
    for element, fraction in elements_x.items():
        try:
            mass = float(var.cPeriodicTable[element.strip()][1])
        except (KeyError, TypeError, ValueError):
            return dict(elements_x)
        weight = max(float(fraction), 0.0) * mass
        weights[element] = weight
        total += weight
    if total <= 0.0 or not np.isfinite(total):
        return dict(elements_x)
    return {element: weight / total for element, weight in weights.items()}


def _phase_name_by_system_id(capture: FortranCaptureState, phase_id: int) -> str:
    """Return the public phase name for a one-based system phase id."""
    if phase_id <= 0 or phase_id > len(capture.phase_names):
        return ""
    return _public_phase_name(capture.phase_names[phase_id - 1]).strip()


def _solution_model_type(system_phase_id: int) -> str:
    """Return the solution model type for a one-based system solution phase id."""
    try:
        database_phase_id = int(var.iSys2DBSoln[system_phase_id - 1])
        return str(var.cSolnPhaseTypeCS[database_phase_id - 1]).strip().upper()
    except (AttributeError, IndexError, TypeError, ValueError):
        return ""


def _solution_sublattice_index(system_phase_id: int) -> int:
    """Return the zero-based sublattice-model index for a system phase."""
    try:
        database_phase_id = int(var.iSys2DBSoln[system_phase_id - 1])
        return int(var.iPhaseSublatticeCS[database_phase_id - 1]) - 1
    except (AttributeError, IndexError, TypeError, ValueError):
        return -1


def _active_slot_value(values: np.ndarray, slot_index: int, default: int) -> int:
    """Return one active-slot integer value with a default fallback."""
    try:
        value = int(values[slot_index])
    except (IndexError, TypeError, ValueError):
        return default
    return value if value > 0 else default


def _active_slot_site_fraction(
    capture: FortranCaptureState,
    slot_index: int | None,
) -> np.ndarray | None:
    """Return one active-slot site-fraction snapshot when available."""
    if slot_index is None:
        return None
    site = capture.active_slot_site_fraction
    if site.ndim != 3 or slot_index < 0 or slot_index >= site.shape[0]:
        return None
    slot_site = np.asarray(site[slot_index], dtype=float)
    if slot_site.size == 0 or float(np.sum(slot_site)) <= 0.0:
        return None
    return slot_site


def _solution_species_slice(
    capture: FortranCaptureState,
    system_phase_id: int,
) -> tuple[slice, np.ndarray] | None:
    """Return selected-system species slice and database species ids."""
    phase_index = system_phase_id - 1
    if phase_index < 0 or phase_index >= capture.solution_count:
        return None
    first = int(capture.species_phase_boundaries[phase_index])
    last = int(capture.species_phase_boundaries[phase_index + 1])
    if last <= first:
        return None
    species_slice = slice(first, last)
    return species_slice, capture.system_to_database_species[species_slice]


def _site_product_fractions(
    system_phase_id: int,
    site_fraction: np.ndarray,
    species_count: int,
) -> np.ndarray | None:
    """Return endmember product fractions from CEF site fractions."""
    sublattice_index = _solution_sublattice_index(system_phase_id)
    if sublattice_index < 0:
        return None

    try:
        n_sublattices = int(var.nSublatticePhaseCS[sublattice_index])
        constituent = np.asarray(var.iConstituentSublatticeCS, dtype=int)
    except (AttributeError, IndexError, TypeError, ValueError):
        return None

    if species_count <= 0 or site_fraction.ndim != 2:
        return None
    if sublattice_index >= constituent.shape[0]:
        return None

    fractions = np.zeros(species_count, dtype=float)
    for local_species in range(species_count):
        product = 1.0
        for sublattice in range(min(n_sublattices, site_fraction.shape[0])):
            if local_species >= constituent.shape[2]:
                product = 0.0
                break
            constituent_index = int(
                constituent[sublattice_index, sublattice, local_species]
            ) - 1
            if constituent_index < 0:
                continue
            if constituent_index >= site_fraction.shape[1]:
                product = 0.0
                break
            product *= max(float(site_fraction[sublattice, constituent_index]), 0.0)
        fractions[local_species] = product

    total = float(np.sum(fractions))
    if total <= 0.0 or not np.isfinite(total):
        return None
    return fractions / total


def _species_fractions_from_element_composition(
    species_indices: np.ndarray,
    elements_x: dict[str, float],
    elements: list[str],
) -> np.ndarray | None:
    """Project an element composition onto one phase's endmember basis."""
    if len(species_indices) == 0:
        return None

    try:
        active_elements = np.asarray(var.iElementDBIndex, dtype=int)
        stoich = np.asarray(var.dStoichSpeciesCS, dtype=float)
    except (AttributeError, TypeError, ValueError):
        return None

    species_indices = np.asarray(species_indices, dtype=int)
    if np.any(species_indices < 0) or np.any(species_indices >= len(stoich)):
        return None
    if active_elements.size != len(elements) or np.any(
        active_elements >= stoich.shape[1]
    ):
        return None

    target = np.asarray(
        [float(elements_x.get(element, 0.0)) for element in elements],
        dtype=float,
    )
    target_total = float(np.sum(target))
    if target_total <= 0.0 or not np.isfinite(target_total):
        return None
    target = target / target_total

    species_stoich = stoich[species_indices][:, active_elements]
    matrix = np.vstack([species_stoich.T, np.ones(len(species_indices))])
    rhs = np.concatenate([target, np.ones(1)])
    try:
        fractions = np.linalg.lstsq(matrix, rhs, rcond=None)[0]
    except np.linalg.LinAlgError:
        return None

    fractions[np.abs(fractions) < 1e-12] = 0.0
    if np.any(fractions < -1e-8):
        return None
    fractions = np.maximum(fractions, 0.0)
    total = float(np.sum(fractions))
    if total <= 0.0 or not np.isfinite(total):
        return None
    fractions = fractions / total

    reconstructed = species_stoich.T @ fractions
    reconstructed_total = float(np.sum(reconstructed))
    if reconstructed_total > 0.0:
        reconstructed = reconstructed / reconstructed_total
    if float(np.linalg.norm(reconstructed - target)) > 1e-8:
        return None
    return fractions


def _active_slot_solution_report_fractions(
    capture: FortranCaptureState,
    system_phase_id: int,
    slot_index: int | None,
    display_species_indices: np.ndarray,
    elements: list[str],
) -> tuple[np.ndarray | None, dict[str, float] | None]:
    """Return slot fractions and element composition for a solution row."""
    site_fraction = _active_slot_site_fraction(capture, slot_index)
    if site_fraction is None:
        return None, None

    parent_phase_id = _active_slot_value(
        capture.active_slot_thermo_phase,
        slot_index if slot_index is not None else -1,
        system_phase_id,
    )
    composition_set_id = _active_slot_value(
        capture.active_slot_identity_ordinal,
        slot_index if slot_index is not None else -1,
        1,
    )
    if parent_phase_id == system_phase_id and composition_set_id <= 1:
        return None, None

    if _solution_model_type(parent_phase_id) not in {"SUBL", "SUBLM", "SUBOM", "SUBM"}:
        return None, None

    parent_slice_info = _solution_species_slice(capture, parent_phase_id)
    if parent_slice_info is None:
        return None, None
    _, parent_species_indices = parent_slice_info
    parent_fractions = _site_product_fractions(
        parent_phase_id,
        site_fraction,
        len(parent_species_indices),
    )
    if parent_fractions is None:
        return None, None

    elements_x = species_indices2elements(
        parent_species_indices, parent_fractions, elements
    )
    if system_phase_id == parent_phase_id:
        return parent_fractions, elements_x

    display_fractions = _species_fractions_from_element_composition(
        display_species_indices,
        elements_x,
        elements,
    )
    return display_fractions, elements_x


def _ordering_degree_from_site(
    system_phase_id: int, site_fraction: np.ndarray | None
) -> float:
    """Return a SUBOM ordering degree from equivalent sublattice differences."""
    if site_fraction is None or _solution_model_type(system_phase_id) != "SUBOM":
        return 0.0

    sublattice_index = _solution_sublattice_index(system_phase_id)
    if sublattice_index < 0:
        return 0.0

    try:
        n_sublattices = int(var.nSublatticePhaseCS[sublattice_index])
        n_constituent = np.asarray(var.nConstituentSublatticeCS, dtype=int)
        constituent = np.asarray(var.iConstituentSublatticeCS, dtype=int)
        stoich = np.asarray(var.dStoichSublatticeCS, dtype=float)
    except (AttributeError, IndexError, TypeError, ValueError):
        return 0.0

    groups: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for sublattice in range(min(n_sublattices, site_fraction.shape[0])):
        count = int(n_constituent[sublattice_index, sublattice])
        if count <= 1:
            continue
        species_signature = tuple(
            int(value)
            for value in constituent[sublattice_index, sublattice, :count]
        )
        stoich_signature = round(float(stoich[sublattice_index, sublattice]), 12)
        groups[(stoich_signature, species_signature)].append(sublattice)

    degree = 0.0
    for sublattice_group in groups.values():
        if len(sublattice_group) < 2:
            continue
        count = int(n_constituent[sublattice_index, sublattice_group[0]])
        for i, sublattice_a in enumerate(sublattice_group[:-1]):
            values_a = site_fraction[sublattice_a, :count]
            for sublattice_b in sublattice_group[i + 1 :]:
                values_b = site_fraction[sublattice_b, :count]
                degree = max(degree, float(np.max(np.abs(values_a - values_b))))
    return degree


def _disordered_helper_name(
    capture: FortranCaptureState,
    system_phase_id: int,
) -> str:
    """Return the public disordered helper name for an ordered phase, if known."""
    try:
        disordered_phase = np.asarray(fort.modulethermo.idisorderedphase, dtype=int)
        helper_id = int(disordered_phase[system_phase_id - 1])
    except (AttributeError, IndexError, TypeError, ValueError):
        helper_id = 0
    helper_name = _phase_name_by_system_id(capture, helper_id)
    return helper_name or _phase_name_by_system_id(capture, system_phase_id)


def _solution_display_label(
    capture: FortranCaptureState,
    parent_phase_id: int,
    display_phase_id: int,
    ordering_degree: float,
) -> str:
    """Return the nonbreaking reporting label for one solution phase slot."""
    base_name = _phase_name_by_system_id(capture, display_phase_id)
    parent_name = _phase_name_by_system_id(capture, parent_phase_id) or base_name
    if _solution_model_type(parent_phase_id) != "SUBOM":
        return base_name
    if ordering_degree <= 1e-6:
        return f"{_disordered_helper_name(capture, parent_phase_id)}-like"
    return f"{parent_name}-like"


def _solution_phase_metadata(
    capture: FortranCaptureState,
    system_phase_id: int,
    slot_index: int | None = None,
) -> dict[str, Any]:
    """Return additive reporting metadata for one solution phase or slot."""
    parent_phase_id = system_phase_id
    display_phase_id = system_phase_id
    composition_set_id = 1
    if slot_index is not None:
        parent_phase_id = _active_slot_value(
            capture.active_slot_thermo_phase,
            slot_index,
            system_phase_id,
        )
        display_phase_id = _active_slot_value(
            capture.active_slot_display_phase,
            slot_index,
            system_phase_id,
        )
        composition_set_id = _active_slot_value(
            capture.active_slot_identity_ordinal,
            slot_index,
            1,
        )

    site_fraction = _active_slot_site_fraction(capture, slot_index)
    ordering_degree = _ordering_degree_from_site(parent_phase_id, site_fraction)
    return {
        "parent_model_id": int(parent_phase_id),
        "parent_model_name": _phase_name_by_system_id(capture, parent_phase_id),
        "composition_set_id": int(composition_set_id),
        "display_label": _solution_display_label(
            capture,
            parent_phase_id,
            display_phase_id,
            ordering_degree,
        ),
        "ordering_degree": float(ordering_degree),
    }


def _stable_phase_metadata(
    capture: FortranCaptureState,
) -> dict[str, np.ndarray]:
    """Return stable-phase reporting metadata aligned to assemblage rows."""
    parent_model_id: list[int] = []
    parent_model_name: list[str] = []
    composition_set_id: list[int] = []
    display_label: list[str] = []
    ordering_degree: list[float] = []

    for slot_index, phase_id_raw in enumerate(capture.assemblage_ids):
        phase_id = int(phase_id_raw)
        if phase_id < 0:
            metadata = _solution_phase_metadata(
                capture,
                system_phase_id=-phase_id,
                slot_index=slot_index,
            )
            parent_model_id.append(metadata["parent_model_id"])
            parent_model_name.append(metadata["parent_model_name"])
            composition_set_id.append(metadata["composition_set_id"])
            display_label.append(metadata["display_label"])
            ordering_degree.append(metadata["ordering_degree"])
        elif phase_id > 0:
            name = get_assemblage_name([phase_id])[0].strip()
            parent_model_id.append(phase_id)
            parent_model_name.append(name)
            composition_set_id.append(1)
            display_label.append(name)
            ordering_degree.append(np.nan)
        else:
            parent_model_id.append(0)
            parent_model_name.append("")
            composition_set_id.append(0)
            display_label.append("")
            ordering_degree.append(np.nan)

    return {
        "parent_model_id": np.asarray(parent_model_id, dtype=int),
        "parent_model_name": np.asarray(parent_model_name),
        "composition_set_id": np.asarray(composition_set_id, dtype=int),
        "display_label": np.asarray(display_label),
        "ordering_degree": np.asarray(ordering_degree, dtype=float),
    }


def _stable_phase_amounts_n(
    capture: FortranCaptureState,
    stable_names: np.ndarray,
    all_phases: dict[str, PhaseResult],
) -> np.ndarray:
    """Return stable phase amounts without collapsing duplicate labels."""
    stripped_names = [str(name).strip() for name in stable_names]
    duplicate_names = {
        name for name in stripped_names if name and stripped_names.count(name) > 1
    }
    amounts: list[float] = []
    for name, raw_amount in zip(stripped_names, capture.phase_amounts_n, strict=False):
        if name in duplicate_names:
            amounts.append(float(raw_amount))
            continue
        phase = all_phases.get(name)
        amounts.append(
            float(phase.amount_n) if phase is not None else float(raw_amount)
        )
    return np.asarray(amounts, dtype=float)


def _stable_phase_amount_basis(
    stable_names: np.ndarray,
    all_phases: dict[str, PhaseResult],
) -> np.ndarray:
    """Return stable phase amount bases aligned with assemblage rows."""
    bases: list[str] = []
    for name in [str(value).strip() for value in stable_names]:
        phase = all_phases.get(name)
        bases.append(phase.amount_n_basis if phase is not None else "unknown")
    return np.asarray(bases)


def _phase_occupied_atom_count(
    first: int,
    last: int,
    fractions: list[float],
) -> float:
    """Return the active occupied-atom denominator for a solution phase."""
    fraction_array = np.asarray(fractions, dtype=float)
    total_fraction = float(np.sum(fraction_array))
    if total_fraction <= 0.0:
        return 1.0

    stoich = np.asarray(fort.modulethermo.dstoichspecies, dtype=float).copy()
    nelements = int(fort.modulethermo.nelements)
    occupied = (fraction_array / total_fraction) @ stoich[first:last, :nelements]
    atom_count = float(np.sum(occupied))
    if atom_count <= 0.0 or not np.isfinite(atom_count):
        return 1.0
    return atom_count


def _solution_reported_phase_amounts(
    raw_amount: float,
    raw_weight: float,
    database_species_indices: np.ndarray,
    fractions: list[float],
    solution_type: str,
) -> tuple[float, float, str]:
    """Return solution phase amount on the active pseudo-component basis.

    MQM/SUBQ phases use internal pair or quadruplet moles.  When the user
    supplied a rank-reduced pseudo-component basis such as CaO-Al2O3-SiO2,
    public solution phase amounts should be reported as pseudo-component
    moles when the phase stoichiometry projects cleanly onto that basis.
    """
    amount_basis = "solution_native"
    if raw_amount <= 0.0:
        return raw_amount, raw_weight, amount_basis
    component_stoich = np.asarray(
        getattr(var, "dPseudoComponentStoichSys", np.empty((0, 0))),
        dtype=float,
    )
    if component_stoich.ndim != 2 or 0 in component_stoich.shape:
        return raw_amount, raw_weight, amount_basis
    if len(database_species_indices) != len(fractions):
        return raw_amount, raw_weight, amount_basis

    full_stoich = np.asarray(var.dStoichSpeciesCS, dtype=float)
    species_indices = np.asarray(database_species_indices, dtype=int)
    if np.any(species_indices < 0) or np.any(species_indices >= len(full_stoich)):
        return raw_amount, raw_weight, amount_basis

    particles = np.asarray(
        getattr(var, "iParticlesPerMoleCS", np.ones(len(full_stoich))),
        dtype=float,
    )
    if np.any(species_indices >= len(particles)):
        return raw_amount, raw_weight, amount_basis
    species_particles = particles[species_indices]
    species_particles[species_particles <= 0.0] = 1.0

    component_matrix = component_stoich.T
    if component_matrix.shape[0] != full_stoich.shape[1]:
        return raw_amount, raw_weight, amount_basis

    fraction_array = np.asarray(fractions, dtype=float)
    if np.any(fraction_array < -1e-10):
        return raw_amount, raw_weight, amount_basis
    fraction_sum = float(np.sum(fraction_array))
    if fraction_sum <= 0.0:
        return raw_amount, raw_weight, amount_basis
    fraction_array = fraction_array / fraction_sum

    species_basis_stoich = full_stoich[species_indices, :].copy()
    if solution_type.strip() not in {"SUBG", "SUBQ"}:
        species_basis_stoich = species_basis_stoich / species_particles[:, None]

    species_component_coeff = np.zeros(
        (len(species_indices), component_stoich.shape[0])
    )
    for i_species, stoich_vector in enumerate(species_basis_stoich):
        coeff = np.linalg.lstsq(component_matrix, stoich_vector, rcond=None)[0]
        reconstructed = component_matrix @ coeff
        scale = max(1.0, float(np.linalg.norm(stoich_vector)))
        residual = float(np.linalg.norm(stoich_vector - reconstructed))
        if residual > 1e-8 * scale or np.any(coeff < -1e-8):
            return raw_amount, raw_weight, amount_basis
        coeff[np.abs(coeff) < 1e-12] = 0.0
        species_component_coeff[i_species, :] = coeff

    component_per_native_mole = fraction_array @ species_component_coeff
    if np.any(component_per_native_mole < -1e-10):
        return raw_amount, raw_weight, amount_basis
    component_per_native_mole[np.abs(component_per_native_mole) < 1e-12] = 0.0

    if solution_type.strip() not in {"SUBG", "SUBQ"}:
        projected = _solution_amount_from_projected_formula_mass(
            raw_amount,
            component_per_native_mole,
            component_stoich,
        )
        if projected is not None:
            projected_amount, projected_weight = projected
            return projected_amount, projected_weight, "pseudo_formula_moles"

    coefficients = raw_amount * component_per_native_mole
    reconstructed = component_matrix @ coefficients
    element_amounts = raw_amount * (fraction_array @ species_basis_stoich)
    scale = max(1.0, float(np.linalg.norm(element_amounts)))
    residual = float(np.linalg.norm(element_amounts - reconstructed))
    if residual > 1e-8 * scale or np.any(coefficients < -1e-8):
        return raw_amount, raw_weight, amount_basis

    coefficients[np.abs(coefficients) < 1e-12] = 0.0
    atomic_masses = np.asarray(getattr(var, "dAtomicMass", []), dtype=float)
    reported_weight = raw_weight
    if len(atomic_masses) == len(element_amounts):
        reported_weight = float(np.dot(element_amounts, atomic_masses))
    return float(np.sum(coefficients)), reported_weight, "pseudo_component_moles"


def _solution_amount_from_projected_formula_mass(
    raw_amount: float,
    component_per_native_mole: np.ndarray,
    component_stoich: np.ndarray,
) -> tuple[float, float] | None:
    """Return formula-mole amount from projected component composition and mass."""
    if raw_amount <= 0.0:
        return None
    component_total = float(np.sum(component_per_native_mole))
    if component_total <= 0.0 or not np.isfinite(component_total):
        return None

    atomic_masses = np.asarray(getattr(var, "dAtomicMass", []), dtype=float)
    if atomic_masses.ndim != 1 or component_stoich.shape[1] != len(atomic_masses):
        return None
    component_masses = component_stoich @ atomic_masses
    if np.any(component_masses <= 0.0):
        return None

    component_fraction = component_per_native_mole / component_total
    formula_molar_mass = float(np.dot(component_fraction, component_masses))
    if formula_molar_mass <= 0.0 or not np.isfinite(formula_molar_mass):
        return None
    formula_amount = float(raw_amount * component_total)
    return formula_amount, float(formula_amount * formula_molar_mass)


def _has_pseudo_component_context() -> bool:
    """Return whether the current input uses a pseudo-component basis."""
    component_stoich = np.asarray(
        getattr(var, "dPseudoComponentStoichSys", np.empty((0, 0))),
        dtype=float,
    )
    return component_stoich.ndim == 2 and all(
        size > 0 for size in component_stoich.shape
    )


def _compound_reported_phase_amount_n(
    raw_amount: float,
    database_species_index: int,
) -> tuple[float, str]:
    """Return public compound amount and its reporting basis.

    In pseudo-component oxide systems, public stoichiometric compound amounts
    follow FactSage workbook semantics: formula moles.  The legacy active-atom
    scaling is retained outside pseudo-component systems to avoid changing
    unrelated result contracts in this slice.
    """
    if raw_amount <= 0.0:
        return raw_amount, "formula_moles"
    if _has_pseudo_component_context():
        return raw_amount, "formula_moles"
    try:
        stoich = np.asarray(var.dStoichSpeciesCS, dtype=float)
        active_indices = np.asarray(var.iElementDBIndex, dtype=int)
        particles = np.asarray(var.iParticlesPerMoleCS, dtype=float)
    except (AttributeError, TypeError, ValueError):
        return raw_amount, "phase_moles"

    if database_species_index < 0 or database_species_index >= stoich.shape[0]:
        return raw_amount, "phase_moles"
    if active_indices.size == 0 or np.any(active_indices >= stoich.shape[1]):
        return raw_amount, "phase_moles"

    particle_count = 1.0
    if (
        database_species_index < len(particles)
        and particles[database_species_index] > 0.0
    ):
        particle_count = float(particles[database_species_index])

    atom_count = (
        float(np.sum(np.abs(stoich[database_species_index, active_indices])))
        / particle_count
    )
    if atom_count <= 0.0 or not np.isfinite(atom_count):
        return raw_amount, "phase_moles"
    return raw_amount * atom_count, "active_atom_moles"


def _phase_amount_from_mass_and_composition(
    phase_mass: float,
    elements_x: dict[str, float],
) -> float | None:
    """Return phase mole amount from mass and active element mole fractions."""
    if phase_mass <= 0.0 or not elements_x:
        return None

    molar_mass = 0.0
    for element_name, fraction in elements_x.items():
        try:
            atomic_mass = float(var.cPeriodicTable[element_name.strip()][1])
        except (KeyError, TypeError, ValueError):
            return None
        if not np.isfinite(fraction) or fraction < -1e-12:
            return None
        molar_mass += max(float(fraction), 0.0) * atomic_mass

    if molar_mass <= 0.0 or not np.isfinite(molar_mass):
        return None
    return phase_mass / molar_mass


def safe_error_component_values(values, indices, size: int) -> np.ndarray:
    """Return indexed component values or NaNs when Fortran output is unavailable."""
    if size <= 0:
        return np.array([])
    try:
        if values is None:
            raise TypeError("Fortran component values are unavailable.")
        indexed_values = np.asarray(values).copy()[indices]
        if len(indexed_values) == size:
            return indexed_values
    except (IndexError, TypeError, ValueError):
        pass
    return np.full(size, np.nan)


def endmembers2elements(
    phase_name: str,
    endmembers: Dict[str, float],
    elements: list,
    unit_out: Optional[List[str]] = None,
    is_solution: bool = True,
) -> Dict[str, float]:
    """Convert endmember amounts to element amounts."""
    if unit_out is None:
        unit_out = ["K", "atm", "moles"]

    # 1. If the elements in endmembers and input condition are same, return
    EndmembersName = list(endmembers.keys())

    if set(EndmembersName) == set(elements):
        return endmembers

    endmemeber_amounts = np.array([endmembers.get(name) for name in EndmembersName])

    is_mass_based = unit_out[2] in [
        "grams",
        "kilograms",
        "pounds",
        "g",
        "kg",
        "lbs",
        "mass fraction",
        "weight fraction",
        "wt%",
        "wt.%",
    ]
    if is_mass_based:
        atomic_masses = np.array([var.cPeriodicTable[el.strip()][1] for el in elements])

    # 2. If not, calulate element fraction from database

    iElementDBIndex = var.iElementDBIndex.copy()
    dStoichSpeciesPhase = var.dStoichSpeciesCS[:, iElementDBIndex]

    if is_solution:
        # 2.1 Get the phase index
        DB_soln_names = var.cSolnPhaseNameCS.copy()
        DB_soln_names = [x.strip() for x in DB_soln_names]

        counts = defaultdict(int)
        all_soln_names = []
        for name in DB_soln_names:
            counts[name] += 1
            count = counts[name]
            if count > 1:
                all_soln_names.append(f"{name}#{count}")
            else:
                all_soln_names.append(name)

        # i_system=var.cSolnPhaseNameCS.index("{:<25}".format(phase_name))

        i_system = all_soln_names.index(phase_name)
        nSpeciesIndex = var.nSpeciesPhaseCS.copy()
        nSpeciesIndex = np.append([0], nSpeciesIndex)
        iFirstSys = nSpeciesIndex[i_system]
        iLastSys = nSpeciesIndex[i_system + 1]
        Temp = [
            str(x).strip() for x in np.array(var.cEndmemberNameCS)[iFirstSys:iLastSys]
        ]

        # iSpeciesDBIndex = np.where(
        #     np.isin(var.iSys2DBSpeciesDefault, np.arange(iFirstSys, iLastSys))
        # )[0]
        iSpeciesDBIndex = np.array([iFirstSys + Temp.index(x) for x in EndmembersName])
        dStoichSpeciesPhase = dStoichSpeciesPhase[iSpeciesDBIndex, :]
        ElementValues = np.matmul(endmemeber_amounts, dStoichSpeciesPhase)

    else:
        stripped_db_names = [name.strip() for name in var.cEndmemberNameCS]
        iSpeciesDBIndex = stripped_db_names.index(EndmembersName[0])
        dStoichSpeciesPhase = dStoichSpeciesPhase[iSpeciesDBIndex, :]
        ElementValues = endmemeber_amounts * dStoichSpeciesPhase

    # 2.3 Refine it for the system

    total_elements = np.sum(ElementValues)
    if np.isclose(total_elements, 0.0):
        return dict(zip(elements, np.zeros(len(elements), dtype=float), strict=False))
    ElementValues = ElementValues / total_elements

    if is_mass_based:
        ElementValues_mass = ElementValues * atomic_masses
        total_mass = np.sum(ElementValues_mass)
        if np.isclose(total_mass, 0.0):
            return dict(
                zip(elements, np.zeros(len(elements), dtype=float), strict=False)
            )
        ElementValues = ElementValues_mass / total_mass

    output = dict(zip(elements, ElementValues, strict=False))
    return output


def species_indices2elements(
    species_indices: np.ndarray,
    species_amounts: list[float],
    elements: list[str],
    unit_out: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Convert selected-system species amounts to element fractions."""
    if unit_out is None:
        unit_out = ["K", "atm", "moles"]

    if len(species_indices) == 0:
        return dict(zip(elements, np.zeros(len(elements), dtype=float), strict=False))

    iElementDBIndex = var.iElementDBIndex.copy()
    dStoichSpeciesPhase = var.dStoichSpeciesCS[:, iElementDBIndex]
    dStoichSpeciesPhase = dStoichSpeciesPhase[np.asarray(species_indices, dtype=int), :]
    element_values = np.matmul(
        np.asarray(species_amounts, dtype=float), dStoichSpeciesPhase
    )

    total_elements = np.sum(element_values)
    if np.isclose(total_elements, 0.0):
        return dict(zip(elements, np.zeros(len(elements), dtype=float), strict=False))
    element_values = element_values / total_elements

    is_mass_based = unit_out[2] in [
        "grams",
        "kilograms",
        "pounds",
        "g",
        "kg",
        "lbs",
        "mass fraction",
        "weight fraction",
        "wt%",
        "wt.%",
    ]
    if is_mass_based:
        atomic_masses = np.array(
            [var.cPeriodicTable[element.strip()][1] for element in elements]
        )
        element_values_mass = element_values * atomic_masses
        total_mass = np.sum(element_values_mass)
        if np.isclose(total_mass, 0.0):
            zeros = np.zeros(len(elements), dtype=float)
            return dict(zip(elements, zeros, strict=False))
        element_values = element_values_mass / total_mass

    return dict(zip(elements, element_values, strict=False))


def create_phase_from_sys(
    i_system: int,
    capture: Optional[FortranCaptureState] = None,
) -> PhaseResult:
    """
    Create a PhaseResult object from the Fortran system state.

    UPDATED: Creates a dictionary for xi.
    """
    if capture is None:
        capture = FortranCaptureState.from_globals()

    name = _public_phase_name(capture.phase_names[i_system])
    elements = list(capture.element_names)

    if i_system < capture.solution_count:
        # Solution PhaseResult
        phase_id = -(i_system + 1)
        assemblage_slot = capture.assemblage_index.get(phase_id)
        phase_metadata = _solution_phase_metadata(
            capture,
            system_phase_id=i_system + 1,
            slot_index=assemblage_slot,
        )
        if capture.phase_is_stable(phase_id):
            amount_n = capture.phase_amount_n(phase_id)
            amount_w = capture.phase_amount_w(phase_id)
            stability = 1.0
        else:
            amount_n = 0.0
            amount_w = 0.0
            stability = 0.0
        amount_n_basis = "solution_native"

        iFirstSys = capture.species_phase_boundaries[i_system]
        iLastSys = capture.species_phase_boundaries[i_system + 1]
        species_slice = slice(iFirstSys, iLastSys)

        idx_species = capture.system_to_database_species[species_slice]

        # --- MODIFIED ---
        # Create the xi dict by zipping names and values
        endmember_names = [capture.endmember_names[index] for index in idx_species]
        xi_values = list(capture.endmember_fractions_n[species_slice])
        wi_values = list(capture.endmember_fractions_w[species_slice])
        solution_type = str(var.cSolnPhaseTypeCS[int(var.iSys2DBSoln[i_system]) - 1])
        slot_elements_x: dict[str, float] | None = None
        if assemblage_slot is not None:
            slot_xi_values, slot_elements_x = _active_slot_solution_report_fractions(
                capture,
                system_phase_id=i_system + 1,
                slot_index=assemblage_slot,
                display_species_indices=idx_species,
                elements=elements,
            )
            if slot_xi_values is not None and len(slot_xi_values) == len(xi_values):
                xi_values = list(slot_xi_values)
                wi_values = _species_mass_fractions(idx_species, xi_values, elements)

        amount_n, amount_w, amount_n_basis = _solution_reported_phase_amounts(
            amount_n,
            amount_w,
            idx_species,
            xi_values,
            solution_type,
        )
        phase_atom_denominator = _phase_occupied_atom_count(
            iFirstSys,
            iLastSys,
            xi_values,
        )
        reported_partial_gibbs = (
            capture.partial_gibbs[species_slice] / phase_atom_denominator
        )
        reported_standard_gibbs = (
            capture.standard_gibbs_energies[species_slice] / phase_atom_denominator
        )
        reported_partial_enthalpy = (
            capture.partial_enthalpies[species_slice] / phase_atom_denominator
        )
        reported_partial_entropy = (
            capture.partial_entropies[species_slice] / phase_atom_denominator
        )
        reported_partial_heat_capacity = (
            capture.partial_heat_capacities[species_slice] / phase_atom_denominator
        )
        endmembers_x = dict(zip(endmember_names, xi_values, strict=False))
        endmembers_w = dict(zip(endmember_names, wi_values, strict=False))
        elements_x = (
            slot_elements_x
            if slot_elements_x is not None
            else species_indices2elements(idx_species, xi_values, elements)
        )
        elements_w = (
            _element_mass_fractions(elements_x)
            if slot_elements_x is not None
            else species_indices2elements(
                idx_species,
                xi_values,
                elements,
                unit_out=["K", "atm", "g"],
            )
        )
        if (
            solution_type.strip() not in {"SUBG", "SUBQ"}
            and amount_n_basis == "solution_native"
            and slot_elements_x is None
        ):
            phase_amount_from_mass = _phase_amount_from_mass_and_composition(
                amount_w,
                elements_x,
            )
            if phase_amount_from_mass is not None:
                amount_n = phase_amount_from_mass
        partial_gibbs = _phase_property_dict(
            endmember_names,
            reported_partial_gibbs,
        )
        standard_gibbs_energy = _phase_property_dict(
            endmember_names,
            reported_standard_gibbs,
        )
        activity = _phase_property_dict(
            endmember_names,
            capture.activities[species_slice],
        )
        partial_enthalpy = _phase_property_dict(
            endmember_names,
            reported_partial_enthalpy,
        )
        partial_entropy = _phase_property_dict(
            endmember_names,
            reported_partial_entropy,
        )
        partial_heat_capacity = _phase_property_dict(
            endmember_names,
            reported_partial_heat_capacity,
        )
        phase_G = _phase_weighted_property(
            xi_values,
            reported_partial_gibbs,
        )
        phase_H = _phase_weighted_property(
            xi_values,
            reported_partial_enthalpy,
        )
        phase_S = _phase_weighted_property(
            xi_values,
            reported_partial_entropy,
        )
        phase_Cp = _phase_weighted_property(
            xi_values,
            reported_partial_heat_capacity,
        )

    else:
        # Compound phase
        compound_species_index = int(
            capture.species_phase_boundaries[capture.solution_count]
            + (i_system - capture.solution_count)
        )
        stable_phase_id = _stable_phase_id_by_name(capture, name)
        if stable_phase_id is not None:
            phase_id = stable_phase_id
        else:
            # Compound assemblage IDs are 1-based system species indices.
            # Phase-row indices diverge after solution species are inserted.
            phase_id = compound_species_index + 1
        if capture.phase_is_stable(phase_id):
            amount_n = capture.phase_amount_n(phase_id)
            database_species_index = int(
                capture.system_to_database_species[compound_species_index],
            )
            amount_n, amount_n_basis = _compound_reported_phase_amount_n(
                amount_n,
                database_species_index,
            )
            amount_w = capture.phase_amount_w(phase_id)
            stability = 1.0
        else:
            amount_n = 0.0
            amount_w = 0.0
            stability = 0.0
            amount_n_basis = "formula_moles"

        # Stoichiometric compounds have fixed composition; do not expose
        # solution-style endmember/element fraction maps in result tables.
        endmembers_x = {}
        endmembers_w = {}
        elements_x = {}
        elements_w = {}
        phase_G = float(capture.partial_gibbs[compound_species_index])
        phase_H = float(capture.partial_enthalpies[compound_species_index])
        phase_S = float(capture.partial_entropies[compound_species_index])
        phase_Cp = float(capture.partial_heat_capacities[compound_species_index])
        partial_gibbs = {}
        standard_gibbs_energy = {}
        activity = {}
        partial_enthalpy = {}
        partial_entropy = {}
        partial_heat_capacity = {}
        phase_metadata = {
            "parent_model_id": None,
            "parent_model_name": None,
            "composition_set_id": 0,
            "display_label": name.strip(),
            "ordering_degree": None,
        }

    return PhaseResult(
        id=phase_id,
        name=name.strip(),
        amount_n=amount_n,
        amount_w=amount_w,
        amount_n_basis=amount_n_basis,
        stability=stability,
        parent_model_id=phase_metadata["parent_model_id"],
        parent_model_name=phase_metadata["parent_model_name"],
        composition_set_id=phase_metadata["composition_set_id"],
        display_label=phase_metadata["display_label"],
        ordering_degree=phase_metadata["ordering_degree"],
        endmembers_x=endmembers_x,
        endmembers_w=endmembers_w,
        elements_x=elements_x,
        elements_w=elements_w,
        G=phase_G,
        H=phase_H,
        S=phase_S,
        Cp=phase_Cp,
        partial_gibbs=partial_gibbs,
        standard_gibbs_energy=standard_gibbs_energy,
        activity=activity,
        partial_enthalpy=partial_enthalpy,
        partial_entropy=partial_entropy,
        partial_heat_capacity=partial_heat_capacity,
    )


class StablePhaseSummaryDict(TypedDict, total=False):
    """Typed dictionary describing stable phase arrays."""

    id: np.ndarray
    name: np.ndarray
    amount_n: np.ndarray
    amount_w: np.ndarray
    amount_n_basis: np.ndarray
    parent_model_id: np.ndarray
    parent_model_name: np.ndarray
    composition_set_id: np.ndarray
    display_label: np.ndarray
    ordering_degree: np.ndarray


# Type alias for the dictionary of phase results.
# Keys are dynamic phase names (e.g. 'LIQUID', 'FCC_A1').
PhaseDict = Dict[str, PhaseResult]


class EquilibPoint(BaseModel):
    """
    Pydantic model for a single equilibrium calculation point.

    This represents one "row" of data.
    """

    T: float
    P: float
    G: Optional[float] = None
    H: Optional[float] = None
    S: Optional[float] = None
    Cp: Optional[float] = None
    n_i: Dict[str, float] = Field(default_factory=dict)
    w_i: Dict[str, float] = Field(default_factory=dict)
    stable_phase_summary: StablePhaseSummaryDict = Field(default_factory=dict)
    phase_map: Dict[str, PhaseResult] = Field(default_factory=dict)

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        extra="forbid",
        populate_by_name=True,
    )

    @model_validator(mode="after")
    def _normalize_amount_maps(self) -> "EquilibPoint":
        """Normalize amount mappings to clean component names."""
        self.n_i = clean_quantity_mapping(self.n_i)
        self.w_i = clean_quantity_mapping(self.w_i)
        return self

    @property
    def phases(self) -> PhaseCollection:
        """Return all phases for this equilibrium point."""
        return PhaseCollection(self.phase_map)

    @property
    def stable_phases(self) -> PhaseCollection:
        """Return stable phases for this equilibrium point."""
        return self.phases.stable()

    def phase(self, name: str) -> PhaseResult:
        """Return one phase by name."""
        return self.phases[name]

    @classmethod
    def from_fortran(cls) -> "EquilibPoint":
        """Create an instance from the Fortran and variable modules."""
        try:
            capture = FortranCaptureState.from_globals()
            G = float(fort.modulethermoio.dgibbsenergysys)
            H = float(fort.modulethermoio.denthalpysys)
            S = float(fort.modulethermoio.dentropysys)
            Cp = float(fort.modulethermoio.dheatcapacitysys)

            var.iPhaseSys = []
            all_phases = {}
            for i in range(capture.phase_count):
                new_phase = create_phase_from_sys(i, capture)
                if _prefer_phase_result(new_phase, all_phases.get(new_phase.name)):
                    all_phases[new_phase.name] = new_phase
                var.iPhaseSys.append(new_phase.id)

            stable_names = np.array(list(get_assemblage_name(capture.assemblage_ids)))
            stable_phase_summary = StablePhaseSummaryDict(
                {
                    "id": capture.assemblage_ids,
                    "name": stable_names,
                    "amount_n": _stable_phase_amounts_n(
                        capture,
                        stable_names,
                        all_phases,
                    ),
                    "amount_w": capture.phase_amounts_w,
                    "amount_n_basis": _stable_phase_amount_basis(
                        stable_names,
                        all_phases,
                    ),
                    **_stable_phase_metadata(capture),
                }
            )

            return cls(
                T=float(fort.modulethermoio.dtemperature),
                P=float(fort.modulethermoio.dpressure),
                n_i=dict(
                    zip(
                        var.cComponentNameSys,
                        fort.modulethermo.dmoleselement.copy()[var.iElementSysIndex],
                        strict=False,
                    )
                ),
                w_i=dict(
                    zip(
                        var.cComponentNameSys,
                        fort.modulethermoio.dgramelement.copy()[var.iElementSysIndex],
                        strict=False,
                    )
                ),
                G=G,
                H=H,
                S=S,
                Cp=Cp,
                stable_phase_summary=stable_phase_summary,
                phase_map=all_phases,
            )
        except Exception as e:
            print(f"Error reading from Fortran, creating error row: {e}")
            return cls.for_error()

    @classmethod
    def for_error(cls) -> "EquilibPoint":
        """
        Create an instance representing a failed calculation.

        The failed calculation row is populated with NaNs.
        """
        stable_phases_dict = {
            "id": np.array([np.nan]),
            "name": np.array(["nan"]),
            "amount_n": np.array([np.nan]),
            "amount_w": np.array([np.nan]),
            "amount_n_basis": np.array(["unknown"]),
            "parent_model_id": np.array([np.nan]),
            "parent_model_name": np.array([""]),
            "composition_set_id": np.array([0]),
            "display_label": np.array(["nan"]),
            "ordering_degree": np.array([np.nan]),
        }
        component_names = list(getattr(var, "cComponentNameSys", []))
        element_index = getattr(var, "iElementSysIndex", [])
        component_count = len(component_names)
        mole_values = safe_error_component_values(
            getattr(fort.modulethermo, "dmoleselement", None),
            element_index,
            component_count,
        )
        mass_values = safe_error_component_values(
            getattr(fort.modulethermoio, "dgramelement", None),
            element_index,
            component_count,
        )

        return cls(
            T=float(fort.modulethermoio.dtemperature),
            P=float(fort.modulethermoio.dpressure),
            n_i=dict(
                zip(
                    component_names,
                    mole_values,
                    strict=False,
                )
            ),
            w_i=dict(
                zip(
                    component_names,
                    mass_values,
                    strict=False,
                )
            ),
            G=np.nan,
            H=np.nan,
            S=np.nan,
            Cp=np.nan,
            stable_phase_summary=stable_phases_dict,
            phase_map={},
        )


class EquilibResult:
    """
    Container for one or more equilibrium calculation points.

    The public API normalizes single and batch results through ``point`` and
    ``points`` while preserving scalar convenience properties for one-point
    results.
    """

    def __init__(self, context: Optional[ResultContext] = None):
        self.context = context
        self._data: Union[None, EquilibPoint, List[EquilibPoint]] = None
        self._phases_cache: PhaseCollection | None = None
        self._stable_phases_cache: PhaseCollection | None = None
        self._phase_species_property_cache: dict[str, PhaseSpeciesProperty] = {}

    @property
    def data(self) -> Union[None, EquilibPoint, List[EquilibPoint]]:
        """Return raw point storage for single or batch equilibrium results."""
        return self._data

    @data.setter
    def data(self, value: Union[None, EquilibPoint, List[EquilibPoint]]) -> None:
        self._data = list(value) if isinstance(value, list) else value
        self._clear_view_cache()

    def _build_phase_collection(self) -> PhaseCollection:
        """Build the materialized phase view from current raw data."""
        if self._data is None:
            return PhaseCollection({})
        if isinstance(self._data, EquilibPoint):
            return self._data.phases
        return PhaseCollection(_batch_phase_result_map(list(self._data)))

    @staticmethod
    def _phase_species_property_from_phases(
        phases: PhaseCollection,
        attribute_name: str,
    ) -> PhaseSpeciesProperty:
        """Build a materialized phase-first view of one species property."""
        return PhaseSpeciesProperty(
            {
                phase_name: getattr(phase, attribute_name, {})
                for phase_name, phase in phases.items()
            }
        )

    def _clear_view_cache(self) -> None:
        """Clear lazily materialized batch phase/property views."""
        self._phases_cache = None
        self._stable_phases_cache = None
        self._phase_species_property_cache = {}

    def _refresh_views(self) -> None:
        """Clear materialized views; public properties rebuild them lazily."""
        self._clear_view_cache()

    @property
    def phases(self) -> PhaseCollection:
        """Return phase results, materialized lazily for batch results."""
        if self._phases_cache is None:
            self._phases_cache = self._build_phase_collection()
        return self._phases_cache

    @property
    def stable_phases(self) -> PhaseCollection:
        """Return stable phase results, materialized lazily for batch results."""
        if self._stable_phases_cache is None:
            self._stable_phases_cache = self.phases.stable()
        return self._stable_phases_cache

    def _phase_species_property(self, attribute_name: str) -> PhaseSpeciesProperty:
        """Return one lazy phase/species thermodynamic-property view."""
        if attribute_name not in self._phase_species_property_cache:
            self._phase_species_property_cache[attribute_name] = (
                self._phase_species_property_from_phases(
                    self.phases,
                    attribute_name,
                )
            )
        return self._phase_species_property_cache[attribute_name]

    @property
    def partial_gibbs(self) -> PhaseSpeciesProperty:
        """Return partial molar Gibbs energies keyed by phase and species."""
        return self._phase_species_property(
            "partial_gibbs",
        )

    @property
    def standard_gibbs_energy(self) -> PhaseSpeciesProperty:
        """Return standard Gibbs energies keyed by phase and species."""
        return self._phase_species_property(
            "standard_gibbs_energy",
        )

    @property
    def activity(self) -> PhaseSpeciesProperty:
        """Return activities keyed by phase and species."""
        return self._phase_species_property(
            "activity",
        )

    @property
    def partial_enthalpy(self) -> PhaseSpeciesProperty:
        """Return partial enthalpies keyed by phase and species."""
        return self._phase_species_property(
            "partial_enthalpy",
        )

    @property
    def partial_entropy(self) -> PhaseSpeciesProperty:
        """Return partial entropies keyed by phase and species."""
        return self._phase_species_property(
            "partial_entropy",
        )

    @property
    def partial_heat_capacity(self) -> PhaseSpeciesProperty:
        """Return partial heat capacities keyed by phase and species."""
        return self._phase_species_property(
            "partial_heat_capacity",
        )

    @property
    def points(self) -> List[EquilibPoint]:
        """Return equilibrium points as a list regardless of internal state."""
        if self.data is None:
            return []
        if isinstance(self.data, EquilibPoint):
            return [self.data]
        return list(self.data)

    @property
    def point(self) -> EquilibPoint:
        """Return the only equilibrium point, or raise for empty/batch results."""
        points = self.points
        if len(points) != 1:
            raise PostProcessError(
                "Expected exactly one equilibrium point; "
                f"found {len(points)}."
            )
        return points[0]

    def phase(self, name: str) -> PhaseResult | BatchPhaseResult:
        """Return one phase by name for single or batch equilibrium results."""
        return self.phases[name]

    def append(self, other: "EquilibResult"):
        """Append another result object into this result."""
        if self.context is None:
            self.context = other.context
        to_append = other.points
        if not to_append:
            return
        if self.data is None:
            self.data = to_append[0] if len(to_append) == 1 else list(to_append)
        elif isinstance(self.data, EquilibPoint):
            self.data = [self.data] + to_append
        elif isinstance(self.data, list):
            self.data = [*self.data, *to_append]

    def append_output(self):
        """Append the current Fortran output state."""
        new_point = EquilibPoint.from_fortran()
        if self.data is None:
            self.data = new_point
        elif isinstance(self.data, EquilibPoint):
            self.data = [self.data, new_point]
        elif isinstance(self.data, list):
            self.data = [*self.data, new_point]

    def append_error(self):
        """Append an error result row."""
        error_point = EquilibPoint.for_error()
        if self.data is None:
            self.data = error_point
        elif isinstance(self.data, EquilibPoint):
            self.data = [self.data, error_point]
        elif isinstance(self.data, list):
            self.data = [*self.data, error_point]

    @property
    def n_i(self) -> Union[Dict[str, float], List[Dict[str, float]]]:
        """Return system component amounts on a mole basis."""
        points = self.points
        if not points:
            return {}
        if len(points) == 1 and isinstance(self.data, EquilibPoint):
            return points[0].n_i
        return [point.n_i for point in points]

    @property
    def w_i(self) -> Union[Dict[str, float], List[Dict[str, float]]]:
        """Return system component amounts on a mass basis."""
        points = self.points
        if not points:
            return {}
        if len(points) == 1 and isinstance(self.data, EquilibPoint):
            return points[0].w_i
        return [point.w_i for point in points]

    @property
    def T(self) -> Union[float, List[float]]:
        """Return temperature values."""
        if isinstance(self.data, EquilibPoint):
            return self.data.T
        if isinstance(self.data, list):
            return [iter.T for iter in self.data]
        return None

    @property
    def P(self) -> Union[float, List[float]]:
        """Return pressure values."""
        if isinstance(self.data, EquilibPoint):
            return self.data.P
        if isinstance(self.data, list):
            return [iter.P for iter in self.data]
        return None

    @property
    def G(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return Gibbs energy values."""
        if isinstance(self.data, EquilibPoint):
            return self.data.G
        if isinstance(self.data, list):
            return [iter.G for iter in self.data]
        return None

    @property
    def H(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return enthalpy values."""
        if isinstance(self.data, EquilibPoint):
            return self.data.H
        if isinstance(self.data, list):
            return [iter.H for iter in self.data]
        return None

    @property
    def S(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return entropy values."""
        if isinstance(self.data, EquilibPoint):
            return self.data.S
        if isinstance(self.data, list):
            return [iter.S for iter in self.data]
        return None

    @property
    def Cp(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return heat capacity values."""
        if isinstance(self.data, EquilibPoint):
            return self.data.Cp
        if isinstance(self.data, list):
            return [iter.Cp for iter in self.data]
        return None

    def _discover_phase_sub_keys(
        self, iterations: List[EquilibPoint], attribute_name: str
    ) -> Dict[str, set]:
        """
        Discover all unique sub-keys for a phase attribute dictionary.

        e.g., for 'endmembers_x', returns:
        {'LIQUID': {'AL', 'FE'}, 'FCC_A1': {'AL', 'VA'}}
        """
        all_sub_keys = defaultdict(set)
        for iter_data in iterations:
            for phase_name, phase_obj in iter_data.phase_map.items():
                attribute_dict = getattr(phase_obj, attribute_name, {})
                if attribute_dict:
                    all_sub_keys[phase_name].update(attribute_dict.keys())
        return all_sub_keys

    def _populate_flattened_phase_attribute(
        self,
        dict_of_lists: defaultdict,
        phase_name: str,
        phase_obj: Optional[PhaseResult],
        sub_keys: set,
        attribute_name: str,
        column_suffix: str,
    ):
        """
        Populate dict_of_lists for one phase's flattened attribute.

        e.g., for phase 'LIQUID', attribute 'xi', suffix '_xi'.
        """
        for component in sub_keys:
            output_key = f"{phase_name}{column_suffix}_{component}"

            if phase_obj:
                attribute_dict = getattr(phase_obj, attribute_name, {})
                value = attribute_dict.get(component, 0.0)
            else:
                value = 0.0

            dict_of_lists[output_key].append(value)

    def _populate_flattened_phase_property(
        self,
        dict_of_lists: defaultdict,
        phase_name: str,
        phase_obj: Optional[PhaseResult],
        sub_keys: set,
        attribute_name: str,
        column_suffix: str,
        unit: str,
    ):
        """Populate a flattened phase species property with a unit suffix."""
        for component in sub_keys:
            output_key = f"{phase_name}{column_suffix}_{component} [{unit}]"
            if phase_obj:
                attribute_dict = getattr(phase_obj, attribute_name, {})
                value = attribute_dict.get(component, 0.0)
            else:
                value = 0.0
            dict_of_lists[output_key].append(value)

    def to_dict(self) -> Dict[str, Any]:
        """
        Export all results to a flattened dictionary.

        Uses reusable helpers to flatten phase attributes.
        """
        points = self.points
        if not points:
            return {}

        iterations = points
        is_single = len(points) == 1 and isinstance(self.data, EquilibPoint)

        dict_of_lists = defaultdict(list)

        all_phase_keys = []
        endmembers_x = self._discover_phase_sub_keys(iterations, "endmembers_x")
        endmembers_w = self._discover_phase_sub_keys(iterations, "endmembers_w")
        elements_x = self._discover_phase_sub_keys(iterations, "elements_x")
        elements_w = self._discover_phase_sub_keys(iterations, "elements_w")
        partial_gibbs = self._discover_phase_sub_keys(
            iterations,
            "partial_gibbs",
        )
        standard_gibbs_energy = self._discover_phase_sub_keys(
            iterations,
            "standard_gibbs_energy",
        )
        activity = self._discover_phase_sub_keys(
            iterations,
            "activity",
        )
        partial_enthalpy = self._discover_phase_sub_keys(
            iterations,
            "partial_enthalpy",
        )
        partial_entropy = self._discover_phase_sub_keys(
            iterations,
            "partial_entropy",
        )
        partial_heat_capacity = self._discover_phase_sub_keys(
            iterations,
            "partial_heat_capacity",
        )

        amount_maps = []
        for iter_data in iterations:
            amount_maps.extend([iter_data.n_i, iter_data.w_i])
            all_phase_keys.extend(iter_data.phase_map.keys())
        all_comp_keys = component_names_from_amount_maps(self.context, *amount_maps)
        all_phase_keys = unique_preserve_order(all_phase_keys)

        for iter_data in iterations:
            dict_of_lists["T [K]"].append(iter_data.T)
            dict_of_lists["P [atm]"].append(iter_data.P)
            for key in all_comp_keys:
                dict_of_lists[f"n_{key} [sp-mol]"].append(
                    iter_data.n_i.get(key, 0.0)
                )
            for key in all_comp_keys:
                dict_of_lists[f"w_{key} [g]"].append(
                    iter_data.w_i.get(key, 0.0)
                )
            dict_of_lists["G [J]"].append(iter_data.G)
            dict_of_lists["H [J]"].append(iter_data.H)
            dict_of_lists["S [J/K]"].append(iter_data.S)
            dict_of_lists["Cp [J/K]"].append(iter_data.Cp)

            stable_summary = iter_data.stable_phase_summary
            dict_of_lists["stable_phase_names"].append(
                str(stable_summary.get("name", np.array([])))
            )
            dict_of_lists["stable_phase_ids"].append(
                str(stable_summary.get("id", np.array([])))
            )
            dict_of_lists["stable_phase_amount_n [sp-mol]"].append(
                str(stable_summary.get("amount_n", np.array([])))
            )
            dict_of_lists["stable_phase_amount_w [g]"].append(
                str(stable_summary.get("amount_w", np.array([])))
            )
            dict_of_lists["stable_phase_amount_n_basis"].append(
                str(stable_summary.get("amount_n_basis", np.array([])))
            )

            for key in all_phase_keys:
                phase = iter_data.phase_map.get(key)

                if phase:
                    dict_of_lists[f"{key}_amount_n [sp-mol]"].append(
                        phase.amount_n
                    )
                    dict_of_lists[f"{key}_amount_w [g]"].append(phase.amount_w)
                    dict_of_lists[f"{key}_amount_n_basis"].append(
                        phase.amount_n_basis
                    )
                    dict_of_lists[f"{key}_stability"].append(phase.stability)
                else:
                    dict_of_lists[f"{key}_amount_n [sp-mol]"].append(0.0)
                    dict_of_lists[f"{key}_amount_w [g]"].append(0.0)
                    dict_of_lists[f"{key}_amount_n_basis"].append("")
                    dict_of_lists[f"{key}_stability"].append(0.0)

                self._populate_flattened_phase_attribute(
                    dict_of_lists,
                    key,
                    phase,
                    endmembers_x[key],
                    "endmembers_x",
                    "_endmembers_x",
                )
                self._populate_flattened_phase_attribute(
                    dict_of_lists,
                    key,
                    phase,
                    endmembers_w[key],
                    "endmembers_w",
                    "_endmembers_w",
                )
                self._populate_flattened_phase_attribute(
                    dict_of_lists,
                    key,
                    phase,
                    elements_x[key],
                    "elements_x",
                    "_elements_x",
                )
                self._populate_flattened_phase_attribute(
                    dict_of_lists,
                    key,
                    phase,
                    elements_w[key],
                    "elements_w",
                    "_elements_w",
                )
                self._populate_flattened_phase_property(
                    dict_of_lists,
                    key,
                    phase,
                    partial_gibbs[key],
                    "partial_gibbs",
                    "_partial_gibbs",
                    "J",
                )
                self._populate_flattened_phase_property(
                    dict_of_lists,
                    key,
                    phase,
                    standard_gibbs_energy[key],
                    "standard_gibbs_energy",
                    "_standard_gibbs_energy",
                    "J",
                )
                self._populate_flattened_phase_property(
                    dict_of_lists,
                    key,
                    phase,
                    activity[key],
                    "activity",
                    "_activity",
                    "-",
                )
                self._populate_flattened_phase_property(
                    dict_of_lists,
                    key,
                    phase,
                    partial_enthalpy[key],
                    "partial_enthalpy",
                    "_partial_enthalpy",
                    "J",
                )
                self._populate_flattened_phase_property(
                    dict_of_lists,
                    key,
                    phase,
                    partial_entropy[key],
                    "partial_entropy",
                    "_partial_entropy",
                    "J/K",
                )
                self._populate_flattened_phase_property(
                    dict_of_lists,
                    key,
                    phase,
                    partial_heat_capacity[key],
                    "partial_heat_capacity",
                    "_partial_heat_capacity",
                    "J/K",
                )

        final_dict = dict(dict_of_lists)

        if is_single:
            for key, value_list in final_dict.items():
                final_dict[key] = value_list[0] if value_list else None

        return final_dict

    def to_table(self) -> ResultTable:
        """Return a selectable/exportable table view of this result."""
        return ResultTable.from_dict(self.to_dict())

    def available_columns(self) -> List[str]:
        """Return flattened result columns available for export or display."""
        return self.to_table().available_columns()

    def to_bundle(self) -> Dict[str, Any]:
        """Return a JSON-safe bundle that can reconstruct this result."""
        points = self.points
        if not points:
            data_state = "none"
        elif isinstance(self.data, EquilibPoint):
            data_state = "single"
        else:
            data_state = "list"

        return {
            "format": "equilipy.result",
            "version": 1,
            "kind": "equilib",
            "context": context_to_payload(self.context),
            "data_state": data_state,
            "data": [single_point_to_payload(point) for point in points],
        }

    @classmethod
    def from_bundle(cls, bundle: Dict[str, Any]) -> "EquilibResult":
        """Reconstruct a result from a JSON-safe bundle."""
        validate_result_bundle(bundle, EQUILIB_BUNDLE_KINDS)
        result = cls(context=context_from_payload(bundle.get("context")))
        points = [
            single_point_from_payload(point)
            for point in bundle.get("data", [])
            if isinstance(point, dict)
        ]
        data_state = bundle.get("data_state")
        if data_state == "single":
            result.data = points[0] if points else None
        elif data_state == "list":
            result.data = points
        elif data_state in {"none", None}:
            result.data = None
        else:
            raise PostProcessError(f"Unsupported result data state: {data_state!r}.")
        return result


__all__ = [
    "BatchPhaseResult",
    "PhaseDict",
    "EquilibPoint",
    "EquilibResult",
    "PhaseResult",
    "StablePhaseSummaryDict",
    "count_unstable_compounds",
    "create_phase_from_sys",
    "endmembers2elements",
    "get_assemblage_name",
    "safe_error_component_values",
]
