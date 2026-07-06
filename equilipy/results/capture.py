"""Snapshot helpers for reading Fortran-backed result state."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var
from equilipy.results.common import normalize_phase_amounts


@dataclass(frozen=True)
class FortranCaptureState:
    """Frozen snapshot of the Fortran and Python global result state."""

    assemblage_ids: np.ndarray
    assemblage_index: dict[int, int]
    phase_names: tuple[str, ...]
    solution_count: int
    compound_count: int
    species_phase_boundaries: np.ndarray
    system_to_database_species: np.ndarray
    endmember_names: tuple[str, ...]
    element_names: tuple[str, ...]
    phase_amounts_n: np.ndarray
    phase_amounts_w: np.ndarray
    endmember_fractions_n: np.ndarray
    endmember_fractions_w: np.ndarray
    partial_gibbs: np.ndarray
    standard_gibbs_energies: np.ndarray
    activities: np.ndarray
    partial_enthalpies: np.ndarray
    partial_entropies: np.ndarray
    partial_heat_capacities: np.ndarray
    active_slot_thermo_phase: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=int),
    )
    active_slot_display_phase: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=int),
    )
    active_slot_identity_ordinal: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=int),
    )
    active_slot_site_fraction: np.ndarray = field(
        default_factory=lambda: np.zeros((0, 0, 0), dtype=float),
    )

    def __post_init__(self) -> None:
        """Normalize tiny Fortran phase amount noise on capture."""
        object.__setattr__(
            self,
            "phase_amounts_n",
            normalize_phase_amounts(self.phase_amounts_n),
        )
        object.__setattr__(
            self,
            "phase_amounts_w",
            normalize_phase_amounts(self.phase_amounts_w),
        )

    @classmethod
    def from_globals(cls) -> "FortranCaptureState":
        """Capture the current global result state after a solver call."""
        assemblage_ids = np.array(fort.modulethermo.iassemblage).copy()
        assemblage_index = {
            int(phase_id): index
            for index, phase_id in enumerate(assemblage_ids)
        }
        element_indices = np.array(var.iElementDBIndex)
        element_names = tuple(
            str(name).strip()
            for name in np.array(var.cElementNameCS)[element_indices]
        )
        ideal_constant = float(fort.modulethermo.didealconstant)
        temperature = float(fort.modulethermoio.dtemperature)
        energy_factor = ideal_constant * temperature
        partial_gibbs = np.array(fort.modulethermo.dchemicalpotential).copy()
        standard_gibbs_energies = np.array(fort.modulethermo.dstdgibbsenergy).copy()
        active_slot_count = len(assemblage_ids)
        try:
            active_slot_thermo_phase = np.array(
                fort.modulegemsolver.iactiveslotthermophase,
                dtype=int,
            ).copy()
        except (AttributeError, ValueError):
            active_slot_thermo_phase = np.zeros(active_slot_count, dtype=int)
        try:
            active_slot_display_phase = np.array(
                fort.modulegemsolver.iactiveslotdisplayphase,
                dtype=int,
            ).copy()
        except (AttributeError, ValueError):
            active_slot_display_phase = np.zeros(active_slot_count, dtype=int)
        try:
            active_slot_identity_ordinal = np.array(
                fort.modulegemsolver.iactiveslotidentityordinal,
                dtype=int,
            ).copy()
        except (AttributeError, ValueError):
            active_slot_identity_ordinal = np.ones(active_slot_count, dtype=int)
        try:
            active_slot_site_fraction = np.array(
                fort.modulegemsolver.dactiveslotsitefraction,
                dtype=float,
            ).copy()
        except (AttributeError, ValueError):
            active_slot_site_fraction = np.zeros((active_slot_count, 0, 0), dtype=float)
        try:
            activities = np.array(fort.modulethermo.dactivity).copy()
        except (AttributeError, ValueError):
            # Compatibility for an extension built before dActivity existed.
            activities = np.exp(
                np.clip(
                    partial_gibbs[: len(standard_gibbs_energies)]
                    - standard_gibbs_energies,
                    -745.0,
                    709.0,
                )
            )

        return cls(
            assemblage_ids=assemblage_ids,
            assemblage_index=assemblage_index,
            phase_names=tuple(str(name).strip() for name in var.cPhaseNameSys),
            solution_count=len(var.iSys2DBSoln),
            compound_count=len(var.iSys2DBComp),
            species_phase_boundaries=np.array(
                fort.modulethermo.nspeciesphase
            ).copy(),
            system_to_database_species=np.array(var.iSys2DBSpecies).copy(),
            endmember_names=tuple(
                str(name).strip() for name in np.array(var.cEndmemberNameCS)
            ),
            element_names=element_names,
            phase_amounts_n=np.array(fort.modulethermo.dmolesphase).copy(),
            phase_amounts_w=np.array(fort.modulethermoio.dgramphase).copy(),
            endmember_fractions_n=np.array(fort.modulethermo.dmolfraction).copy(),
            endmember_fractions_w=np.array(
                fort.modulethermoio.dgramfraction
            ).copy(),
            partial_gibbs=partial_gibbs * energy_factor,
            standard_gibbs_energies=standard_gibbs_energies * energy_factor,
            activities=activities,
            partial_enthalpies=np.array(fort.modulethermo.dpartialenthalpy).copy()
            * energy_factor,
            partial_entropies=np.array(fort.modulethermo.dpartialentropy).copy()
            * ideal_constant,
            partial_heat_capacities=np.array(
                fort.modulethermo.dpartialheatcapacity
            ).copy()
            * ideal_constant,
            active_slot_thermo_phase=active_slot_thermo_phase,
            active_slot_display_phase=active_slot_display_phase,
            active_slot_identity_ordinal=active_slot_identity_ordinal,
            active_slot_site_fraction=active_slot_site_fraction,
        )

    @property
    def phase_count(self) -> int:
        """Return the number of phases in the selected system."""
        return self.solution_count + self.compound_count

    def phase_is_stable(self, phase_id: int) -> bool:
        """Return whether a phase ID is in the stable assemblage."""
        return int(phase_id) in self.assemblage_index

    def phase_amount_n(self, phase_id: int) -> float:
        """Return stable phase amount on a mole basis."""
        index = self.assemblage_index.get(int(phase_id))
        if index is None:
            return 0.0
        return float(self.phase_amounts_n[index])

    def phase_amount_w(self, phase_id: int) -> float:
        """Return stable phase amount on a mass basis."""
        index = self.assemblage_index.get(int(phase_id))
        if index is None:
            return 0.0
        return float(self.phase_amounts_w[index])
