"""Phase collection helpers for result objects."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class CompositionFractions:
    """Mole- and mass-basis composition fractions for a phase surface."""

    x_i: dict[str, float]
    w_i: dict[str, float]


class PhaseCollection(dict):
    """Ordered dictionary of phase names to phase-like result objects."""

    def __init__(self, phases: Mapping[str, Any] | None = None) -> None:
        super().__init__(dict(phases or {}))

    @property
    def names(self) -> list[str]:
        """Return phase names in deterministic insertion order."""
        return list(self.keys())

    @property
    def stable_names(self) -> list[str]:
        """Return names for phases whose stability flag is positive."""
        return [
            name
            for name, phase in self.items()
            if _has_positive_stability(getattr(phase, "stability", 0.0))
        ]

    def stable(self) -> "PhaseCollection":
        """Return a collection containing only stable phases."""
        return PhaseCollection(
            {
                name: self[name]
                for name in self.stable_names
            }
        )

    def __repr__(self) -> str:
        """Return a readable phase-name summary for interactive sessions."""
        return f"{type(self).__name__}(names={self.names!r})"


class PhaseSpeciesProperty(dict):
    """Phase-first view of one species-level property."""

    def __init__(self, phases: Mapping[str, Mapping[str, Any]] | None = None) -> None:
        super().__init__(
            {
                phase_name: dict(values)
                for phase_name, values in dict(phases or {}).items()
            }
        )

    @property
    def names(self) -> list[str]:
        """Return phase names in deterministic insertion order."""
        return list(self.keys())

    @property
    def species_names(self) -> list[str]:
        """Return species names in deterministic insertion order across phases."""
        names: list[str] = []
        seen: set[str] = set()
        for values in self.values():
            for species_name in values:
                if species_name not in seen:
                    seen.add(species_name)
                    names.append(species_name)
        return names

    def to_numpy(self) -> np.ndarray:
        """
        Return a dense numerical view with missing entries as NaN.

        Shape is ``(n_phases, n_species)`` for scalar values and
        ``(n_points, n_phases, n_species)`` when stored values are vectors.
        """
        phase_names = self.names
        species_names = self.species_names
        if not phase_names or not species_names:
            return np.empty((0, 0), dtype=float)

        sample_shape: tuple[int, ...] = ()
        for values in self.values():
            for value in values.values():
                array = np.asarray(value, dtype=float)
                sample_shape = array.shape
                break
            if sample_shape:
                break

        if sample_shape:
            output = np.full(
                (*sample_shape, len(phase_names), len(species_names)),
                np.nan,
                dtype=float,
            )
        else:
            output = np.full(
                (len(phase_names), len(species_names)),
                np.nan,
                dtype=float,
            )

        for phase_index, phase_name in enumerate(phase_names):
            values = self[phase_name]
            for species_index, species_name in enumerate(species_names):
                if species_name not in values:
                    continue
                value = np.asarray(values[species_name], dtype=float)
                if sample_shape:
                    output[..., phase_index, species_index] = value
                else:
                    output[phase_index, species_index] = float(value)
        return output

    def __getitem__(self, key: str) -> Any:
        """Return a phase property map or a flat ``species@phase`` value."""
        if "@" in key:
            species_name, phase_name = key.split("@", 1)
            return self[phase_name][species_name]
        return super().__getitem__(key)

    def __repr__(self) -> str:
        """Return a readable phase-name summary for interactive sessions."""
        return f"{type(self).__name__}(names={self.names!r})"


def _has_positive_stability(value: Any) -> bool:
    if isinstance(value, (list, tuple)):
        return any(float(item or 0.0) > 0.0 for item in value)
    return float(value or 0.0) > 0.0
