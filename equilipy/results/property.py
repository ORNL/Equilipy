"""Property-result models for fixed or minimized single phases."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from .context import ResultContext
from .phase import PhaseCollection


class SpeciesProperty(Mapping[str, Any]):
    """Species-keyed vector or matrix with a stable numerical view."""

    def __init__(self, names: list[str] | tuple[str, ...], values: Any) -> None:
        self._names = tuple(str(name).strip() for name in names)
        self._values = np.asarray(values, dtype=float).copy()
        if self._values.ndim == 0:
            self._values = self._values.reshape(1)
        if self._values.shape[-1] != len(self._names):
            raise ValueError(
                "SpeciesProperty values must have the species dimension last."
            )
        self._index = {name: index for index, name in enumerate(self._names)}

    @property
    def names(self) -> list[str]:
        """Return species names in column order."""
        return list(self._names)

    def to_numpy(self) -> np.ndarray:
        """Return values as a NumPy vector or matrix."""
        return self._values.copy()

    def to_dict(self) -> dict[str, Any]:
        """Return a regular dictionary preserving scalar/vector shape."""
        return {name: self[name] for name in self._names}

    def __getitem__(self, key: str) -> Any:
        """Return one species column as a scalar or vector."""
        values = self._values[..., self._index[str(key).strip()]]
        if np.ndim(values) == 0:
            return float(values)
        return np.asarray(values, dtype=float).copy()

    def __iter__(self) -> Iterator[str]:
        """Iterate over species names."""
        return iter(self._names)

    def __len__(self) -> int:
        """Return the number of species columns."""
        return len(self._names)

    def __repr__(self) -> str:
        """Return a compact interactive representation."""
        return (
            f"{type(self).__name__}(names={self.names!r}, "
            f"shape={self._values.shape!r})"
        )


@dataclass(frozen=True)
class SlnPropertyPhase:
    """Phase-like view for one solution-property result."""

    name: str
    stability: float
    endmembers_x: dict[str, float]
    endmembers_w: dict[str, float]
    elements_x: dict[str, float]
    elements_w: dict[str, float]
    G: float | None
    H: float | None
    S: float | None
    Cp: float | None
    partial_gibbs: SpeciesProperty
    partial_enthalpy: SpeciesProperty
    partial_entropy: SpeciesProperty
    activity: SpeciesProperty
    partial_heat_capacity: SpeciesProperty


@dataclass(frozen=True)
class SlnPropertyPoint:
    """One solution-property calculation point."""

    T: float
    P: float
    phase_name: str
    G: float | None
    H: float | None
    S: float | None
    Cp: float | None
    n_i: dict[str, float]
    w_i: dict[str, float]
    endmembers_x: dict[str, float]
    endmembers_w: dict[str, float]
    elements_x: dict[str, float]
    elements_w: dict[str, float]
    partial_gibbs: SpeciesProperty
    partial_enthalpy: SpeciesProperty
    partial_entropy: SpeciesProperty
    activity: SpeciesProperty
    partial_heat_capacity: SpeciesProperty
    dg_nn: np.ndarray | None = None

    def phase(self, name: str) -> SlnPropertyPhase:
        """Return a lightweight phase view for compatibility."""
        if str(name).strip() != self.phase_name:
            raise KeyError(name)
        return SlnPropertyPhase(
            name=self.phase_name,
            stability=1.0,
            endmembers_x=dict(self.endmembers_x),
            endmembers_w=dict(self.endmembers_w),
            elements_x=dict(self.elements_x),
            elements_w=dict(self.elements_w),
            G=self.G,
            H=self.H,
            S=self.S,
            Cp=self.Cp,
            partial_gibbs=self.partial_gibbs,
            partial_enthalpy=self.partial_enthalpy,
            partial_entropy=self.partial_entropy,
            activity=self.activity,
            partial_heat_capacity=self.partial_heat_capacity,
        )


class SlnPropertyResult:
    """Container for one or more solution-property calculation points."""

    def __init__(
        self,
        points: SlnPropertyPoint | list[SlnPropertyPoint] | None = None,
        *,
        context: ResultContext | None = None,
    ) -> None:
        self.context = context
        self.data = points

    @property
    def data(self) -> SlnPropertyPoint | list[SlnPropertyPoint] | None:
        """Return raw point storage."""
        return self._data

    @data.setter
    def data(
        self,
        value: SlnPropertyPoint | list[SlnPropertyPoint] | None,
    ) -> None:
        self._data = list(value) if isinstance(value, list) else value

    @property
    def points(self) -> list[SlnPropertyPoint]:
        """Return calculation points as a list."""
        if self.data is None:
            return []
        if isinstance(self.data, SlnPropertyPoint):
            return [self.data]
        return list(self.data)

    @property
    def point(self) -> SlnPropertyPoint:
        """Return the only point, or raise for empty/batch results."""
        points = self.points
        if len(points) != 1:
            raise ValueError(
                f"Expected exactly one solution-property point; found {len(points)}."
            )
        return points[0]

    @property
    def T(self) -> float | np.ndarray | None:
        """Return temperature as a scalar or vector."""
        return self._scalar_or_vector("T")

    @property
    def P(self) -> float | np.ndarray | None:
        """Return pressure as a scalar or vector."""
        return self._scalar_or_vector("P")

    @property
    def G(self) -> float | np.ndarray | None:
        """Return Gibbs energy as a scalar or vector."""
        return self._scalar_or_vector("G")

    @property
    def H(self) -> float | np.ndarray | None:
        """Return enthalpy as a scalar or vector."""
        return self._scalar_or_vector("H")

    @property
    def S(self) -> float | np.ndarray | None:
        """Return entropy as a scalar or vector."""
        return self._scalar_or_vector("S")

    @property
    def Cp(self) -> float | np.ndarray | None:
        """Return heat capacity as a scalar or vector."""
        return self._scalar_or_vector("Cp")

    @property
    def n_i(self) -> dict[str, float] | list[dict[str, float]]:
        """Return input component amounts on a mole basis."""
        return self._mapping_or_list("n_i")

    @property
    def w_i(self) -> dict[str, float] | list[dict[str, float]]:
        """Return input component amounts on a mass basis."""
        return self._mapping_or_list("w_i")

    @property
    def endmembers_x(self) -> dict[str, float] | list[dict[str, float]]:
        """Return endmember mole fractions."""
        return self._mapping_or_list("endmembers_x")

    @property
    def elements_x(self) -> dict[str, float] | list[dict[str, float]]:
        """Return element mole fractions."""
        return self._mapping_or_list("elements_x")

    @property
    def partial_gibbs(self) -> SpeciesProperty:
        """Return partial molar Gibbs energies by endmember."""
        return self._species_property("partial_gibbs")

    @property
    def partial_enthalpy(self) -> SpeciesProperty:
        """Return partial molar enthalpies by endmember."""
        return self._species_property("partial_enthalpy")

    @property
    def partial_entropy(self) -> SpeciesProperty:
        """Return partial molar entropies by endmember."""
        return self._species_property("partial_entropy")

    @property
    def activity(self) -> SpeciesProperty:
        """Return endmember activities."""
        return self._species_property("activity")

    @property
    def partial_heat_capacity(self) -> SpeciesProperty:
        """Return partial molar heat capacities by endmember."""
        return self._species_property("partial_heat_capacity")

    @property
    def dg_nn(self) -> np.ndarray | None:
        """Return Hessian data when it is explicitly evaluated."""
        points = self.points
        if not points:
            return None
        values = [point.dg_nn for point in points]
        if any(value is None for value in values):
            return None
        if len(values) == 1:
            return np.asarray(values[0], dtype=float).copy()
        return np.stack([np.asarray(value, dtype=float) for value in values])

    @property
    def phases(self) -> PhaseCollection:
        """Return a phase collection for compatibility with result displays."""
        points = self.points
        if len(points) == 1:
            point = points[0]
            return PhaseCollection({point.phase_name: point.phase(point.phase_name)})
        return PhaseCollection(
            {
                point.phase_name: point.phase(point.phase_name)
                for point in points
            }
        )

    @property
    def stable_phases(self) -> PhaseCollection:
        """Return the evaluated solution phase as stable."""
        return self.phases

    def phase(self, name: str) -> SlnPropertyPhase:
        """Return the evaluated phase for a single-point result."""
        return self.point.phase(name)

    def to_dict(self) -> dict[str, Any]:
        """Return a compact flat dictionary for table/export use."""
        points = self.points
        if not points:
            return {}
        data: dict[str, Any] = {
            "T [K]": self.T,
            "P [atm]": self.P,
            "G [J]": self.G,
            "H [J]": self.H,
            "S [J/K]": self.S,
            "Cp [J/K]": self.Cp,
        }
        for prefix, prop in (
            ("partial_gibbs", self.partial_gibbs),
            ("partial_enthalpy", self.partial_enthalpy),
            ("partial_entropy", self.partial_entropy),
            ("activity", self.activity),
            ("partial_heat_capacity", self.partial_heat_capacity),
        ):
            unit = {
                "partial_gibbs": "J",
                "partial_enthalpy": "J",
                "partial_entropy": "J/K",
                "activity": "-",
                "partial_heat_capacity": "J/K",
            }[prefix]
            for name in prop:
                data[f"{prefix}_{name} [{unit}]"] = prop[name]
        return data

    def _scalar_or_vector(self, attribute_name: str) -> float | np.ndarray | None:
        points = self.points
        if not points:
            return None
        values = np.asarray(
            [getattr(point, attribute_name) for point in points],
            dtype=float,
        )
        if len(values) == 1:
            return float(values[0])
        return values

    def _mapping_or_list(
        self,
        attribute_name: str,
    ) -> dict[str, float] | list[dict[str, float]]:
        points = self.points
        if not points:
            return {}
        values = [dict(getattr(point, attribute_name)) for point in points]
        return values[0] if len(values) == 1 else values

    def _species_property(self, attribute_name: str) -> SpeciesProperty:
        points = self.points
        if not points:
            return SpeciesProperty([], np.empty((0,), dtype=float))
        names: list[str] = []
        seen: set[str] = set()
        for point in points:
            prop = getattr(point, attribute_name)
            for name in prop:
                if name not in seen:
                    seen.add(name)
                    names.append(name)
        rows = []
        for point in points:
            prop = getattr(point, attribute_name)
            rows.append([prop[name] if name in prop else np.nan for name in names])
        values = np.asarray(rows, dtype=float)
        if len(points) == 1:
            values = values[0]
        return SpeciesProperty(names, values)
