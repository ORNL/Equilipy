"""Scheil-Gulliver result models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Union

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

import equilipy.equilifort as fort
from equilipy.exceptions import EquilibError, PostProcessError

from .common import (
    apply_temperature_axis_unit as _apply_temperature_axis_unit,
)
from .common import (
    as_list_for_segments as _as_list_for_segments,
)
from .common import (
    clean_quantity_mapping as _clean_quantity_mapping,
)
from .common import (
    component_names_from_amount_maps as _component_names_from_amount_maps,
)
from .common import (
    pad_missing_cumulative_phase_amounts as _pad_missing_cumulative_phase_amounts,
)
from .common import (
    phase_is_liquid as _phase_is_liquid,
)
from .common import (
    temperature_from_kelvin as _temperature_from_kelvin,
)
from .common import (
    temperature_unit_label as _temperature_unit_label,
)
from .context import ResultContext
from .equilib import EquilibResult
from .serialization import (
    SCHEIL_BUNDLE_KINDS,
)
from .serialization import (
    context_from_payload as _context_from_payload,
)
from .serialization import (
    context_to_payload as _context_to_payload,
)
from .serialization import (
    json_safe_mapping as _json_safe_mapping,
)
from .serialization import (
    scheil_point_from_payload as _scheil_row_from_payload,
)
from .serialization import (
    scheil_point_to_payload as _scheil_row_to_payload,
)
from .serialization import (
    validate_result_bundle as _validate_result_bundle,
)
from .table import ResultTable

_PHASE_ENDMEMBER_EXPORT_MARKERS = (
    "_endmembers_x_",
    "_endmembers_w_",
    "_partial_gibbs_",
    "_standard_gibbs_energy_",
    "_activity_",
    "_partial_enthalpy_",
    "_partial_entropy_",
    "_partial_heat_capacity_",
)

__all__ = [
    "ScheilConstituentSegment",
    "ScheilPoint",
    "ScheilResult",
]


@dataclass(frozen=True)
class ScheilConstituentSegment:
    """Contiguous Scheil path segment with the same stable phase label."""

    label: str
    x: List[Any]
    y: List[Any]


class ScheilPoint(BaseModel):
    """Single step in a Scheil-Gulliver simulation."""

    T: float
    P: float
    G: Optional[float]
    H: Optional[float]
    S: Optional[float]
    Cp: Optional[float]
    fl: float = Field(ge=0.0)
    fs: float = Field(ge=0.0)
    fl_w: float = Field(ge=0.0)
    fs_w: float = Field(ge=0.0)
    label: str
    cumulative_phases_n: Dict[str, float] = Field(default_factory=dict)
    cumulative_phases_w: Dict[str, float] = Field(default_factory=dict)

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        extra="forbid",
        populate_by_name=True,
    )


class ScheilResult:
    """Scheil-Gulliver cooling result object."""

    def __init__(
        self,
        liquid_phase_name: str | None = None,
        input_unit: Optional[List[str]] = None,
        context: Optional[ResultContext] = None,
    ):
        self.context = context
        self.equilib_result = EquilibResult(context=context)
        self.liquid_phase_name = liquid_phase_name
        self.input_unit = input_unit.copy() if input_unit else ["K", "atm", "moles"]
        self.data: Union[None, ScheilPoint, List[ScheilPoint]] = None
        self._n_i: Dict[str, float] = {}
        self._w_i: Dict[str, float] = {}
        self._scheil_constituents: Dict[str, float] = {}
        self.warnings: List[str] = []

        try:
            self.equilib_result.append_output()
            self.n_i = self.equilib_result.n_i
            self.w_i = self.equilib_result.w_i
            self.data = self._scheil_point_from_equilib_step(
                self.equilib_result.point,
                previous_point=None,
            )
        except Exception as exc:
            print(f"Error during ScheilResult initialization: {exc}")
            raise
        finally:
            fort.resetthermo()

    @property
    def points(self) -> List[ScheilPoint]:
        """Return Scheil points as a list regardless of internal storage."""
        if self.data is None:
            return []
        if isinstance(self.data, ScheilPoint):
            return [self.data]
        return list(self.data)

    @property
    def point(self) -> ScheilPoint:
        """Return the only Scheil point, or raise for empty/multi-step results."""
        points = self.points
        if len(points) != 1:
            raise PostProcessError(
                "Expected exactly one Scheil point; "
                f"found {len(points)}."
            )
        return points[0]

    @property
    def n_i(self) -> Dict[str, float]:
        """Return initial system component amounts on a mole basis."""
        return dict(getattr(self, "_n_i", {}))

    @n_i.setter
    def n_i(self, values: Mapping[str, Any]) -> None:
        self._n_i = _clean_quantity_mapping(values)

    @property
    def w_i(self) -> Dict[str, float]:
        """Return initial system component amounts on a mass basis."""
        return dict(getattr(self, "_w_i", {}))

    @w_i.setter
    def w_i(self, values: Mapping[str, Any]) -> None:
        self._w_i = _clean_quantity_mapping(values)

    @property
    def Q(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return enthalpy change relative to the first Scheil point."""
        enthalpy = self.H
        if enthalpy is None:
            return None
        if isinstance(enthalpy, list):
            if not enthalpy:
                return []
            h0 = enthalpy[0]
            return [
                None if h0 is None or value is None else value - h0
                for value in enthalpy
            ]
        return 0.0

    @property
    def fl(self) -> Union[float, List[float], None]:
        """Return liquid fraction on a mole basis."""
        if isinstance(self.data, ScheilPoint):
            return self.data.fl
        if isinstance(self.data, list):
            return [row.fl for row in self.data]
        return None

    @property
    def fs(self) -> Union[float, List[float], None]:
        """Return solid fraction on a mole basis."""
        if isinstance(self.data, ScheilPoint):
            return self.data.fs
        if isinstance(self.data, list):
            return [row.fs for row in self.data]
        return None

    @property
    def fl_w(self) -> Union[float, List[float], None]:
        """Return liquid fraction on a mass basis."""
        if isinstance(self.data, ScheilPoint):
            return self.data.fl_w
        if isinstance(self.data, list):
            return [row.fl_w for row in self.data]
        return None

    @property
    def fs_w(self) -> Union[float, List[float], None]:
        """Return solid fraction on a mass basis."""
        if isinstance(self.data, ScheilPoint):
            return self.data.fs_w
        if isinstance(self.data, list):
            return [row.fs_w for row in self.data]
        return None

    @property
    def phase_amounts_n(self) -> Dict[str, List[float]]:
        """Return cumulative phase amount histories on a mole basis."""
        points = self.points
        phase_names = set()
        for row in points:
            phase_names.update(row.cumulative_phases_n)
        return {
            phase_name: [
                row.cumulative_phases_n.get(phase_name, 0.0)
                for row in points
            ]
            for phase_name in phase_names
        }

    @property
    def phase_amounts_w(self) -> Dict[str, List[float]]:
        """Return cumulative phase amount histories on a mass basis."""
        points = self.points
        phase_names = set()
        for row in points:
            phase_names.update(row.cumulative_phases_w)
        return {
            phase_name: [
                row.cumulative_phases_w.get(phase_name, 0.0)
                for row in points
            ]
            for phase_name in phase_names
        }

    def phase_amount_n(self, name: str) -> List[float]:
        """Return one cumulative phase amount history on a mole basis."""
        return self.phase_amounts_n[name]

    def phase_amount_w(self, name: str) -> List[float]:
        """Return one cumulative phase amount history on a mass basis."""
        return self.phase_amounts_w[name]

    @property
    def phase_labels(self) -> Union[str, List[str]]:
        """Return stable phase labels along the Scheil path."""
        if isinstance(self.data, ScheilPoint):
            return self.data.label
        if isinstance(self.data, list):
            return [row.label for row in self.data]
        return []

    def scheil_constituent_segments(
        self,
        x: str = "fs_w",
        y: str = "T",
        *,
        connect: bool = True,
        x_unit: Optional[str] = None,
        y_unit: Optional[str] = None,
    ) -> List[ScheilConstituentSegment]:
        """
        Return Scheil constituent segments grouped by contiguous stable labels.

        ``x`` and ``y`` are names of per-step ScheilResult properties such as
        ``"T"``, ``"Q"``, ``"fs"``, ``"fs_w"``, ``"fl"``, ``"fl_w"``, ``"G"``,
        ``"H"``, ``"S"``, and ``"Cp"``.  Use ``x_unit`` or ``y_unit`` only when
        that axis is ``"T"``; accepted values are ``"K"``, ``"C"``, ``"F"``,
        ``"R"``, or ``"input"``.  When ``connect`` is true, each new segment
        starts from the previous point so phase-change lines remain visually
        connected.
        """
        labels = _as_list_for_segments(self.phase_labels)
        try:
            x_values = _as_list_for_segments(getattr(self, x))
            y_values = _as_list_for_segments(getattr(self, y))
        except AttributeError as exc:
            raise PostProcessError(
                "scheil_constituent_segments x and y must be names of "
                f"ScheilResult properties. Received x={x!r}, y={y!r}."
            ) from exc

        input_unit = getattr(self, "input_unit", ["K", "atm", "moles"])
        x_values = _apply_temperature_axis_unit(
            x,
            x_values,
            x_unit,
            input_unit,
            "x_unit",
        )
        y_values = _apply_temperature_axis_unit(
            y,
            y_values,
            y_unit,
            input_unit,
            "y_unit",
        )

        if not labels:
            return []
        if len(labels) != len(x_values) or len(labels) != len(y_values):
            raise PostProcessError(
                "Cannot build phase-label segments because selected result "
                f"series have different lengths: labels={len(labels)}, "
                f"{x}={len(x_values)}, {y}={len(y_values)}."
            )

        segments: List[ScheilConstituentSegment] = []
        current_label = labels[0]
        current_x = [x_values[0]]
        current_y = [y_values[0]]

        for index in range(1, len(labels)):
            label = labels[index]
            if label != current_label:
                segments.append(
                    ScheilConstituentSegment(
                        label=str(current_label),
                        x=current_x,
                        y=current_y,
                    )
                )
                if connect:
                    current_x = [x_values[index - 1], x_values[index]]
                    current_y = [y_values[index - 1], y_values[index]]
                else:
                    current_x = [x_values[index]]
                    current_y = [y_values[index]]
                current_label = label
                continue

            current_x.append(x_values[index])
            current_y.append(y_values[index])

        segments.append(
            ScheilConstituentSegment(
                label=str(current_label),
                x=current_x,
                y=current_y,
            )
        )
        return segments

    def append_output(self) -> None:
        """Append the next Scheil step's results from the current Fortran state."""
        try:
            self.equilib_result.append_output()
            points = self.points
            if not points:
                raise EquilibError(
                    "Cannot append_output, ScheilResult is not initialized."
                )
            new_row = self._scheil_point_from_equilib_step(
                self.equilib_result.points[-1],
                previous_point=points[-1],
            )
            if isinstance(self.data, ScheilPoint):
                self.data = [self.data, new_row]
            elif isinstance(self.data, list):
                self.data.append(new_row)
            else:
                self.data = new_row
        except Exception as exc:
            print(f"Error during ScheilResult.append_output: {exc}")
            raise
        finally:
            fort.resetthermo()

    def _scheil_point_from_equilib_step(
        self,
        step: Any,
        *,
        previous_point: ScheilPoint | None,
    ) -> ScheilPoint:
        """Build a Scheil point from one equilibrium point."""
        stable_summary = step.stable_phase_summary
        phase_names = stable_summary.get("name", np.array([]))
        amounts_n = stable_summary.get("amount_n", np.array([]))
        amounts_w = stable_summary.get("amount_w", np.array([]))

        current_fl_n = 0.0
        current_fs_n = 0.0
        current_fl_w = 0.0
        current_fs_w = 0.0
        label_phases: List[str] = []
        cumulative_n: Dict[str, float] = {}
        cumulative_w: Dict[str, float] = {}

        last_cumulative_n = (
            previous_point.cumulative_phases_n if previous_point is not None else {}
        )
        last_cumulative_w = (
            previous_point.cumulative_phases_w if previous_point is not None else {}
        )

        for index, name in enumerate(phase_names):
            if name == " " or name == "nan":
                continue

            amount_n = float(np.round(amounts_n[index], 10))
            amount_w = float(np.round(amounts_w[index], 10))
            label_phases.append(name)

            if _phase_is_liquid(name, self.liquid_phase_name):
                cumulative_n[name] = amount_n
                cumulative_w[name] = amount_w
                current_fl_n += amount_n
                current_fl_w += amount_w
                continue

            previous_n = last_cumulative_n.get(name, 0.0)
            previous_w = last_cumulative_w.get(name, 0.0)
            current_fs_n += amount_n + previous_n
            current_fs_w += amount_w + previous_w
            cumulative_n[name] = previous_n + amount_n
            cumulative_w[name] = previous_w + amount_w

        if previous_point is not None:
            for name in last_cumulative_n:
                if name in phase_names:
                    continue
                previous_n = last_cumulative_n.get(name, 0.0)
                previous_w = last_cumulative_w.get(name, 0.0)
                current_fs_n += previous_n
                current_fs_w += previous_w
                if _phase_is_liquid(name, self.liquid_phase_name):
                    cumulative_n[name] = 0.0
                    cumulative_w[name] = 0.0
                else:
                    cumulative_n[name] = previous_n
                    cumulative_w[name] = previous_w

        total_n = current_fl_n + current_fs_n
        total_w = current_fl_w + current_fs_w
        fl = min(current_fl_n / total_n, 1.0) if total_n > 0 else 0.0
        fl_w = min(current_fl_w / total_w, 1.0) if total_w > 0 else 0.0

        if not label_phases:
            raise EquilibError("No phase appears stable during Scheil simulation")

        if previous_point is not None:
            _pad_missing_cumulative_phase_amounts(
                last_cumulative_n,
                last_cumulative_w,
                cumulative_n,
                cumulative_w,
                self.liquid_phase_name,
            )

        return ScheilPoint(
            T=step.T,
            P=step.P,
            G=step.G,
            H=step.H,
            S=step.S,
            Cp=step.Cp,
            fl=fl,
            fs=1.0 - fl,
            fl_w=fl_w,
            fs_w=1.0 - fl_w,
            label="+".join(sorted(label_phases)),
            cumulative_phases_n=cumulative_n,
            cumulative_phases_w=cumulative_w,
        )

    @property
    def T(self) -> Union[float, List[float], None]:
        """Return Scheil temperatures."""
        if isinstance(self.data, ScheilPoint):
            return self.data.T
        if isinstance(self.data, list):
            return [row.T for row in self.data]
        return None

    @property
    def P(self) -> Union[float, List[float], None]:
        """Return Scheil pressures."""
        if isinstance(self.data, ScheilPoint):
            return self.data.P
        if isinstance(self.data, list):
            return [row.P for row in self.data]
        return None

    @property
    def G(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return Scheil Gibbs energy values."""
        if isinstance(self.data, ScheilPoint):
            return self.data.G
        if isinstance(self.data, list):
            return [row.G for row in self.data]
        return None

    @property
    def H(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return Scheil enthalpy values."""
        if isinstance(self.data, ScheilPoint):
            return self.data.H
        if isinstance(self.data, list):
            return [row.H for row in self.data]
        return None

    @property
    def S(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return Scheil entropy values."""
        if isinstance(self.data, ScheilPoint):
            return self.data.S
        if isinstance(self.data, list):
            return [row.S for row in self.data]
        return None

    @property
    def Cp(self) -> Union[Optional[float], List[Optional[float]]]:
        """Return Scheil heat capacity values."""
        if isinstance(self.data, ScheilPoint):
            return self.data.Cp
        if isinstance(self.data, list):
            return [row.Cp for row in self.data]
        return None

    @property
    def scheil_constituents(self) -> Dict[str, float]:
        """Return Scheil constituent summary values."""
        return dict(self._scheil_constituents)

    @scheil_constituents.setter
    def scheil_constituents(self, values: Dict[str, float]) -> None:
        """Set Scheil constituent summary values."""
        self._scheil_constituents = dict(values)

    def update_scheil_constituents(self) -> Dict[str, float]:
        """Rebuild and return Scheil constituent summary values."""
        self._update_scheil_constituents_from_segments()
        return self.scheil_constituents

    def _update_scheil_constituents_from_segments(self) -> None:
        """Update final Scheil constituent amounts from phase-label segments."""
        self._scheil_constituents = {}
        points = self.points
        if len(points) < 2:
            return None

        temperatures = [row.T for row in points]
        labels = [row.label for row in points]
        fs = [row.fs for row in points]
        fs_w = [row.fs_w for row in points]
        component_names = _component_names_from_amount_maps(
            self.context or self.equilib_result.context,
            self.n_i,
            self.w_i,
        )

        self._scheil_constituents["T_Delta [K]"] = temperatures[0] - temperatures[-1]
        for component in component_names:
            self._scheil_constituents[f"{component} [sp-mol]"] = self.n_i.get(
                component, 0.0
            )
        for component in component_names:
            self._scheil_constituents[f"{component} [g]"] = self.w_i.get(
                component, 0.0
            )

        solid_labels: list[str] = []
        for label in labels:
            phases = [
                phase
                for phase in label.split("+")
                if not _phase_is_liquid(phase, self.liquid_phase_name) and phase != ""
            ]
            solid_labels.append("+".join(phases))

        def add_segment(phase_name: str, start_index: int, end_index: int) -> None:
            if not phase_name:
                return
            base_index = start_index - 1
            base_fs = fs[base_index] if base_index >= 0 else 0.0
            base_fs_w = fs_w[base_index] if base_index >= 0 else 0.0
            self._scheil_constituents[f"{phase_name} [sp-mol]"] = (
                self._scheil_constituents.get(f"{phase_name} [sp-mol]", 0.0)
                + (fs[end_index] - base_fs)
            )
            self._scheil_constituents[f"{phase_name} [g]"] = (
                self._scheil_constituents.get(f"{phase_name} [g]", 0.0)
                + (fs_w[end_index] - base_fs_w)
            )

        start_index = 0
        current_phase = solid_labels[0]
        for index in range(1, len(solid_labels)):
            if solid_labels[index] == current_phase:
                continue
            add_segment(current_phase, start_index, index - 1)
            start_index = index
            current_phase = solid_labels[index]

        add_segment(current_phase, start_index, len(solid_labels) - 1)

        return None

    def to_dict(self) -> Dict[str, Any]:
        """Export this Scheil result to a flattened dictionary."""
        points = self.points
        if not points:
            return {}

        df: Dict[str, Any] = {}
        df_eq = self.equilib_result.to_dict()
        temperature_unit = _temperature_unit_label(self.input_unit[0])
        df[f"T [{temperature_unit}]"] = [
            _temperature_from_kelvin(row.T, temperature_unit) for row in points
        ]
        if temperature_unit != "K":
            df["T [K]"] = [row.T for row in points]
        df["P [atm]"] = [row.P for row in points]

        component_names = _component_names_from_amount_maps(
            self.context or self.equilib_result.context,
            self.n_i,
            self.w_i,
        )
        for component in component_names:
            df[f"n_{component} [sp-mol]"] = [
                self.n_i.get(component, 0.0) for _ in points
            ]
        for component in component_names:
            df[f"w_{component} [g]"] = [
                self.w_i.get(component, 0.0) for _ in points
            ]

        df["label"] = [row.label for row in points]
        df["fl"] = [row.fl for row in points]
        df["fs"] = [row.fs for row in points]
        df["fl_w"] = [row.fl_w for row in points]
        df["fs_w"] = [row.fs_w for row in points]
        df["G [J]"] = [row.G for row in points]
        df["H [J]"] = [row.H for row in points]
        q_values = self.Q
        df["Q [J]"] = q_values if isinstance(q_values, list) else [q_values]
        df["S [J/K]"] = [row.S for row in points]
        df["Cp [J/K]"] = [row.Cp for row in points]

        phase_amounts_n = self.phase_amounts_n
        phase_amounts_w = self.phase_amounts_w
        phase_names = list(phase_amounts_n)
        for phase_name in phase_names:
            df[f"{phase_name}_amount_n [sp-mol]"] = phase_amounts_n[phase_name]
            df[f"{phase_name}_amount_w [g]"] = phase_amounts_w.get(
                phase_name,
                [0.0 for _ in points],
            )
            phase_prefix = f"{phase_name}_"
            for key, value in df_eq.items():
                if key.startswith(phase_prefix) and any(
                    marker in key for marker in _PHASE_ENDMEMBER_EXPORT_MARKERS
                ):
                    df[key] = value

        return df

    def to_table(self) -> ResultTable:
        """Return a selectable/exportable table view of this Scheil result."""
        return ResultTable.from_dict(self.to_dict())

    def available_columns(self) -> List[str]:
        """Return flattened Scheil result columns available for export/display."""
        return self.to_table().available_columns()

    def to_bundle(self) -> Dict[str, Any]:
        """Return a JSON-safe bundle that can reconstruct this Scheil result."""
        points = self.points
        if not points:
            data_state = "none"
        elif isinstance(self.data, ScheilPoint):
            data_state = "single"
        else:
            data_state = "list"

        return {
            "format": "equilipy.result",
            "version": 1,
            "kind": "scheil",
            "context": _context_to_payload(self.context),
            "liquid_phase_name": self.liquid_phase_name,
            "input_unit": list(self.input_unit),
            "n_i": _json_safe_mapping(self.n_i),
            "w_i": _json_safe_mapping(self.w_i),
            "scheil_constituents": _json_safe_mapping(self.scheil_constituents),
            "equilib_result": self.equilib_result.to_bundle(),
            "data_state": data_state,
            "data": [_scheil_row_to_payload(point) for point in points],
        }

    @classmethod
    def from_bundle(cls, bundle: Dict[str, Any]) -> "ScheilResult":
        """Reconstruct a Scheil result without calling the Fortran initializer."""
        _validate_result_bundle(bundle, SCHEIL_BUNDLE_KINDS)
        result = cls.__new__(cls)
        result.context = _context_from_payload(bundle.get("context"))
        equilib_bundle = bundle["equilib_result"]
        result.equilib_result = EquilibResult.from_bundle(equilib_bundle)
        if result.equilib_result.context is None and result.context is not None:
            result.equilib_result.context = result.context
        result.liquid_phase_name = bundle.get("liquid_phase_name")
        result.input_unit = list(bundle.get("input_unit") or ["K", "atm", "moles"])
        result.n_i = bundle.get("n_i") or {}
        result.w_i = bundle.get("w_i") or {}
        result.scheil_constituents = dict(bundle.get("scheil_constituents") or {})

        rows = [
            _scheil_row_from_payload(row)
            for row in bundle.get("data", [])
            if isinstance(row, dict)
        ]
        data_state = bundle.get("data_state")
        if data_state == "single":
            result.data = rows[0] if rows else None
        elif data_state == "list":
            result.data = rows
        elif data_state in {"none", None}:
            result.data = None
        else:
            raise PostProcessError(f"Unsupported Scheil data state: {data_state!r}.")
        return result
