"""Versioned result bundle serialization helpers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np

from equilipy.exceptions import PostProcessError

EQUILIB_BUNDLE_KINDS = frozenset({"equilib"})
SCHEIL_BUNDLE_KINDS = frozenset({"scheil"})


def json_safe_mapping(values: Mapping[Any, Any]) -> dict[str, Any]:
    """Return a JSON-safe mapping with string keys."""
    return {str(key): json_safe_value(value) for key, value in values.items()}


def json_safe_value(value: Any) -> Any:
    """Return a JSON-safe representation of common scientific Python values."""
    if isinstance(value, np.ndarray):
        return json_safe_value(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Mapping):
        return json_safe_mapping(value)
    if isinstance(value, (list, tuple)):
        return [json_safe_value(item) for item in value]
    return value


def context_to_payload(context: Any | None) -> dict[str, Any] | None:
    """Return a JSON-safe result context payload."""
    if context is None:
        return None
    return context.model_dump(mode="json")


def context_from_payload(payload: Any) -> Any | None:
    """Rebuild a result context from a JSON bundle payload."""
    if payload is None:
        return None
    if not isinstance(payload, dict):
        raise PostProcessError("Result bundle context must be a dictionary.")
    from equilipy.results.context import ResultContext

    return ResultContext(**payload)


def single_point_to_payload(point: Any) -> dict[str, Any]:
    """Return a JSON-safe payload for one Equilib point."""
    return json_safe_value(point.model_dump(mode="python"))


def single_point_from_payload(payload: dict[str, Any]) -> Any:
    """Rebuild one Equilib point from a JSON bundle payload."""
    from equilipy.results.equilib import EquilibPoint, PhaseResult

    payload = dict(payload)
    stable_phases = payload.get("stable_phase_summary", {})
    if isinstance(stable_phases, dict):
        payload["stable_phase_summary"] = {
            key: np.asarray(value) for key, value in stable_phases.items()
        }

    phases = payload.get("phase_map", {})
    if isinstance(phases, dict):
        payload["phase_map"] = {
            name: phase if isinstance(phase, PhaseResult) else PhaseResult(**phase)
            for name, phase in phases.items()
        }

    return EquilibPoint(**payload)


def scheil_point_to_payload(point: Any) -> dict[str, Any]:
    """Return a JSON-safe payload for one Scheil point."""
    return json_safe_value(point.model_dump(mode="python"))


def scheil_point_from_payload(payload: dict[str, Any]) -> Any:
    """Rebuild one Scheil point from a JSON bundle payload."""
    from equilipy.results.scheil import ScheilPoint

    return ScheilPoint(**payload)


def validate_result_bundle(bundle: Any, kinds: str | Iterable[str]) -> None:
    """Validate common result bundle metadata."""
    if isinstance(kinds, str):
        accepted_kinds = frozenset({kinds})
    else:
        accepted_kinds = frozenset(kinds)

    if not isinstance(bundle, dict):
        raise PostProcessError("Result bundle must be a dictionary.")
    if bundle.get("format") != "equilipy.result":
        raise PostProcessError("Unsupported result bundle format.")
    if int(bundle.get("version", 0)) != 1:
        raise PostProcessError("Unsupported result bundle version.")
    if bundle.get("kind") not in accepted_kinds:
        expected = ", ".join(sorted(repr(kind) for kind in accepted_kinds))
        raise PostProcessError(
            f"Expected result bundle kind {expected}, got {bundle.get('kind')!r}."
        )
