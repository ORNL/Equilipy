"""Serializable result context metadata."""

from __future__ import annotations

import hashlib
from typing import Any, Dict, List, Optional

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from .serialization import json_safe_mapping as _json_safe_mapping


class ResultContext(BaseModel):
    """Serializable metadata needed to post-process a result after reload."""

    database_name: Optional[str] = None
    database_path: Optional[str] = None
    database_hash: Optional[str] = None
    input_condition: Dict[str, Any] = Field(default_factory=dict)
    unit: List[str] = Field(default_factory=lambda: ["K", "atm", "moles"])
    selected_phases: Optional[List[str]] = None
    component_names: List[str] = Field(default_factory=list)
    element_names: List[str] = Field(default_factory=list)
    phase_names: List[str] = Field(default_factory=list)
    solution_phase_names: List[str] = Field(default_factory=list)
    compound_phase_names: List[str] = Field(default_factory=list)
    species_names: List[str] = Field(default_factory=list)
    endmember_names: List[str] = Field(default_factory=list)
    atomic_masses: Dict[str, float] = Field(default_factory=dict)

    model_config = ConfigDict(arbitrary_types_allowed=True)


def capture_result_context(
    database: Any,
    condition: Dict[str, Any],
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
) -> ResultContext:
    """Capture result metadata without requiring live Fortran state later."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    safe_condition = _json_safe_mapping(condition)
    component_names = [key for key in safe_condition if key not in {"T", "P"}]

    return ResultContext(
        database_name=_database_metadata_text(database, "name"),
        database_path=_database_metadata_text(database, "path"),
        database_hash=_database_hash(database),
        input_condition=safe_condition,
        unit=list(unit),
        selected_phases=list(phases) if phases is not None else None,
        component_names=component_names,
        element_names=_database_element_names(database),
        phase_names=_database_phase_names(database),
        solution_phase_names=_database_solution_phase_names(database),
        compound_phase_names=_database_compound_phase_names(database),
        species_names=_database_species_names(database),
        endmember_names=_database_endmember_names(database),
        atomic_masses=_database_atomic_masses(database),
    )


def _database_metadata_text(database: Any, key: str) -> Optional[str]:
    if not isinstance(database, dict):
        return None
    candidates = (
        key,
        key.capitalize(),
        f"database_{key}",
        f"Database{key.capitalize()}",
    )
    for candidate in candidates:
        value = database.get(candidate)
        if value:
            return str(value)
    return None


def _database_hash(database: Any) -> Optional[str]:
    if not isinstance(database, dict):
        return None
    digest = hashlib.sha256()
    for key in sorted(str(item) for item in database):
        value = database.get(key)
        digest.update(key.encode("utf-8", errors="ignore"))
        if hasattr(value, "shape"):
            digest.update(str(getattr(value, "shape", "")).encode())
            digest.update(str(getattr(value, "dtype", "")).encode())
        elif isinstance(value, (list, tuple)):
            digest.update(str(len(value)).encode())
        else:
            digest.update(str(value).encode("utf-8", errors="ignore"))
    return digest.hexdigest()


def _database_element_names(database: Any) -> List[str]:
    if not isinstance(database, dict):
        return []
    return _clean_string_list(database.get("cElementNameCS"))


def _database_solution_phase_names(database: Any) -> List[str]:
    if not isinstance(database, dict):
        return []
    return _clean_string_list(database.get("cSolnPhaseNameCS"))


def _database_compound_phase_names(database: Any) -> List[str]:
    if not isinstance(database, dict):
        return []
    species_names = _database_species_names(database)
    try:
        n_pure_species = int(database.get("nPureSpeciesCS", 0) or 0)
    except (TypeError, ValueError):
        return []
    if n_pure_species <= 0:
        return []
    return species_names[max(0, len(species_names) - n_pure_species) :]


def _database_phase_names(database: Any) -> List[str]:
    if not isinstance(database, dict):
        return []
    names = _clean_string_list(database.get("cPhaseNames"))
    if names:
        return names
    return _unique_preserve_order(
        _database_solution_phase_names(database)
        + _database_compound_phase_names(database)
    )


def _database_species_names(database: Any) -> List[str]:
    if not isinstance(database, dict):
        return []
    return _clean_string_list(database.get("cSpeciesNameCS"))


def _database_endmember_names(database: Any) -> List[str]:
    if not isinstance(database, dict):
        return []
    names = _clean_string_list(database.get("cEndmemberNameCS"))
    return names or _database_species_names(database)


def _database_atomic_masses(database: Any) -> Dict[str, float]:
    if not isinstance(database, dict):
        return {}
    elements = _database_element_names(database)
    masses = database.get("dAtomicMass")
    if masses is None:
        return {}
    atomic_masses: Dict[str, float] = {}
    for element, mass in zip(elements, masses, strict=False):
        try:
            atomic_masses[element] = float(mass)
        except (TypeError, ValueError):
            continue
    return atomic_masses


def _clean_string_list(values: Any) -> List[str]:
    if values is None:
        return []
    if isinstance(values, np.ndarray):
        values = values.tolist()
    cleaned = []
    for value in values:
        text = _text_from_database_value(value)
        if text:
            cleaned.append(text)
    return cleaned


def _text_from_database_value(value: Any) -> str:
    if hasattr(value, "tobytes"):
        try:
            return value.tobytes().decode(errors="ignore").strip()
        except Exception:
            pass
    return str(value).strip()


def _unique_preserve_order(values: List[str]) -> List[str]:
    unique = []
    for value in values:
        if value not in unique:
            unique.append(value)
    return unique
