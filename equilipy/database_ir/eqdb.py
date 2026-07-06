"""Equilipy DatabaseIR JSON save/load helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .model import DatabaseIR


def dumps_eqdb(database: DatabaseIR) -> str:
    """Serialize a DatabaseIR as readable Equilipy database JSON."""
    return json.dumps(_clean_eqdb_payload(database.to_dict()), indent=2) + "\n"


def write_eqdb(database: DatabaseIR, path: str | Path) -> Path:
    """Write a DatabaseIR to an ``.eqdb`` JSON file."""
    output_path = Path(path)
    output_path.write_text(dumps_eqdb(database), encoding="utf-8")
    return output_path


def _clean_eqdb_payload(value: Any, *, parent_key: str = "") -> Any:
    """Remove noisy source-path details from the saved native DatabaseIR JSON."""
    if isinstance(value, dict):
        if _looks_like_function_payload(value):
            value = _clean_function_payload(value)
        cleaned: dict[str, Any] = {}
        for key, item in value.items():
            if parent_key == "metadata" and key == "_loaded_path":
                continue
            if parent_key == "metadata" and key == "source_file":
                cleaned[key] = Path(str(item)).name
                continue
            if key == "source":
                continue
            cleaned[key] = _clean_eqdb_payload(item, parent_key=key)
        return cleaned
    if isinstance(value, list):
        return [_clean_eqdb_payload(item) for item in value]
    return value


def _looks_like_function_payload(value: dict[str, Any]) -> bool:
    """Return whether a dict is a serialized FunctionDefinition."""
    return "gibbs_ranges" in value and "tdb_src" in value


def _clean_function_payload(value: dict[str, Any]) -> dict[str, Any]:
    """Keep the native function payload non-duplicative."""
    cleaned = dict(value)
    expression = str(cleaned.pop("expression", ""))
    cleaned.pop("temperature_ranges", None)
    tdb_src = dict(cleaned.get("tdb_src") or {})
    if expression and not tdb_src.get("source_expression"):
        tdb_src["source_expression"] = expression
    cleaned["tdb_src"] = tdb_src
    return cleaned
