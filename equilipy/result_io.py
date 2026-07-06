"""Safe JSON result persistence helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .exceptions import PostProcessError
from .results.equilib import EquilibResult
from .results.scheil import ScheilResult
from .results.serialization import EQUILIB_BUNDLE_KINDS, SCHEIL_BUNDLE_KINDS


def save_result(result: Any, path: str | Path) -> Path:
    """Save an Equilipy result object as a JSON `.eqres` bundle."""
    if not hasattr(result, "to_bundle"):
        raise PostProcessError(
            "Unsupported result object. Expected an Equilipy result with "
            "to_bundle()."
        )

    output_path = Path(path)
    bundle = result.to_bundle()
    output_path.write_text(
        json.dumps(bundle, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return output_path


def load_result(path: str | Path) -> EquilibResult | ScheilResult:
    """Load an Equilipy result object from a JSON `.eqres` bundle."""
    input_path = Path(path)
    try:
        bundle = json.loads(input_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PostProcessError(f"Invalid result bundle JSON: {input_path}") from exc
    except OSError as exc:
        raise PostProcessError(f"Cannot read result bundle: {input_path}") from exc

    if not isinstance(bundle, dict) or bundle.get("format") != "equilipy.result":
        raise PostProcessError("Unsupported result bundle format.")

    kind = bundle.get("kind")
    if kind in EQUILIB_BUNDLE_KINDS:
        return EquilibResult.from_bundle(bundle)
    if kind in SCHEIL_BUNDLE_KINDS:
        return ScheilResult.from_bundle(bundle)
    raise PostProcessError(f"Unsupported result bundle kind: {kind!r}.")
