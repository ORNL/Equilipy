"""Public result-model helpers."""

from __future__ import annotations

from .capture import FortranCaptureState
from .context import ResultContext, capture_result_context
from .phase import (
    CompositionFractions,
    PhaseCollection,
    PhaseSpeciesProperty,
)
from .table import ResultColumn, ResultTable, columns_from_dict

_MODEL_EXPORTS = {
    "EquilibPoint": ("equilipy.results.equilib", "EquilibPoint"),
    "EquilibResult": ("equilipy.results.equilib", "EquilibResult"),
    "PhaseResult": ("equilipy.results.equilib", "PhaseResult"),
    "SlnPropertyPhase": ("equilipy.results.property", "SlnPropertyPhase"),
    "SlnPropertyPoint": ("equilipy.results.property", "SlnPropertyPoint"),
    "SlnPropertyResult": ("equilipy.results.property", "SlnPropertyResult"),
    "SpeciesProperty": ("equilipy.results.property", "SpeciesProperty"),
    "ScheilConstituentSegment": (
        "equilipy.results.scheil",
        "ScheilConstituentSegment",
    ),
    "ScheilPoint": ("equilipy.results.scheil", "ScheilPoint"),
    "ScheilResult": ("equilipy.results.scheil", "ScheilResult"),
}

__all__ = [
    "PhaseCollection",
    "PhaseSpeciesProperty",
    "CompositionFractions",
    "EquilibPoint",
    "EquilibResult",
    "FortranCaptureState",
    "PhaseResult",
    "ResultColumn",
    "ResultTable",
    "ResultContext",
    "SlnPropertyPhase",
    "SlnPropertyPoint",
    "SlnPropertyResult",
    "SpeciesProperty",
    "ScheilConstituentSegment",
    "ScheilPoint",
    "ScheilResult",
    "capture_result_context",
    "columns_from_dict",
]


def __getattr__(name: str):
    """Lazy-load model names without creating post_process import cycles."""
    if name not in _MODEL_EXPORTS:
        raise AttributeError(name)
    module_name, attribute_name = _MODEL_EXPORTS[name]
    module = __import__(module_name, fromlist=[attribute_name])
    return getattr(module, attribute_name)
