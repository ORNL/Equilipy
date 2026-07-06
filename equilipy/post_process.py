"""Post-processing public exports.

Canonical result model implementations live under :mod:`equilipy.results`.
"""

from __future__ import annotations

import equilipy.equilifort as fort  # noqa: F401
import equilipy.variables as var  # noqa: F401

from .exceptions import EquilibError as EquilibError
from .exceptions import PostProcessError as PostProcessError
from .results import ResultTable as ResultTable
from .results.common import combine_dicts as combine_dicts
from .results.common import (
    pad_missing_cumulative_phase_amounts as _pad_missing_cumulative_phase_amounts,
)
from .results.common import phase_is_liquid as _phase_is_liquid
from .results.context import ResultContext as ResultContext
from .results.context import capture_result_context as capture_result_context
from .results.equilib import EquilibPoint as EquilibPoint
from .results.equilib import EquilibResult as EquilibResult
from .results.equilib import PhaseResult as PhaseResult
from .results.equilib import create_phase_from_sys as create_phase_from_sys
from .results.equilib import endmembers2elements as endmembers2elements
from .results.scheil import ScheilConstituentSegment as ScheilConstituentSegment
from .results.scheil import ScheilPoint as ScheilPoint
from .results.scheil import ScheilResult as ScheilResult

__all__ = [
    "EquilibError",
    "EquilibPoint",
    "EquilibResult",
    "PhaseResult",
    "PostProcessError",
    "ResultContext",
    "ResultTable",
    "ScheilConstituentSegment",
    "ScheilPoint",
    "ScheilResult",
    "_pad_missing_cumulative_phase_amounts",
    "_phase_is_liquid",
    "capture_result_context",
    "combine_dicts",
    "create_phase_from_sys",
    "endmembers2elements",
    "fort",
    "var",
]
