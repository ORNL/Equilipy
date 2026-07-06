"""Temporary compatibility re-exports for moved GUI modules."""

from __future__ import annotations

# ruff: noqa: F401,F403,I001

from .calculation import (
    input_condition as _input_condition_module,
    results as _results_module,
    session_serialization as _session_serialization_module,
)
from .calculation.input_condition import *
from .calculation.session_serialization import *
from .calculation.results import *
from .database import (
    cp_editor as _cp_editor_module,
    overview as _database_overview_module,
)
from .database.cp_editor import *
from .database.overview import *

__all__ = [
    name
    for name in globals()
    if name.startswith("_") and not name.startswith("__")
]

_compat_modules = (
    _input_condition_module,
    _database_overview_module,
    _session_serialization_module,
    _results_module,
    _cp_editor_module,
)
_compat_exports = {name: globals()[name] for name in __all__}
for _compat_module in _compat_modules:
    _compat_module.__dict__.update(_compat_exports)
