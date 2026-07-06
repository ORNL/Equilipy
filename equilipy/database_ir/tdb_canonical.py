"""Canonical phase names for TDB import/export.

The IR should use these canonical names internally.  The map exists at TDB
boundaries so older/vendor spellings can be read, calculation can treat helper
names as the same physical phase, and the writer can choose either canonical
or source-helper-preserving export dialects.
"""

from __future__ import annotations

DISORDERED_PHASE_CANONICAL_NAMES: dict[str, str] = {
    "A1_FCC": "FCC_A1",
    "A2_BCC": "BCC_A2",
}


def canonical_disordered_phase_name(name: str) -> str:
    """Return the canonical name for a known disordered helper phase."""
    normalized = str(name or "").strip().upper()
    return DISORDERED_PHASE_CANONICAL_NAMES.get(normalized, normalized)
