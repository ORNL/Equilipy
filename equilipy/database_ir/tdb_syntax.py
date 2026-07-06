"""Shared TDB command spelling and syntax helpers."""

from __future__ import annotations

import re

COMMAND_RE = re.compile(r"^\s*([A-Z_-]+)\b\s*(.*)$", re.IGNORECASE | re.DOTALL)

FUNCTION_KEYWORDS = frozenset({"FUNCTION", "FUNCT", "FUNC"})
ELEMENT_KEYWORDS = frozenset({"ELEMENT", "ELEM"})
SPECIES_KEYWORDS = frozenset({"SPECIES", "SPECIE"})
PHASE_KEYWORDS = frozenset({"PHASE"})
CONSTITUENT_KEYWORDS = frozenset({"CONSTITUENT", "CONST"})
PARAMETER_KEYWORDS = frozenset({"PARAMETER", "PARAM", "PAR"})
REFERENCE_KEYWORDS = frozenset(
    {
        "ADD-REFERENCES",
        "ADD_REFERENCES",
        "LIST-OF-REFERENCE",
        "LIST-OF-REFERENCES",
        "LIST_OF_REFERENCE",
        "LIST_OF_REFERENCES",
    }
)
SYSTEM_DEFAULT_KEYWORDS = frozenset(
    {
        "DEFINE-SYSTEM-DEFAULT",
        "DEFINE_SYSTEM_DEFAULT",
    }
)
DEFAULT_COMMAND_KEYWORDS = frozenset(
    {
        "DEFAULT-COM",
        "DEFAULT_COM",
        "DEFAULT-COMMAND",
        "DEFAULT_COMMAND",
    }
)
TYPE_DEFINITION_KEYWORDS = frozenset({"TYPE-DEF", "TYPE_DEF", "TYPE_DEFINITION"})
PRESERVED_COMMENT_COMMANDS = (
    *SYSTEM_DEFAULT_KEYWORDS,
    *DEFAULT_COMMAND_KEYWORDS,
    *TYPE_DEFINITION_KEYWORDS,
)


def command_keyword_pattern(keywords: frozenset[str]) -> str:
    """Return an alternation pattern for TDB command aliases."""
    return "|".join(
        re.escape(keyword)
        for keyword in sorted(keywords, key=len, reverse=True)
    )


_PARAMETER_COMMAND_PREFIX_RE = re.compile(
    rf"(?is)^(\s*(?:{command_keyword_pattern(PARAMETER_KEYWORDS)})"
    r"\s+\S+\(.*;\s*[-+]?\d+\)\s*)"
)


def parameter_command_prefix_match(command: str) -> re.Match[str] | None:
    """Return the TDB PARAMETER command prefix before the temperature expression."""
    return _PARAMETER_COMMAND_PREFIX_RE.match(str(command or "").strip())
