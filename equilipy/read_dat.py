"""Read and parse ChemSage DAT database files."""

from __future__ import annotations

import re
from pathlib import Path

import equilipy.variables as var

from .exceptions import DatabaseParsingError
from .parse_chemsage_data_block import ParseCSDataBlock
from .parse_chemsage_header import ParseCSHeader
from .parse_HSCp_functions import ParseHSCpFunctions

_LEGACY_DATABASE_NAME_PATTERNS = (
    re.compile(r"(?<![a-z0-9])fs[\s_.-]*7(?:[\s_.-]*[0-9x]+)?(?=[^a-z0-9]|$)"),
    re.compile(
        r"fact[\s_.-]*sage[\s_.-]*(?:v(?:ersion)?[\s_.-]*)?"
        r"7(?:[\s_.-]*[0-9x]+)?(?=[^a-z0-9]|$)"
    ),
    re.compile(r"(?:chem[\s_.-]*sage|cs)[\s_.-]*(?:old|legacy|classic)(?=[^a-z0-9]|$)"),
    re.compile(
        r"(?<![a-z0-9])(?:old|legacy|classic)[\s_.-]*"
        r"(?:chem[\s_.-]*sage|cs)(?=[^a-z0-9]|$)"
    ),
    re.compile(
        r"(?<![a-z0-9])pre[\s_.-]*(?:fs|fact[\s_.-]*sage)"
        r"[\s_.-]*8(?=[^a-z0-9]|$)"
    ),
)


def _infer_factsage_8_plus(file_name: str) -> bool:
    """Infer whether a DAT file should use the FactSage 8+ parser mode."""
    normalized_name = Path(file_name).name.lower()
    return not any(
        pattern.search(normalized_name)
        for pattern in _LEGACY_DATABASE_NAME_PATTERNS
    )


def read_dat(file_name, factsage_8_plus=None):
    """Read a DAT database file and return a parsed database dictionary."""
    assert isinstance(file_name, str), "Error: File name must be a string"
    if factsage_8_plus is None:
        factsage_8_plus = _infer_factsage_8_plus(file_name)

    var.cSolnPhaseTypeSupport = [
        "IDMX    ",
        "QKTO    ",
        "RKMP    ",
        "RKMPM   ",
        "QKTOM   ",  # Parsing of this model has not been tested
        "SUBLM   ",
        "SUBOM   ",
        "SUBL    ",
        "SUBG    ",
        "SUBQ    ",
    ]

    var.FactSage8Plus = factsage_8_plus

    # Initialize variables:
    infothermo = 0
    var.INFO = 0

    # Attempt to open datafile
    try:
        datafile = open(file_name, "rt")
        lines = datafile.readlines()
        var.DataBase = []
        for i in range(1, len(lines)):
            var.DataBase = var.DataBase + lines[i].split()
    except IOError as e:
        raise DatabaseParsingError(
            "Database file not found or the path is incorrect"
        ) from e

    if factsage_8_plus:
        # Process function header
        ParseHSCpFunctions()

    # Parse the "header section" of the data-file:
    if infothermo == 0:
        ParseCSHeader()
    if var.INFO != 0:
        infothermo = var.INFO

    # Parse the "data block section" of the data-file:
    if infothermo == 0:
        ParseCSDataBlock()
    if var.INFO != 0:
        infothermo = var.INFO

    # All phase names
    var.cPhaseNames = var.cSolnPhaseNameCS + var.cSpeciesNameCS[-var.nPureSpeciesCS :]
    var.cPhaseNames = [x.strip() for x in var.cPhaseNames]
    if "" in var.cPhaseNames:
        var.cPhaseNames.remove("")

    return var.to_dict()
