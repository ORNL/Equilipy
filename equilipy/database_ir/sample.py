"""Small sample DatabaseIR used before parsers are implemented."""

from __future__ import annotations

from .model import (
    ConstituentSet,
    DatabaseIR,
    Diagnostic,
    Element,
    FunctionDefinition,
    Parameter,
    Phase,
    SourceRef,
    Species,
)


def sample_database() -> DatabaseIR:
    """Return a tiny Al-Cu database for GUI and model smoke tests."""
    source = SourceRef(file="mock_alcu.tdb", line=1, column=1, command="MOCK")
    al2cu_g_source = SourceRef(
        file="mock_alcu.tdb",
        line=15,
        column=1,
        command=(
            "PARAMETER G(AL2CU(s),AL2CU;0) 298.15 "
            "2*GHSERAL+GHSERCU-12000+0.5*T; 6000 N !"
        ),
    )
    al2cu_cp_source = SourceRef(
        file="mock_alcu.tdb",
        line=16,
        column=1,
        command="PARAMETER CP(AL2CU(s),AL2CU;0) 298.15 65.0+0.012*T; 6000 N !",
    )
    return DatabaseIR(
        name="Mock Al-Cu Database",
        metadata={
            "source": "mock",
            "purpose": "PySide6 DatabaseIR GUI draft",
        },
        elements=[
            Element("Al", 26.9815384, "FCC_A1", source),
            Element("Cu", 63.546, "FCC_A1", source),
        ],
        species=[
            Species("Al", {"Al": 1.0}, source=source),
            Species("Cu", {"Cu": 1.0}, source=source),
            Species("AL2CU", {"Al": 2.0, "Cu": 1.0}, source=source),
        ],
        functions=[
            FunctionDefinition(
                "GHSERAL",
                "298.15  -7976.15+137.0715*T-24.3671976*T*LN(T); 700 N",
                [(298.15, 700.0)],
                source,
            ),
            FunctionDefinition(
                "GHSERCU",
                "298.15  -7770.458+130.485235*T-24.112392*T*LN(T); 1358 N",
                [(298.15, 1358.0)],
                source,
            ),
        ],
        phases=[
            Phase(
                "LIQUID",
                "IDMX",
                [ConstituentSet(1, ["Al", "Cu"], 1.0)],
                source,
            ),
            Phase(
                "FCC_A1",
                "RKMP",
                [ConstituentSet(1, ["Al", "Cu"], 1.0)],
                source,
            ),
            Phase(
                "AL2CU(s)",
                "COMP",
                [ConstituentSet(1, ["AL2CU"], 1.0)],
                source,
            ),
        ],
        parameters=[
            Parameter(
                "FCC_A1",
                "G",
                ["Al", "Cu"],
                0,
                "-53520+2*T",
                source,
            ),
            Parameter(
                "AL2CU(s)",
                "G",
                ["AL2CU"],
                0,
                "2*GHSERAL+GHSERCU-12000+0.5*T",
                al2cu_g_source,
            ),
            Parameter(
                "AL2CU(s)",
                "Cp",
                ["AL2CU"],
                0,
                "65.0+0.012*T",
                al2cu_cp_source,
            ),
        ],
        diagnostics=[
            Diagnostic(
                "info",
                "Mock DatabaseIR loaded. Parser integration is not active yet.",
                source,
            )
        ],
    )
