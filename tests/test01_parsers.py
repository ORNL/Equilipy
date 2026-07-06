"""Public parser smoke tests for DAT and TDB database readers.

Author:
    S.Y. Kwon
Date:
    06/30/2026
Tests:
    1. test_read_dat_parses_fs73_database_in_legacy_mode: Parse a small FS7.3
       DAT database and confirm legacy mode.
    2. test_read_dat_parses_fs83_database_in_modern_mode: Parse a small FS8.3
       DAT database and confirm modern mode.
    3. test_read_tdb_parses_alfesi_liu_database: Parse AlFeSi_liu.tdb and
       confirm order/disorder phase metadata.
    4. test_tdb_writer_keeps_type_definitions_in_one_setup_block: Keep
       magnetic and DIS_PART commands together on export.
    5. test_tdb_writer_separates_command_terminators: Keep exported command
       terminators separated from the final token.
    6. test_read_dat_infers_factsage_version_from_filename: Infer DAT parser
       mode from representative FactSage/ChemSage filenames.
    7. test_phase_selection_adds_canonical_disordered_partner: Add canonical
       disordered partners for selected ordered phases.
"""

import pytest

import equilipy as eq
from equilipy.database_ir import (
    ConstituentSet,
    DatabaseIR,
    Element,
    FunctionDefinition,
    Parameter,
    Phase,
    TdbCommand,
    dumps_tdb,
    read_tdb,
)
from equilipy.read_dat import _infer_factsage_8_plus


def _phase_names(database: dict) -> set[str]:
    return {str(name).strip() for name in database["cSolnPhaseNameCS"]}


def _element_names(database: dict) -> set[str]:
    return {str(name).strip() for name in database["cElementNameCS"]}


def _solution_name_to_id(database: dict) -> dict[str, int]:
    return {
        str(name).strip(): index
        for index, name in enumerate(database["cSolnPhaseNameCS"], start=1)
    }


def test_read_dat_parses_fs73_database_in_legacy_mode(
    alcumgsi_fs73_database_path,
) -> None:
    """Legacy-mode read_dat parses the FS73 database."""
    database = eq.read_dat(str(alcumgsi_fs73_database_path))

    assert database["INFO"] == 0
    assert database["nElementsCS"] == 4
    assert database["nSolnPhasesSysCS"] == 42
    assert database["nPureSpeciesCS"] == 29
    assert {"Al", "Cu", "Mg", "Si"} == _element_names(database)
    assert {"LIQUID", "FCC_A1", "HCP_A3", "BCC_A2"} <= _phase_names(database)
    assert not eq.variables.FactSage8Plus


def test_read_dat_parses_fs83_database_in_modern_mode(
    alcumgsi_fs83_database_path,
) -> None:
    """Modern-mode read_dat parses the FS83 database."""
    database = eq.read_dat(str(alcumgsi_fs83_database_path))

    assert database["INFO"] == 0
    assert database["nElementsCS"] == 4
    assert database["nSolnPhasesSysCS"] == 40
    assert database["nPureSpeciesCS"] == 25
    assert {"Al", "Cu", "Mg", "Si"} == _element_names(database)
    assert {"LIQUID", "FCC_A1", "HCP_A3", "GAMMA_H"} <= _phase_names(database)
    assert eq.variables.FactSage8Plus


def test_read_tdb_parses_alfesi_liu_database(liu_alfesi_tdb_database_path) -> None:
    """read_tdb parses the Al-Fe-Si 99Liu TDB."""
    database = eq.read_tdb(liu_alfesi_tdb_database_path)
    solution_id = _solution_name_to_id(database)
    disordered_phase = database["iDisorderedPhaseCS"]

    assert database["nElementsCS"] == 3
    assert database["nSolnPhasesSysCS"] == 20
    assert database["nPureSpeciesCS"] == 9
    assert {"Al", "Fe", "Si"} == _element_names(database)
    assert {"FCC_A1", "FCC_4SL", "BCC_A2", "BCC_B2"} <= set(solution_id)
    assert {"A1_FCC", "A2_BCC"} <= set(solution_id)
    assert database["cOrderDisorderHelperPhaseNames"] == ["A1_FCC", "A2_BCC"]
    assert disordered_phase[solution_id["FCC_4SL"] - 1] == solution_id["A1_FCC"]
    assert disordered_phase[solution_id["BCC_B2"] - 1] == solution_id["A2_BCC"]

    canonicalized = eq.read_tdb(
        liu_alfesi_tdb_database_path,
        remove_redundant_phases=True,
    )
    canonical_solution_id = _solution_name_to_id(canonicalized)
    canonical_disordered_phase = canonicalized["iDisorderedPhaseCS"]

    assert canonicalized["nSolnPhasesSysCS"] == 18
    assert "A1_FCC" not in canonical_solution_id
    assert "A2_BCC" not in canonical_solution_id
    assert (
        canonical_disordered_phase[canonical_solution_id["FCC_4SL"] - 1]
        == canonical_solution_id["FCC_A1"]
    )
    assert (
        canonical_disordered_phase[canonical_solution_id["BCC_B2"] - 1]
        == canonical_solution_id["BCC_A2"]
    )


def test_tdb_writer_keeps_type_definitions_in_one_setup_block(
    liu_alfesi_tdb_database_path,
) -> None:
    """The TDB writer keeps TYPE_DEFINITIONs in one setup block."""
    database = read_tdb(
        liu_alfesi_tdb_database_path,
        remove_redundant_phases=False,
    )
    text = dumps_tdb(database, order_disorder_alias_mode="preserve")

    setup_index = text.index(
        "Define SYSTEM DEFAULTS, PHASE TYPE DEFINITIONS, AND ORDRER-DISORDER"
    )
    phase_index = text.index("Define PHASES")
    setup_block = text[setup_index:phase_index]

    assert setup_index < phase_index
    assert "Define ORDER-DISORDER PHASE DEFINITIONS" not in text
    assert (
        "TYPE_DEF C GES AMEND_PHASE_DESCRIPTION BCC_B2 DIS_PART A2_BCC !"
        in setup_block
    )
    assert (
        "TYPE_DEF D GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART A1_FCC !"
        in setup_block
    )


def test_tdb_writer_separates_command_terminators(
    liu_alfesi_tdb_database_path,
) -> None:
    """The TDB writer separates command terminators."""
    database = read_tdb(
        liu_alfesi_tdb_database_path,
        remove_redundant_phases=False,
    )
    text = dumps_tdb(database, order_disorder_alias_mode="preserve")

    assert "DEFAULT_COMMAND DEF_SYS_ELEMENT VA /- !" in text
    assert "CONST A2_BCC : Al Fe Si VA : VA : !" in text
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("CONST "):
            assert "," not in stripped
    for line in text.splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("$") and stripped.endswith("!"):
            assert stripped.endswith(" !")


def _bcc_order_disorder_helper_database(helper_expression: str) -> DatabaseIR:
    return DatabaseIR(
        elements=[
            Element(symbol="AL", atomic_mass=26.98154, reference_state="FCC_A1"),
            Element(symbol="VA", reference_state="VACUUM"),
        ],
        functions=[
            FunctionDefinition(
                name="GBCCAL",
                expression="10083-4.813*T",
                temperature_ranges=[(298.15, 2900.0)],
            ),
        ],
        phases=[
            Phase(
                name="BCC_A2",
                constituents=[
                    ConstituentSet(sublattice=1, species=["AL", "VA"], site_ratio=1),
                    ConstituentSet(sublattice=2, species=["VA"], site_ratio=3),
                ],
            ),
            Phase(
                name="A2_BCC",
                constituents=[
                    ConstituentSet(sublattice=1, species=["AL", "VA"], site_ratio=1),
                    ConstituentSet(sublattice=2, species=["VA"], site_ratio=3),
                ],
            ),
            Phase(
                name="B2_BCC",
                constituents=[
                    ConstituentSet(sublattice=1, species=["AL", "VA"], site_ratio=0.5),
                    ConstituentSet(sublattice=2, species=["AL", "VA"], site_ratio=0.5),
                    ConstituentSet(sublattice=3, species=["VA"], site_ratio=3),
                ],
            ),
        ],
        parameters=[
            Parameter(
                phase_name="BCC_A2",
                parameter_type="G",
                target=["AL:VA"],
                expression="GBCCAL#",
            ),
            Parameter(
                phase_name="A2_BCC",
                parameter_type="G",
                target=["AL:VA"],
                expression=helper_expression,
            ),
        ],
        tdb_commands=[
            TdbCommand(
                kind="type_definition",
                command="TYPE_DEF",
                args=[
                    "O",
                    "GES",
                    "AMEND_PHASE_DESCRIPTION",
                    "B2_BCC",
                    "DIS_PART",
                    "A2_BCC",
                ],
                parsed={
                    "action": "DIS_PART",
                    "ordered_phase": "B2_BCC",
                    "disordered_phase": "A2_BCC",
                },
            ),
        ],
    )


def test_tdb_writer_preserves_existing_order_disorder_helper_parameters() -> None:
    """Existing order/disorder helper parameters survive TDB export."""
    database = _bcc_order_disorder_helper_database("GBCCAL#+0.1")

    text = dumps_tdb(database, order_disorder_alias_mode="preserve")

    assert "GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART A2_BCC !" in text
    assert "PAR G(BCC_A2,Al:VA;0) 298.15 GBCCAL#;" in text
    assert "PAR G(A2_BCC,Al:VA;0) 298.15 0.1+GBCCAL#;" in text


def test_tdb_writer_offsets_degenerate_order_disorder_helper_parameters() -> None:
    """Degenerate order/disorder helper parameters get the +0.1 offset."""
    database = _bcc_order_disorder_helper_database("GBCCAL#")

    text = dumps_tdb(database, order_disorder_alias_mode="preserve")

    assert "PAR G(BCC_A2,Al:VA;0) 298.15 GBCCAL#;" in text
    assert "PAR G(A2_BCC,Al:VA;0) 298.15 0.1+GBCCAL#;" in text


def test_tdb_writer_exports_factsage_style_database() -> None:
    """The TDB writer exports a FactSage-style database."""
    database = DatabaseIR(
        elements=[
            Element(symbol="AL", atomic_mass=26.98154, reference_state="FCC_A1"),
            Element(symbol="FE", atomic_mass=55.847, reference_state="BCC_A2"),
            Element(symbol="SI", atomic_mass=28.086, reference_state="DIAMOND_A4"),
            Element(symbol="VA", reference_state="VACUUM"),
        ],
        functions=[
            FunctionDefinition(name="GHSERAL", expression="0"),
            FunctionDefinition(name="GHSERFE", expression="0"),
            FunctionDefinition(name="GHSERSI", expression="0"),
        ],
        phases=[
            Phase(
                name="BCC_A2",
                constituents=[
                    ConstituentSet(1, ["AL", "FE"], 1),
                    ConstituentSet(2, ["VA"], 3),
                ],
            ),
            Phase(
                name="A2_BCC",
                constituents=[
                    ConstituentSet(1, ["AL", "FE"], 1),
                    ConstituentSet(2, ["VA"], 3),
                ],
            ),
            Phase(
                name="B2_BCC",
                constituents=[
                    ConstituentSet(1, ["AL", "FE"], 0.5),
                    ConstituentSet(2, ["AL", "FE"], 0.5),
                    ConstituentSet(3, ["VA"], 3),
                ],
            ),
            Phase(
                name="LT-FE1SI2",
                constituents=[
                    ConstituentSet(1, ["FE"], 0.33),
                    ConstituentSet(2, ["SI"], 0.67),
                ],
            ),
            Phase(
                name="T11_ALFESI",
                constituents=[
                    ConstituentSet(1, ["AL"], 85),
                    ConstituentSet(2, ["FE"], 30),
                    ConstituentSet(3, ["SI"], 15),
                ],
            ),
        ],
        parameters=[
            Parameter("BCC_A2", "G", ["AL:VA"], expression="0"),
            Parameter("BCC_A2", "G", ["FE:VA"], expression="0"),
            Parameter("A2_BCC", "G", ["AL:VA"], expression="0.1"),
            Parameter("A2_BCC", "G", ["FE:VA"], expression="0.1"),
            Parameter("B2_BCC", "G", ["AL:FE:VA"], expression="-1"),
            Parameter(
                "LT-FE1SI2",
                "G",
                ["FE:SI"],
                expression="-10+GHSERFE+2*GHSERSI",
            ),
            Parameter(
                "T11_ALFESI",
                "G",
                ["AL:FE:SI"],
                expression="-100+85*GHSERAL+30*GHSERFE+15*GHSERSI",
            ),
        ],
        tdb_commands=[
            TdbCommand(
                kind="type_definition",
                command="TYPE_DEF",
                args=[
                    "O",
                    "GES",
                    "AMEND_PHASE_DESCRIPTION",
                    "B2_BCC",
                    "DIS_PART",
                    "A2_BCC",
                ],
                parsed={
                    "action": "DIS_PART",
                    "ordered_phase": "B2_BCC",
                    "disordered_phase": "A2_BCC",
                },
            )
        ],
    )

    text = dumps_tdb(database, export_style="factsage")

    assert (
        "TYPE_DEFINITION A GES AMEND_PHASE_DESCRIPTION "
        "B2_BCC DISORDER_PART BCC_A2 !"
    ) in text
    assert "PHASE A2_BCC" not in text
    assert "PAR G(A2_BCC," not in text
    assert "LT-FE1SI2" not in text
    assert "PHASE B2_BCC % 3 1 1 6 !" in text
    assert "PHASE LTFE1SI2 % 2 1 2 !" in text
    assert "PHASE T11_ALFESI % 3 17 6 3 !" in text
    assert "PHASE B2_BCC_AL_FE_VA_ENDMBR % 2 1 1 !" in text
    assert "PAR G(B2_BCC,Al:Fe:VA;0) 298.15 -2;" in text
    assert "PAR G(LTFE1SI2,Fe:Si;0) 298.15 -30+3*GHSERFE#" in text
    assert "PAR G(T11_ALFESI,Al:Fe:Si;0) 298.15 -20+17*GHSERAL#" in text


@pytest.mark.parametrize(
    ("file_name", "expected_factsage_8_plus"),
    [
        ("AlCuMgSi_ORNL_FS73.dat", False),
        ("database_ChemSage-legacy.dat", False),
        ("AlCuMgSi_ORNL_FS83.dat", True),
        ("modern_database.dat", True),
    ],
)
def test_read_dat_infers_factsage_version_from_filename(
    file_name: str,
    expected_factsage_8_plus: bool,
) -> None:
    """read_dat infers the FactSage version from the file name."""
    assert _infer_factsage_8_plus(file_name) is expected_factsage_8_plus


def test_phase_selection_adds_canonical_disordered_partner(
    liu_alfesi_tdb_database: dict,
) -> None:
    """phase_selection adds the canonical disordered partner phase."""
    eq.list_phases(liu_alfesi_tdb_database, ["Al", "Fe", "Si"])

    selected = eq.phase_selection(["FCC_4SL"])

    assert "FCC_4SL" in selected
    assert "FCC_A1" in selected
    assert "A1_FCC" not in selected
