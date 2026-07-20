"""Phase-selection helpers for system preparation."""

from __future__ import annotations

from typing import Any

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

from .exceptions import InputConditionError


def ordered_solution_phase_names(
    database: dict[str, Any],
    phases: list[str] | None = None,
) -> list[str]:
    """Return ordered SUBOM phases using database order/disorder metadata.

    A phase is treated as ordered when it is a SUBOM solution phase and the
    parsed database maps it to a disordered companion through
    ``iDisorderedPhaseCS``.  This is the same source metadata used later when
    phase selection expands ordered phases to include their companions.
    """
    names = [str(name).strip() for name in database.get("cSolnPhaseNameCS", [])]
    models = [
        str(model).strip().upper()
        for model in database.get("cSolnPhaseTypeCS", [])
    ]
    disordered = list(database.get("iDisorderedPhaseCS", []))
    ordered_bases: set[str] = set()
    for index, name in enumerate(names):
        if not name:
            continue
        model = models[index] if index < len(models) else ""
        companion = int(disordered[index]) if index < len(disordered) else 0
        if model == "SUBOM" and companion > 0:
            ordered_bases.add(name.upper())

    candidates = phases if phases is not None else names
    ordered: list[str] = []
    for phase in candidates:
        phase_name = str(phase).strip()
        base_name = phase_name.split("#", 1)[0].upper()
        if base_name in ordered_bases:
            ordered.append(phase_name)
    return ordered


def default_scheil_phase_selection(
    database: dict[str, Any],
    phases: list[str],
) -> tuple[list[str], list[str]]:
    """Return Scheil's default selected phases and excluded ordered phases."""
    ordered = set(ordered_solution_phase_names(database, phases))
    selected = [phase for phase in phases if phase not in ordered]
    excluded = [phase for phase in phases if phase in ordered]
    return selected, excluded


def scheil_ordered_phase_exclusion_notice(excluded_phases: list[str]) -> str:
    """Return the user-visible notice for Scheil's ordered-phase default."""
    phase_text = ", ".join(excluded_phases)
    return (
        "Notice: Scheil cooling default phase selection excluded ordered "
        f"phase(s) {phase_text} for performance. Pass include_ordered=True or "
        "an explicit phases list to include them."
    )


def phase_selection(phases):
    """Revise variables passed to equilifort for phase selection.

    Description
    ===========
    This function revise variables to pass to equilifort module. Note that
    this function is active when system elements are defined.

    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     02/06/2024      S.Y. KWON       Original code


    Variables
    =========

    Input
    DBobject : Global variable object that includes parsed database
    phases : A list of phases selected by the user

    Output
    Revised DBobject
    """
    iPhasePS = np.ones(len(var.iPhaseCS), dtype=int) * (-1)
    iSolnPS = np.ones(var.nSolnPhasesSysCS, dtype=int) * (-1)

    # Check if phases are list
    if not isinstance(phases, list):
        raise InputConditionError("Make sure if the data type of phases is list.")
    unknown_phases = [
        phase for phase in phases if not _phase_name_is_available(str(phase))
    ]
    if unknown_phases:
        unknown_text = ", ".join(str(phase) for phase in unknown_phases)
        raise InputConditionError(
            f"Unknown phase(s) for selected system: {unknown_text}"
        )

    phases = _expand_composition_set_aliases(phases)
    phases = _canonical_phase_selection_order(
        _expand_required_disordered_phases(phases)
    )

    # Revise system-to-database maps and phase-name variables.
    iSys2DBSoln = list([])
    iSys2DBComp = list([])
    iSys2DBSpecies = list([])
    cPhaseNameSys = list([])
    PhaseTypeIDSys = list([])
    NewPhaseList = list([])
    iElementCS = np.arange(len(var.cElementNameCS))
    allowed_element_indices = set(int(index) for index in var.iElementDBIndex)
    allowed_element_indices.update(
        int(index)
        for index in getattr(var, "iPseudoComponentDependentElementDBIndex", [])
    )
    iElementDBIndex_reverse = [
        x for x in iElementCS if x not in allowed_element_indices
    ]

    # 1. Add each phase information into PS variables
    for phase in phases:
        # Check the type of phase name
        if not isinstance(phase, str):
            raise InputConditionError("Make sure if all elements in phases are string.")

        # For immiscible phases, the second and third phase must involve the
        # first phase.
        if "#" in phase:
            fphase = phase.split("#")[0]
            if fphase not in phases:
                raise InputConditionError(
                    f"Immiscibile phase {phase} was selected without "
                    f"involving {fphase}."
                )

        # Get the type of the phase and corresponding id in database
        p_type, p_id = var.PhaseNameSys[phase]
        NewPhaseList.append(phase)
        PhaseTypeIDSys.append([p_type, p_id])
        cPhaseNameSys.append(phase)

        if p_type == "soln":
            if p_id == 1:
                iFirst = 0
            else:
                iFirst = var.nSpeciesPhaseCS[int(p_id - 2)]
            iLast = var.nSpeciesPhaseCS[int(p_id - 1)]
            iPhasePS[iFirst:iLast] = p_id
            iSolnPS[int(p_id - 1)] = p_id
            iSys2DBSoln.append(p_id)
            for i in range(iFirst, iLast):
                if sum(var.dStoichSpeciesCS[i, iElementDBIndex_reverse]) < 1e-15:
                    iSys2DBSpecies.append(i)
        elif p_type == "compd":
            iPhasePS[p_id] = 0
            iSys2DBComp.append(p_id)
            iSys2DBSpecies.append(p_id)
        else:
            raise NameError(f"Error: PhaseSelection cannot identify type of {phase}")
    # print(iSys2DBSpecies)
    public_phase_list = _public_phase_selection_names(NewPhaseList)

    # Allocate python variables:
    var.iSys2DBSoln = iSys2DBSoln
    var.iSys2DBComp = iSys2DBComp
    var.iSys2DBSpecies = iSys2DBSpecies

    var.cPhaseNameSys = cPhaseNameSys
    var.PhaseNameSys = dict(zip(cPhaseNameSys, PhaseTypeIDSys, strict=False))
    var.iPhaseCS = iPhasePS

    # Allocate to fortran variables
    fort.moduleparsecs.iphasecs = iPhasePS
    fort.modulethermo.isolnps = iSolnPS
    return public_phase_list


def _expand_required_disordered_phases(phases: list[str]) -> list[str]:
    """Include DIS_PART helpers and physical disordered phases for SUBOM."""
    if not hasattr(var, "iDisorderedPhaseCS"):
        return phases

    expanded = list(phases)
    selected = set(expanded)
    name_by_soln_id: dict[int, str] = {}
    for phase_name, phase_info in var.PhaseNameSys.items():
        if len(phase_info) != 2:
            continue
        phase_type, phase_id = phase_info
        if phase_type == "soln":
            name_by_soln_id[int(phase_id)] = phase_name

    for phase in list(expanded):
        phase_type, phase_id = var.PhaseNameSys[phase]
        if phase_type != "soln":
            continue
        index = int(phase_id) - 1
        if index < 0 or index >= len(var.iDisorderedPhaseCS):
            continue
        disordered_phase_id = int(var.iDisorderedPhaseCS[index])
        if disordered_phase_id <= 0:
            continue
        disordered_phase = name_by_soln_id.get(disordered_phase_id)
        if disordered_phase is None:
            raise InputConditionError(
                f"Ordered phase {phase} requires disordered phase "
                f"#{disordered_phase_id}, but it is not available for the "
                "selected system."
            )
        if disordered_phase not in selected:
            expanded.append(disordered_phase)
            selected.add(disordered_phase)
        standalone_ids = getattr(var, "iOrderDisorderStandalonePhaseCS", [])
        standalone_id = (
            int(standalone_ids[index]) if index < len(standalone_ids) else 0
        )
        standalone_phase = name_by_soln_id.get(standalone_id)
        if standalone_phase is not None and standalone_phase not in selected:
            expanded.append(standalone_phase)
            selected.add(standalone_phase)

    return expanded


def _phase_name_is_available(phase: str) -> bool:
    """Return True when a phase or one of its composition-set aliases exists."""
    if phase in var.PhaseNameSys:
        return True
    first_composition_set = _first_composition_set_alias_base(phase)
    if first_composition_set is not None:
        return True
    prefix = f"{phase}#"
    return any(str(name).startswith(prefix) for name in var.PhaseNameSys)


def _expand_composition_set_aliases(phases: list[str]) -> list[str]:
    """Select every composition set when the base phase name is requested."""
    expanded: list[str] = []
    selected: set[str] = set()
    requested = set(phases)

    for phase in phases:
        first_composition_set = _first_composition_set_alias_base(phase)
        if first_composition_set is not None:
            if first_composition_set not in selected:
                expanded.append(first_composition_set)
                selected.add(first_composition_set)
            continue

        if phase not in var.PhaseNameSys:
            continue
        if "#" in phase:
            if phase not in selected:
                expanded.append(phase)
                selected.add(phase)
            continue

        prefix = f"{phase}#"
        aliases = [
            candidate
            for candidate in var.PhaseNameSys
            if candidate == phase or str(candidate).startswith(prefix)
        ]
        if len(aliases) == 1:
            aliases = [phase]
        for alias in aliases:
            if alias in selected:
                continue
            expanded.append(alias)
            selected.add(alias)

    for phase in requested:
        if phase in selected or phase not in var.PhaseNameSys:
            continue
        expanded.append(phase)
        selected.add(phase)

    return expanded


def _first_composition_set_alias_base(phase: str) -> str | None:
    """Return the base phase name for an explicit ``#1`` composition-set alias."""
    if "#" not in phase:
        return None
    base_phase, composition_set = phase.rsplit("#", 1)
    if composition_set != "1":
        return None
    if base_phase not in var.PhaseNameSys:
        return None
    return base_phase


def _canonical_phase_selection_order(phases: list[str]) -> list[str]:
    """Return selected phases in the same order used by the Fortran backend."""
    selected = set(phases)
    return [phase for phase in var.PhaseNameSys if phase in selected]


def _public_phase_selection_names(phases: list[str]) -> list[str]:
    """Return selected phase names without leaking hidden DIS_PART helpers."""
    hidden_helpers = {
        str(name).strip()
        for name in getattr(var, "cOrderDisorderHelperPhaseNames", [])
    }
    name_by_soln_id = {
        int(phase_info[1]): str(phase_name)
        for phase_name, phase_info in var.PhaseNameSys.items()
        if len(phase_info) == 2 and phase_info[0] == "soln"
    }
    helper_public_names: dict[str, str | None] = {}
    disordered = getattr(var, "iDisorderedPhaseCS", [])
    standalone = getattr(var, "iOrderDisorderStandalonePhaseCS", [])
    for ordered_index, helper_id_raw in enumerate(disordered):
        helper_id = int(helper_id_raw)
        if helper_id <= 0:
            continue
        helper_name = name_by_soln_id.get(helper_id)
        if helper_name is None or helper_name not in hidden_helpers:
            continue
        standalone_id = (
            int(standalone[ordered_index]) if ordered_index < len(standalone) else 0
        )
        helper_public_names[helper_name] = name_by_soln_id.get(standalone_id)

    public_names: list[str] = []
    selected: set[str] = set()
    for phase in phases:
        public_phase = phase
        if phase in hidden_helpers:
            public_phase = helper_public_names.get(phase)
            if public_phase is None:
                continue
        if public_phase in selected:
            continue
        public_names.append(public_phase)
        selected.add(public_phase)
    return public_names
