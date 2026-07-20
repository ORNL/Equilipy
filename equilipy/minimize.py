"""Wrappers around Fortran minimization routines."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

from .exceptions import EquilibError, SolverFailureReason

_DEBUG_MODE = False


def _thermo_info() -> int:
    """Return the scalar Fortran thermodynamic status code."""
    value = fort.modulethermoio.infothermo
    if hasattr(value, "item"):
        value = value.item()
    return int(value)


def _gem_exit_status() -> int:
    """Return the scalar GEM exit status exposed by the Fortran solver."""
    value = getattr(fort.modulegemsolver, "igemexitstatus", 0)
    if hasattr(value, "item"):
        value = value.item()
    return int(value)


@dataclass(frozen=True, slots=True)
class SolverStatus:
    """Raw terminal solver facts and their named public failure reason."""

    converged: bool
    thermo_status: int
    gem_exit_status: int
    phase_change_reason: int
    pea_exit_status: int
    pea_exit_reason: int
    grid_fallback_reason: int
    residual_norm: float
    failure_reason: SolverFailureReason | None

    @property
    def failed(self) -> bool:
        """Return whether the raw facts describe a failed terminal state."""
        return self.failure_reason is not None


def classify_solver_failure(
    *,
    stage: str,
    converged: bool,
    thermo_status: int,
    gem_exit_status: int,
    pea_exit_status: int,
    pea_exit_reason: int,
    grid_fallback_reason: int,
    residual_norm: float,
) -> SolverFailureReason | None:
    """Classify raw solver facts without changing any of them.

    Parameters
    ----------
    stage
        Name of the completed Fortran stage.
    converged
        Raw Lagrangian convergence flag.
    thermo_status
        Raw ``INFOThermo`` value.
    gem_exit_status
        Raw GEM exit status.
    pea_exit_status
        Raw strict PEA certificate status.
    pea_exit_reason
        Raw PEA exit-path reason.
    grid_fallback_reason
        Raw static-grid failure reason.
    residual_norm
        Final GEM residual norm.

    Returns
    -------
    SolverFailureReason or None
        Named failure reason, or ``None`` for a coherent successful state.
    """
    gem = fort.modulegemsolver
    failure_present = (
        thermo_status != 0
        or gem_exit_status != int(getattr(gem, "gem_exit_status_ok", 0))
        or not converged
    )

    if failure_present and grid_fallback_reason != int(
        getattr(gem, "grid_fallback_none", 0)
    ):
        if grid_fallback_reason == int(
            getattr(gem, "grid_fallback_lagrangian", 2)
        ):
            return SolverFailureReason.PEALG_FAILED
        if grid_fallback_reason == int(
            getattr(gem, "grid_fallback_certificate", 3)
        ):
            return SolverFailureReason.DF_SWEEP_UNKNOWN
        return SolverFailureReason.GRID_SEED_INVALID

    if failure_present and pea_exit_reason == int(
        getattr(gem, "pea_exit_reason_max_iter", 8)
    ):
        return SolverFailureReason.PEA_MAX

    if stage in {"pealg", "pea_lagrangian"} and failure_present:
        return SolverFailureReason.PEALG_FAILED
    if stage in {"df_sweep", "candidate_sweep"} and failure_present:
        return SolverFailureReason.DF_SWEEP_UNKNOWN
    if stage == "postprocess" and failure_present:
        return SolverFailureReason.POSTPROCESS_FAILED

    if failure_present and gem_exit_status == int(
        getattr(gem, "gem_exit_status_pea_uncertified_handoff", 2)
    ):
        error_reasons = {
            int(getattr(gem, "pea_exit_reason_error_direct_handoff", 1)),
            int(getattr(gem, "pea_exit_reason_error_leveling", 2)),
        }
        if pea_exit_reason in error_reasons:
            return SolverFailureReason.DF_SWEEP_UNKNOWN
        return SolverFailureReason.PEA_UNCERTIFIED

    if thermo_status != 0:
        return SolverFailureReason.FORTRAN_STAGE_FAILED
    if not converged and gem_exit_status == int(
        getattr(gem, "gem_exit_status_ok", 0)
    ):
        return SolverFailureReason.STATUS_INCONSISTENT
    if not converged or gem_exit_status != int(
        getattr(gem, "gem_exit_status_ok", 0)
    ):
        return SolverFailureReason.LAGRANGIAN_UNCONVERGED
    if residual_norm > 1e-5:
        return SolverFailureReason.RESIDUAL_NOT_CONVERGED
    return None


def capture_solver_status(stage: str = "gemsolver") -> SolverStatus:
    """Capture immutable raw solver facts after a Fortran stage."""
    gem = fort.modulegemsolver
    converged = bool(getattr(gem, "lconverged", False))
    thermo_status = _thermo_info()
    gem_exit_status = _gem_exit_status()
    phase_change_reason = int(getattr(gem, "iphasechangereason", 0))
    pea_exit_status = int(getattr(gem, "ipeaexitstatus", 0))
    pea_exit_reason = int(getattr(gem, "ipeaexitreason", 0))
    grid_fallback_reason = int(getattr(gem, "igridfallbackreason", 0))
    residual_norm = float(getattr(gem, "dgemfunctionnorm", np.inf))
    failure_reason = classify_solver_failure(
        stage=stage,
        converged=converged,
        thermo_status=thermo_status,
        gem_exit_status=gem_exit_status,
        pea_exit_status=pea_exit_status,
        pea_exit_reason=pea_exit_reason,
        grid_fallback_reason=grid_fallback_reason,
        residual_norm=residual_norm,
    )
    return SolverStatus(
        converged=converged,
        thermo_status=thermo_status,
        gem_exit_status=gem_exit_status,
        phase_change_reason=phase_change_reason,
        pea_exit_status=pea_exit_status,
        pea_exit_reason=pea_exit_reason,
        grid_fallback_reason=grid_fallback_reason,
        residual_norm=residual_norm,
        failure_reason=failure_reason,
    )


def set_debug_mode(enabled: bool = True) -> None:
    """Enable or disable Fortran minimizer debug diagnostics."""
    global _DEBUG_MODE
    _DEBUG_MODE = bool(enabled)
    _apply_debug_mode()


def debug_mode_enabled() -> bool:
    """Return whether Fortran minimizer debug diagnostics are requested."""
    return _DEBUG_MODE


def _apply_debug_mode() -> None:
    """Synchronize the Python debug switch to the Fortran solver module."""
    try:
        fort.modulegemsolver.ldebugmode = bool(_DEBUG_MODE)
    except AttributeError:
        return


def _raise_fortran_status(stage: str) -> None:
    """Raise if a Fortran stage set INFOThermo without throwing."""
    info = _thermo_info()
    if info != 0:
        raise EquilibError(
            f"Equilifort {stage} failed: infothermo={info}",
            reason=SolverFailureReason.FORTRAN_STAGE_FAILED,
        )


def _raise_solver_status(stage: str) -> None:
    """Raise with immutable raw facts when a terminal solver stage failed."""
    status = capture_solver_status(stage)
    if not status.failed:
        return
    raise EquilibError(
        "Equilibrium calculation failed: "
        f"reason={status.failure_reason.value}; "
        f"converged={status.converged}; "
        f"infothermo={status.thermo_status}; "
        f"gem_exit_status={status.gem_exit_status}; "
        f"phase_change_reason={status.phase_change_reason}; "
        f"pea_exit_status={status.pea_exit_status}; "
        f"pea_exit_reason={status.pea_exit_reason}; "
        f"grid_fallback_reason={status.grid_fallback_reason}; "
        f"dGEMFunctionNorm={status.residual_norm}",
        reason=status.failure_reason,
        solver_status=status,
    )


def _run_fortran_stage(stage: str) -> None:
    """Run a Fortran stage and surface either Python or INFOThermo failures."""
    _apply_debug_mode()
    try:
        getattr(fort, stage)()
    except Exception as e:
        info = _thermo_info()
        raise EquilibError(
            f"Equilifort {stage} failed: infothermo={info}, Error: {e}",
            reason=SolverFailureReason.FORTRAN_STAGE_FAILED,
        ) from e
    if stage in {"gemsolver", "postprocess", "runlagrangiangem"}:
        _raise_solver_status(stage)
    else:
        _raise_fortran_status(stage)


def prepare_minimization() -> None:
    """Prepare Fortran thermodynamic state for staged minimization calls."""
    _run_fortran_stage("checkthermoinput")
    _run_fortran_stage("initthermo")
    _run_fortran_stage("checksystem")
    _sync_python_system_metadata_after_checksystem()
    _run_fortran_stage("compthermodata")
    _run_fortran_stage("checkthermodata")


def run_leveling() -> None:
    """Run the Leveling stage after thermodynamic state has been prepared."""
    _run_fortran_stage("runleveling")


def check_phase_assemblage() -> None:
    """Run the PEA phase-assemblage check stage."""
    _run_fortran_stage("checkphaseassemblage")


def run_lagrangian_gem() -> None:
    """Run Lagrangian GEM for the current phase assemblage only."""
    _run_fortran_stage("runlagrangiangem")


def minimize():
    """Run the full equilibrium minimization workflow."""
    prepare_minimization()

    # Estimate the equilibrium phase assemblage and other important properties
    # using the Leveling algorithm:
    _run_fortran_stage("gemsolver")
    equilibrium_status = capture_solver_status("gemsolver")

    # Perform post-processing calculations of results:
    try:
        _run_fortran_stage("postprocess")
    except EquilibError as error:
        # Postprocess may fail after GEM has already certified the requested-T
        # equilibrium.  Preserve that independent fact so public result
        # construction can retain the valid equilibrium as a warning result.
        error.equilibrium_status = equilibrium_status
        raise

    return None


def _sync_python_system_metadata_after_checksystem() -> None:
    """Mirror Fortran's screened system definition back to Python globals."""
    try:
        species_pass = np.asarray(fort.modulethermo.ispeciespass, dtype=int)
        element_system = np.asarray(fort.modulethermo.ielementsystem, dtype=int)
    except (AttributeError, ValueError):
        return

    if species_pass.size == 0 or element_system.size == 0:
        return

    active_element_indices = [
        index
        for index, value in enumerate(element_system[: var.nElementsCS])
        if int(value) > 0 or int(value) == -1
    ]
    var.iElementDBIndex = np.asarray(active_element_indices, dtype=int)
    var.iElementSysIndex = np.arange(len(active_element_indices), dtype=int)
    if active_element_indices:
        var.iElementSys = np.asarray(
            [
                var.cPeriodicTable[str(var.cElementNameCS[index]).strip()][0]
                for index in active_element_indices
            ],
            dtype=int,
        )
        var.cComponentNameSys = [
            str(var.cElementNameCS[index]).strip()
            for index in active_element_indices
        ]

    species_indices: list[int] = []
    solution_indices: list[int] = []
    phase_names: list[str] = []
    solution_count = int(var.nSolnPhasesSysCS)
    species_phase_cs = np.asarray(var.nSpeciesPhaseCS, dtype=int)

    for phase_index in range(solution_count):
        first = 0 if phase_index == 0 else int(species_phase_cs[phase_index - 1])
        last = int(species_phase_cs[phase_index])
        passed = [
            species_index
            for species_index in range(first, last)
            if int(species_pass[species_index]) > 0
        ]
        if not passed:
            continue
        solution_indices.append(phase_index + 1)
        species_indices.extend(passed)
        phase_names.append(str(var.cSolnPhaseNameCS[phase_index]).strip())

    compound_indices: list[int] = []
    pure_start = int(species_phase_cs[solution_count - 1]) if solution_count > 0 else 0
    for species_index in range(pure_start, len(species_pass)):
        if int(species_pass[species_index]) <= 0:
            continue
        compound_indices.append(species_index)
        species_indices.append(species_index)
        phase_names.append(str(var.cSpeciesNameCS[species_index]).strip())

    var.iSys2DBSoln = np.asarray(solution_indices, dtype=int)
    var.iSys2DBComp = np.asarray(compound_indices, dtype=int)
    var.iSys2DBSpecies = np.asarray(species_indices, dtype=int)
    var.cSpeciesNameSys = np.asarray(var.cSpeciesNameCS)[var.iSys2DBSpecies]
    var.cPhaseNameSys = _with_composition_set_aliases(phase_names)

    var.iSys2DBSolnDefault = var.iSys2DBSoln
    var.cPhaseNameSysDefault = var.cPhaseNameSys
    var.iSys2DBCompDefault = var.iSys2DBComp
    var.iSys2DBSpeciesDefault = var.iSys2DBSpecies


def _with_composition_set_aliases(phase_names: list[str]) -> list[str]:
    """Append #n aliases to repeated phase names while preserving order."""
    counts: dict[str, int] = {}
    output: list[str] = []
    for phase_name in phase_names:
        counts[phase_name] = counts.get(phase_name, 0) + 1
        if counts[phase_name] > 1:
            output.append(f"{phase_name}#{counts[phase_name]}")
        else:
            output.append(phase_name)
    return output
