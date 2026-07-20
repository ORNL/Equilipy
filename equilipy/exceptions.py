"""Custom exception types used by Equilipy."""

from __future__ import annotations

from enum import Enum


class SolverFailureReason(str, Enum):
    """Named public reason for a failed equilibrium solver stage."""

    GRID_SEED_INVALID = "grid_seed_invalid"
    PEALG_FAILED = "pealg_failed"
    DF_SWEEP_UNKNOWN = "df_sweep_unknown"
    PEA_MAX = "pea_max"
    PEA_UNCERTIFIED = "pea_uncertified"
    LAGRANGIAN_UNCONVERGED = "lagrangian_unconverged"
    POSTPROCESS_FAILED = "postprocess_failed"
    RESIDUAL_NOT_CONVERGED = "residual_not_converged"
    FORTRAN_STAGE_FAILED = "fortran_stage_failed"
    STATUS_INCONSISTENT = "status_inconsistent"
    EQUILIBRIUM_FAILED = "equilibrium_failed"


class SolidificationWarningCode(str, Enum):
    """Named public warning emitted by a solidification calculation."""

    LIQUIDUS_SEARCH_FAILED = "LIQUIDUS_SEARCH_FAILED"


def liquidus_search_failure_warning(exc: Exception) -> str:
    """Return the report-ready warning for a failed liquidus search."""
    return (
        f"Warning [{SolidificationWarningCode.LIQUIDUS_SEARCH_FAILED.value}]: "
        "The liquidus search failed, so cooling is starting from the input "
        f"temperature. Cause: {exc}"
    )


class EquilibError(Exception):
    """Raise when equilibrium calculation fails."""

    def __init__(
        self,
        message="Equilibrium calculation failed.",
        *,
        reason: SolverFailureReason | None = None,
        solver_status=None,
        equilibrium_status=None,
    ):
        self.message = message
        self.reason = reason
        self.solver_status = solver_status
        self.equilibrium_status = equilibrium_status
        super().__init__(self.message)


class DatabaseParsingError(Exception):
    """Raise when database parsing fails."""

    def __init__(self, message="Failed to parse the database"):
        self.message = message
        super().__init__(self.message)


class DatabaseLoadError(Exception):
    """Raise when a parsed database cannot be loaded into global state."""

    def __init__(self, message="Failed to load the database"):
        self.message = message
        super().__init__(self.message)


class InputConditionError(Exception):
    """Raise when input conditions cannot be interpreted."""

    def __init__(self, message="Input condition is not recognized properly"):
        self.message = message
        super().__init__(self.message)


class PostProcessError(Exception):
    """Raise when post-processing fails."""

    def __init__(self, message="Failed processing calcualtion result."):
        self.message = message
        super().__init__(self.message)


class TransitionError(Exception):
    """Raise when transition search cannot produce a valid transition."""

    def __init__(self, message="Failed to find a phase transition."):
        self.message = message
        super().__init__(self.message)
