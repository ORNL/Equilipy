"""PyInstaller entry point for the Equilipy GUI."""

from __future__ import annotations

import sys
from multiprocessing import freeze_support
from pathlib import Path

_SMOKE_ARG = "--equilipy-smoke-joblib"
_CALCULATION_SMOKE_ARG = "--equilipy-smoke-calculation"


def _run_joblib_smoke() -> int:
    """Run a tiny frozen-app joblib task to verify loky worker startup."""
    from equilipy._parallel import starmap_joblib

    updates: list[tuple[int, int, str]] = []
    result = starmap_joblib(
        pow,
        [(2, 1), (2, 2), (2, 3)],
        2,
        progress_callback=lambda current, total, label: updates.append(
            (current, total, label)
        ),
        desc="Frozen joblib smoke",
    )
    if result != [2, 4, 8]:
        raise RuntimeError(f"unexpected joblib smoke result: {result!r}")
    if not updates or updates[-1] != (3, 3, "Frozen joblib smoke"):
        raise RuntimeError(f"unexpected joblib smoke progress: {updates!r}")
    print("Equilipy frozen joblib smoke passed.", flush=True)
    return 0


def _smoke_database_path() -> Path:
    """Return the bundled or source-tree database used by frozen smoke tests."""
    candidates: list[Path] = []
    bundle_root = getattr(sys, "_MEIPASS", None)
    if bundle_root:
        candidates.append(Path(bundle_root) / "database" / "AlFeSi_99Liu.dat")
    candidates.append(
        Path(__file__).resolve().parents[2] / "database" / "AlFeSi_99Liu.dat"
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("AlFeSi_99Liu.dat was not bundled for smoke testing.")


def _run_calculation_smoke() -> int:
    """Run frozen-app parallel equilibrium and Scheil calculations."""
    import math

    import equilipy as eq

    database = eq.read_dat(str(_smoke_database_path()))
    phases = list(eq.list_phases(database, ["Al", "Fe"]))
    condition = {
        "T": [1600, 1500],
        "P": [1, 1],
        "Al": [0.75, 0.80],
        "Fe": [0.25, 0.20],
    }

    equilib_result = eq.equilib_batch(
        database,
        condition,
        phases=phases,
        n_cpu=2,
    )
    if len(equilib_result.points) != 2:
        raise RuntimeError(
            f"unexpected equilibrium smoke point count: {len(equilib_result.points)}"
        )
    if any(point.G is None or math.isnan(point.G) for point in equilib_result.points):
        raise RuntimeError("equilibrium smoke returned failed points.")

    scheil_result = eq.scheil_batch(
        "LIQUID",
        database,
        condition,
        delta_T=500,
        phases=phases,
        n_cpu=2,
    )
    if "Error" in scheil_result and any(
        str(error) != "0.0" for error in scheil_result["Error"]
    ):
        raise RuntimeError(f"Scheil smoke returned errors: {scheil_result['Error']!r}")
    if "T [K]" not in scheil_result or len(scheil_result["T [K]"]) < 2:
        raise RuntimeError("Scheil smoke did not produce result rows.")

    print("Equilipy frozen calculation smoke passed.", flush=True)
    return 0


if __name__ == "__main__":
    from equilipy.gui._frozen_loky import loky_freeze_support

    loky_freeze_support()
    freeze_support()
    if _SMOKE_ARG in sys.argv[1:]:
        raise SystemExit(_run_joblib_smoke())
    if _CALCULATION_SMOKE_ARG in sys.argv[1:]:
        raise SystemExit(_run_calculation_smoke())

    from equilipy.gui.app import main

    raise SystemExit(main())
