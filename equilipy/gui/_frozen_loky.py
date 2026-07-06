"""Frozen-app support for joblib/loky worker processes."""

from __future__ import annotations

import re
import runpy
import sys

_LOKY_POSIX_WORKER_MODULE = "joblib.externals.loky.backend.popen_loky_posix"
_RESOURCE_TRACKER_PATTERN = re.compile(
    r"\s*from\s+joblib\.externals\.loky\.backend\.resource_tracker"
    r"\s+import\s+main\s*;\s*"
    r"main\s*\(\s*(?P<fd>\d+)\s*,\s*(?P<verbose>\d+|True|False)\s*\)\s*"
)


def loky_freeze_support() -> None:
    """Route loky subprocess command lines before the GUI parser sees them.

    PyInstaller's ``multiprocessing.freeze_support()`` handles standard
    multiprocessing workers, but joblib's loky backend launches its own worker
    and resource-tracker command lines. In a frozen GUI app those arguments are
    sent back to the app executable and must be dispatched before ``argparse``.
    """
    if not getattr(sys, "frozen", False):
        return

    args = sys.argv[1:]
    if _is_loky_posix_worker_args(args):
        _run_loky_posix_worker(args)
    if _is_loky_win32_worker_args(args):
        _run_loky_win32_worker(args)

    tracker_args = _parse_loky_resource_tracker_args(args)
    if tracker_args is not None:
        _run_loky_resource_tracker(*tracker_args)


def _is_loky_posix_worker_args(args: list[str]) -> bool:
    """Return whether args are for a POSIX loky worker process."""
    return "--pipe" in args and "--process-name" in args


def _is_loky_win32_worker_args(args: list[str]) -> bool:
    """Return whether args are for a frozen Windows loky worker process."""
    return (
        len(args) >= 2
        and args[0] == "--multiprocessing-fork"
        and args[1].isdigit()
    )


def _filter_loky_posix_worker_args(args: list[str]) -> list[str]:
    """Strip the optional Python ``-m module`` prefix from loky worker args."""
    filtered_args: list[str] = []
    index = 0
    while index < len(args):
        arg = args[index]
        if arg == "-m":
            index += 1
            if index < len(args) and args[index] == _LOKY_POSIX_WORKER_MODULE:
                index += 1
            continue
        filtered_args.append(arg)
        index += 1
    return filtered_args


def _run_loky_posix_worker(args: list[str]) -> None:
    filtered_args = _filter_loky_posix_worker_args(args)
    sys.argv = [_LOKY_POSIX_WORKER_MODULE, *filtered_args]
    runpy.run_module(_LOKY_POSIX_WORKER_MODULE, run_name="__main__", alter_sys=True)
    raise SystemExit(0)


def _run_loky_win32_worker(args: list[str]) -> None:
    from joblib.externals.loky.backend.popen_loky_win32 import main

    sys.argv = [sys.argv[0], *args]
    main(pipe_handle=int(args[1]))
    raise SystemExit(0)


def _parse_loky_resource_tracker_args(args: list[str]) -> tuple[int, bool] | None:
    """Parse a frozen loky resource-tracker command, if present."""
    if "-c" not in args:
        return None

    command_index = args.index("-c") + 1
    if command_index >= len(args):
        return None

    match = _RESOURCE_TRACKER_PATTERN.fullmatch(args[command_index])
    if match is None:
        return None

    return int(match.group("fd")), _parse_verbose(match.group("verbose"))


def _parse_verbose(value: str) -> bool:
    """Parse the resource tracker verbose flag from loky's command string."""
    if value == "True":
        return True
    if value == "False":
        return False
    return bool(int(value))


def _run_loky_resource_tracker(fd: int, verbose: bool) -> None:
    from joblib.externals.loky.backend.resource_tracker import main

    main(fd, verbose)
    raise SystemExit(0)
