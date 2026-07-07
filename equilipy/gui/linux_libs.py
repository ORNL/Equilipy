"""Make PySide6's xcb platform plugin self-contained on Linux.

PySide6 wheels do not vendor a handful of X11 client libraries that Qt's
xcb platform plugin needs (Qt >= 6.5 aborts with "xcb-cursor0 is needed
to load the Qt xcb platform plugin" when they are absent).  Instead of
requiring every machine to export LD_LIBRARY_PATH, ``equilipy.gui
--setup-linux-libs`` copies the missing libraries into PySide6's own
``Qt/lib`` directory, which the plugin already searches through its
built-in RPATH.  The fix is per-environment and persists until PySide6
is reinstalled.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

# SONAMEs Qt's xcb platform plugin may need but PySide6 wheels do not
# vendor. Present on desktop distros; typically absent on minimal or
# server images.
_XCB_SONAMES = (
    "libxcb-cursor.so.0",
    "libxcb-icccm.so.4",
    "libxcb-image.so.0",
    "libxcb-keysyms.so.1",
    "libxcb-render-util.so.0",
    "libxcb-shape.so.0",
    "libxcb-xkb.so.1",
    "libxkbcommon-x11.so.0",
)

_INSTALL_HINTS = """\
Install the missing libraries with one of:
  Ubuntu/Debian:  sudo apt install libxcb-cursor0 libxcb-icccm4 libxcb-image0 \
libxcb-keysyms1 libxcb-render-util0 libxcb-shape0 libxcb-xkb1 libxkbcommon-x11-0
  RHEL/Fedora:    sudo dnf install xcb-util-cursor xcb-util-wm xcb-util-image \
xcb-util-keysyms xcb-util-renderutil libxkbcommon-x11
  no root:        conda install -c conda-forge xcb-util-cursor xcb-util-wm \
xcb-util-image xcb-util-keysyms xcb-util-renderutil libxkbcommon
then re-run: equilipy.gui --setup-linux-libs"""


def _ldconfig_sonames() -> set[str]:
    """Return SONAMEs the system dynamic loader already resolves."""
    for command in (["/sbin/ldconfig", "-p"], ["ldconfig", "-p"]):
        try:
            output = subprocess.run(
                command, capture_output=True, text=True, check=False
            ).stdout
        except OSError:
            continue
        if output:
            return {
                line.strip().split()[0]
                for line in output.splitlines()
                if "=>" in line
            }
    return set()


def _candidate_dirs() -> list[Path]:
    """Directories that may hold the libraries (conda, LD_LIBRARY_PATH)."""
    candidates: list[Path] = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        candidates.append(Path(conda_prefix) / "lib")
    conda = shutil.which("conda")
    if conda:
        try:
            base = subprocess.run(
                [conda, "info", "--base"],
                capture_output=True,
                text=True,
                check=False,
            ).stdout.strip()
        except OSError:
            base = ""
        if base:
            candidates.append(Path(base) / "lib")
    for entry in os.environ.get("LD_LIBRARY_PATH", "").split(":"):
        if entry:
            candidates.append(Path(entry))
    seen: set[Path] = set()
    unique: list[Path] = []
    for candidate in candidates:
        if candidate not in seen and candidate.is_dir():
            seen.add(candidate)
            unique.append(candidate)
    return unique


def ensure_xcb_libs(
    qt_lib_dir: Path,
    system_sonames: set[str],
    candidate_dirs: list[Path],
) -> tuple[list[str], list[str]]:
    """Copy missing xcb SONAMEs into qt_lib_dir; return (actions, missing)."""
    actions: list[str] = []
    missing: list[str] = []
    for soname in _XCB_SONAMES:
        if (qt_lib_dir / soname).exists():
            actions.append(f"{soname}: already in {qt_lib_dir}")
            continue
        if soname in system_sonames:
            actions.append(f"{soname}: provided by the system")
            continue
        source = next(
            (
                directory / soname
                for directory in candidate_dirs
                if (directory / soname).exists()
            ),
            None,
        )
        if source is None:
            missing.append(soname)
            continue
        # copy() follows symlinks, so the SONAME file lands as a regular
        # file regardless of how the source packaging links it.
        shutil.copy(source, qt_lib_dir / soname)
        actions.append(f"{soname}: copied from {source.parent}")
    return actions, missing


def setup_linux_libs() -> int:
    """Entry point for ``equilipy.gui --setup-linux-libs``."""
    if not sys.platform.startswith("linux"):
        print("--setup-linux-libs is only needed on Linux; nothing to do.")
        return 0
    try:
        import PySide6
    except ModuleNotFoundError:
        print(
            "PySide6 is not installed. Install the GUI extra first: "
            "pip install 'equilipy[gui]'"
        )
        return 1
    qt_lib_dir = Path(PySide6.__file__).resolve().parent / "Qt" / "lib"
    if not qt_lib_dir.is_dir():
        print(f"PySide6 Qt library directory not found: {qt_lib_dir}")
        return 1
    if not os.access(qt_lib_dir, os.W_OK):
        print(
            f"No write permission for {qt_lib_dir}.\n"
            "Re-run with sufficient permissions for this Python environment."
        )
        return 1

    try:
        actions, missing = ensure_xcb_libs(
            qt_lib_dir, _ldconfig_sonames(), _candidate_dirs()
        )
    except OSError as exc:
        print(f"Failed while copying libraries: {exc}")
        return 1
    for action in actions:
        print(action)
    if missing:
        print()
        print("Could not locate: " + ", ".join(missing))
        print(_INSTALL_HINTS)
        return 1
    print()
    print("Done. The GUI is ready: run `equilipy.gui`.")
    return 0
