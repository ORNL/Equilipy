"""PySide6 application entry point for the Equilipy database editor."""

from __future__ import annotations

import argparse
import faulthandler
import os
import sys
from datetime import datetime
from multiprocessing import freeze_support
from pathlib import Path

from equilipy.database_ir import DatabaseIR, load_database_ir
from equilipy.gui._frozen_loky import loky_freeze_support

# Keep the crash-log handle alive for the process lifetime; faulthandler
# writes to the raw file descriptor when a native crash occurs.
_CRASH_LOG_FILE = None


def _enable_crash_log() -> None:
    """Route native-crash tracebacks to a persistent log file.

    The Fortran solver runs inside the GUI process, so a solver segfault
    kills the app with no dialog.  faulthandler cannot prevent that, but it
    writes the Python stack of every thread to this file at crash time, so
    the next crash is diagnosable instead of silent.
    """
    global _CRASH_LOG_FILE
    try:
        crash_path = Path.home() / ".equilipy" / "gui_crash.log"
        crash_path.parent.mkdir(parents=True, exist_ok=True)
        _CRASH_LOG_FILE = crash_path.open("a", encoding="utf-8")
        _CRASH_LOG_FILE.write(
            f"\n=== Equilipy GUI session {datetime.now().isoformat()} ===\n"
        )
        _CRASH_LOG_FILE.flush()
        faulthandler.enable(file=_CRASH_LOG_FILE)
    except OSError:
        faulthandler.enable()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Launch the Equilipy database editor.")
    parser.add_argument(
        "database_ir_json",
        nargs="?",
        help="Optional DatabaseIR JSON or TDB file to open.",
    )
    parser.add_argument(
        "--setup-linux-libs",
        action="store_true",
        help=(
            "Copy X11 client libraries missing from PySide6 into its Qt/lib "
            "directory so the GUI runs without system packages (Linux only)."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Launch the PySide6 database editor."""
    loky_freeze_support()
    freeze_support()
    _enable_crash_log()
    parser = _build_parser()
    args = _parse_launch_args(parser, argv)

    if getattr(args, "setup_linux_libs", False):
        from .linux_libs import setup_linux_libs

        return setup_linux_libs()

    try:
        from PySide6.QtCore import Qt, QTimer
        from PySide6.QtGui import QColor, QIcon, QPalette
        from PySide6.QtWidgets import QApplication
    except ModuleNotFoundError:
        print(
            "PySide6 is required for the GUI. Install it with "
            "`pip install -e '.[gui]'` or `pip install PySide6`.",
            file=sys.stderr,
        )
        return 1

    from .assets import equilipy_app_icon_path
    from .log_capture import GuiLogCapture
    from .main_window import DatabaseEditorWindow

    database = (
        load_database_ir(args.database_ir_json)
        if args.database_ir_json
        else DatabaseIR(name="")
    )

    app = QApplication(sys.argv[:1])
    log_capture = GuiLogCapture()
    log_capture.start()
    app.aboutToQuit.connect(log_capture.stop)
    _apply_dark_appearance(app, Qt, QColor, QPalette)
    app_icon_path = equilipy_app_icon_path()
    if app_icon_path.exists():
        app.setWindowIcon(QIcon(str(app_icon_path)))
    window = DatabaseEditorWindow(database)
    window.attach_log_capture(log_capture)
    if app_icon_path.exists():
        window.setWindowIcon(QIcon(str(app_icon_path)))
    window.resize(1200, 760)
    window.show()
    if os.environ.get("EQUILIPY_GUI_SMOKE") == "1":
        # Smoke mode (CI): run one event-loop pass through the full
        # startup path, then quit so the process exits with app.exec()'s
        # return code.
        QTimer.singleShot(0, app.quit)
    return app.exec()


def _parse_launch_args(
    parser: argparse.ArgumentParser,
    argv: list[str] | None,
) -> argparse.Namespace:
    """Parse GUI launch args, tolerating frozen-app launcher metadata."""
    if getattr(sys, "frozen", False):
        args, _unknown = parser.parse_known_args(_strip_frozen_launcher_args(argv))
        return args
    return parser.parse_args(argv)


def _strip_frozen_launcher_args(argv: list[str] | None) -> list[str] | None:
    """Remove unknown frozen-app launcher options before parsing positionals."""
    if argv is None:
        argv = sys.argv[1:]

    stripped: list[str] = []
    skip_next = False
    for index, arg in enumerate(argv):
        if skip_next:
            skip_next = False
            continue
        if arg.startswith("-"):
            if index + 1 < len(argv) and not argv[index + 1].startswith("-"):
                skip_next = True
            continue
        stripped.append(arg)
    return stripped


def _apply_dark_appearance(app, qt, qcolor, qpalette) -> None:
    """Request native dark chrome and install a matching widget palette."""
    style_hints = app.styleHints()
    if hasattr(style_hints, "setColorScheme"):
        style_hints.setColorScheme(qt.ColorScheme.Dark)

    palette = qpalette()
    palette.setColor(qpalette.ColorRole.Window, qcolor("#1f2526"))
    palette.setColor(qpalette.ColorRole.WindowText, qcolor("#eceff0"))
    palette.setColor(qpalette.ColorRole.Base, qcolor("#202627"))
    palette.setColor(qpalette.ColorRole.AlternateBase, qcolor("#252b2c"))
    palette.setColor(qpalette.ColorRole.ToolTipBase, qcolor("#30383a"))
    palette.setColor(qpalette.ColorRole.ToolTipText, qcolor("#eceff0"))
    palette.setColor(qpalette.ColorRole.Text, qcolor("#eceff0"))
    palette.setColor(qpalette.ColorRole.Button, qcolor("#2a3132"))
    palette.setColor(qpalette.ColorRole.ButtonText, qcolor("#eceff0"))
    palette.setColor(qpalette.ColorRole.BrightText, qcolor("#ffffff"))
    palette.setColor(qpalette.ColorRole.Highlight, qcolor("#1687ff"))
    palette.setColor(qpalette.ColorRole.HighlightedText, qcolor("#ffffff"))
    palette.setColor(qpalette.ColorRole.Link, qcolor("#1687ff"))
    app.setPalette(palette)


if __name__ == "__main__":
    raise SystemExit(main())
