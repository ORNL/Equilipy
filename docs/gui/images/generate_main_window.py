"""Regenerate the GUI screenshots for the docs.

Renders the real GUI offscreen with the production dark palette, so the
screenshots can be refreshed for each release without a display:

    python docs/gui/images/generate_main_window.py

Produces:
- main_window.png                      plain main window (gui/index.md)
- session_overview_session_controls.png
      Session Overview with example/GUI_example as the active directory and
      AlCuMgSi_ORNL_FS83.dat attached; the sidebar Open/New Session/Save
      buttons are highlighted (gui/calculation.md).
- equilibrium_module.png               Equilibrium#1 module page with a
      single Al-balance/Mg/Si wt% condition (gui/calculation.md).
- equilibrium_type_dialog.png          the calculation type dialog opened
      by the Type button, Batch condition selected (gui/calculation.md).
- solidification_module.png            Solidification#1 module page loaded
      from example/GUI_example/Solidification_1.eq (gui/calculation.md).
- assistant_panel.png                  wide capture with the Assistant panel
      open, showing the confirmed Codex conversation that created and ran
      Solidification#2 (gui/calculation.md, AI-assist section).
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# Strip inherited Claude/Anthropic session variables (present when this
# script is launched from inside another agent session) so the Claude CLI
# provider authenticates with the user's own login.
for _key in list(os.environ):
    if _key.startswith(("CLAUDE", "ANTHROPIC")):
        os.environ.pop(_key)

from PySide6.QtCore import QPoint, QRect, Qt
from PySide6.QtGui import QColor, QPainter, QPalette, QPen
from PySide6.QtWidgets import QApplication, QPushButton

from equilipy.database_ir import DatabaseIR
from equilipy.gui.app import _apply_dark_appearance
from equilipy.gui.calculation.state import CalculationDatabase
from equilipy.gui.main_window import DatabaseEditorWindow

HIGHLIGHT_BLUE = "#1687ff"  # matches the app's highlight color
SIDEBAR_MAX_X = 330  # sidebar width; excludes same-named buttons elsewhere


def _sidebar_session_buttons(window) -> list[QPushButton]:
    buttons = []
    for button in window.findChildren(QPushButton):
        if button.text() not in {"Open", "New Session", "Save"}:
            continue
        if not button.isVisible():
            continue
        if button.mapTo(window, QPoint(0, 0)).x() >= SIDEBAR_MAX_X:
            continue
        buttons.append(button)
    return buttons


def _union_rect_in_window(window, widgets) -> QRect:
    rect = QRect()
    for widget in widgets:
        top_left = widget.mapTo(window, QPoint(0, 0))
        rect = rect.united(QRect(top_left, widget.size()))
    return rect


def _highlight(pixmap, rect: QRect, padding: int = 6):
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    pen = QPen(QColor(HIGHLIGHT_BLUE))
    pen.setWidth(3)
    painter.setPen(pen)
    painter.drawRoundedRect(
        rect.adjusted(-padding, -padding, padding, padding), 10, 10
    )
    painter.end()
    return pixmap


def _attach_session_database(window, session, database_path: Path) -> None:
    """Attach a database to the session the way _load_session_database does,
    minus the file dialog."""
    database = window._load_runtime_database_for_calculation(database_path)
    session.database_counter += 1
    database_id = f"{session.id}:database:{session.database_counter}"
    for loaded in session.databases.values():
        loaded.selected = False
    session.databases[database_id] = CalculationDatabase(
        id=database_id,
        name=database_path.name,
        path=str(database_path),
        database=database,
        selected=True,
    )
    window._show_calculation_overview(session)
    window._update_current_calculation_heading(session)


def main() -> None:
    repo = Path(__file__).resolve().parents[3]
    gui_example = repo / "example" / "GUI_example"
    database_path = gui_example / "databases" / "AlCuMgSi_ORNL_FS83.dat"

    app = QApplication([])
    _apply_dark_appearance(app, Qt, QColor, QPalette)
    window = DatabaseEditorWindow(DatabaseIR(name=""))
    window.resize(1200, 760)
    window.show()
    app.processEvents()
    here = Path(__file__).parent

    window.grab().save(str(here / "main_window.png"), "PNG")
    print(f"Saved {here / 'main_window.png'}")

    # Session Overview screenshot: GUI_example active directory, one session
    # with the example database attached.
    window.active_calculation_directory = gui_example
    window._update_active_directory_button()
    session = window._active_calculation_session(create_if_missing=True)
    _attach_session_database(window, session, database_path)
    app.processEvents()

    rect = _union_rect_in_window(window, _sidebar_session_buttons(window))
    annotated = _highlight(window.grab(), rect)
    out = here / "session_overview_session_controls.png"
    annotated.save(str(out), "PNG")
    print(f"Saved {out}")

    # Equilibrium module page: single condition Al (balance), Mg 0.4, Si 6
    # in C / atm / wt%, matching the docs walkthrough.
    from equilipy.gui.calculation.session_serialization import (
        _apply_module_state_to_runtime,
    )

    session.temperature_unit = "C"
    session.pressure_unit = "atm"
    session.amount_unit = "wt%"
    module = window._add_calculation_module(session, "equilibrium")
    window._show_calculation_module(session, module)
    runtime = window._ensure_calculation_module_runtime(session, module)
    _apply_module_state_to_runtime(
        module,
        runtime,
        {
            "composition": [
                {"species": "Al", "amount": "?100"},
                {"species": "Mg", "amount": "0.4"},
                {"species": "Si", "amount": "6"},
            ],
            "temperature": "1600",
            "pressure": "1",
        },
    )
    window._confirm_calculation_phases(
        session.id,
        module.id,
        runtime["composition_table"],
        runtime["phase_tree"],
        runtime["status"],
    )
    window._show_calculation_module(session, module)
    app.processEvents()
    out = here / "equilibrium_module.png"
    window.grab().save(str(out), "PNG")
    print(f"Saved {out}")

    # Calculation type dialog (opened by the Type button), Batch selected.
    from PySide6.QtWidgets import QDialog

    module.calculation_type = "batch"
    dialog_out = here / "equilibrium_type_dialog.png"

    def grab_dialog_exec(dialog):
        dialog.show()
        app.processEvents()
        dialog.grab().save(str(dialog_out), "PNG")
        print(f"Saved {dialog_out}")
        dialog.close()
        return 0

    original_exec = QDialog.exec
    QDialog.exec = grab_dialog_exec
    try:
        window._open_calculation_type_dialog(module)
    finally:
        QDialog.exec = original_exec

    # Solidification module page, state loaded from the saved example module.
    from equilipy.gui.calculation.session_files import _read_eq_payload

    state = _read_eq_payload(gui_example / "Solidification_1.eq")["modules"][0]
    solidification = window._add_calculation_module(session, "solidification")
    solidification.solidification_model = state.get(
        "solidification_model", "scheil"
    )
    solidification.transition_search = bool(state.get("transition_search", True))
    solidification.start_from_liquidus = bool(
        state.get("start_from_liquidus", True)
    )
    window._show_calculation_module(session, solidification)
    runtime = window._ensure_calculation_module_runtime(session, solidification)
    _apply_module_state_to_runtime(solidification, runtime, state)
    window._confirm_calculation_phases(
        session.id,
        solidification.id,
        runtime["composition_table"],
        runtime["phase_tree"],
        runtime["status"],
    )
    window._show_calculation_module(session, solidification)
    app.processEvents()
    out = here / "solidification_module.png"
    window.grab().save(str(out), "PNG")
    print(f"Saved {out}")

    # Assistant panel: a REAL provider run. The prompt is sent to the Claude
    # CLI through the panel; the assistant replies with an action block and
    # Equilipy executes it, so the capture shows a genuine transcript.
    # Requires a logged-in Claude CLI; the screenshot is skipped (kept
    # as-is) when the provider is unavailable.
    import time

    # The docs' equilibrium example prompt verbatim: no "run it", so the
    # assistant prepares the module for review instead of launching the
    # solver.
    ASSISTANT_PROMPT = (
        "Using AlCuMgSi_ORNL_FS83.dat, add an equilibrium module in wt% "
        "with Al as balance, Mg 0.4, Si 6 at 700 C and 1 atm."
    )

    provider_label = os.environ.get("EQUILIPY_DOCS_ASSISTANT_PROVIDER", "Claude")
    panel = window.calculation_chat_sidebar
    provider_index = next(
        (
            index
            for index in range(panel.provider_select.count())
            if panel.provider_select.itemText(index) == provider_label
        ),
        None,
    )
    if provider_index is None:
        print(
            f"{provider_label} CLI not available; "
            "assistant_panel.png not regenerated."
        )
        return

    panel.provider_select.setCurrentIndex(provider_index)
    panel.setVisible(True)
    # Same 1200x760 aspect ratio as the other window captures, scaled up
    # to give the assistant panel room.
    window.resize(1440, 912)
    app.processEvents()

    def pump_until(condition, timeout_s: float, label: str) -> bool:
        deadline = time.monotonic() + timeout_s
        while time.monotonic() < deadline:
            app.processEvents()
            if condition():
                return True
            time.sleep(0.05)
        print(f"Timed out waiting for {label}.")
        return False

    panel.prompt.setText(ASSISTANT_PROMPT)
    panel.send_message()

    def _request_settled() -> bool:
        text = panel.transcript.toPlainText()
        if panel._thread is not None:
            return False
        # Success appends the executed-action summary; failures append an
        # error block. Either way the request is settled.
        return "Equilipy:" in text or "exited with status" in text or "[" in text

    # Provider reply arrives on a worker thread; actions execute on the main
    # thread once the reply completes and append the "Equilipy: ..." summary.
    pump_until(_request_settled, 480, "assistant reply and action execution")
    if "Equilipy:" not in panel.transcript.toPlainText():
        print("Assistant request did not complete cleanly; see capture.")
    # The requested calculation runs through the asynchronous solver lane.
    pump_until(
        lambda: not window._calculation_in_progress
        and not getattr(window, "_calculation_thread_refs", []),
        300,
        "assistant-launched calculation",
    )
    app.processEvents()
    out = here / "assistant_panel.png"
    window.grab().save(str(out), "PNG")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
