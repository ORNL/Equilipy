"""Run calculations and render calculation outputs."""

from __future__ import annotations

# ruff: noqa: F401,F403,F405,I001,E501

import copy
import json
import os
import pickle
import shutil
from pathlib import Path
from typing import Any

from PySide6.QtCore import QModelIndex, QSize, Qt, QThread, QTimer
from PySide6.QtGui import QAction, QIcon
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPlainTextEdit,
    QProgressDialog,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QStyle,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from equilipy.composition import expand_condition_species as _expand_condition_species
from equilipy.phase_selection import scheil_ordered_phase_exclusion_notice
from equilipy.results import ResultTable
from equilipy.results.equilib import EquilibResult
from equilipy.results.scheil import ScheilResult

from equilipy.gui.compat_exports import *
from equilipy.gui.widgets import forms as _forms
from equilipy.gui.widgets.composition import (
    COMPOSITION_INPUT_MIN_HEIGHT,
    COMPOSITION_ROW_HEIGHT,
    CompositionEditor,
)
from equilipy.gui.calculation.state import CalculationDatabase, CalculationModule, CalculationSession
from equilipy.gui.calculation.workers import CalculationCanceled, CalculationResultReceiver, CalculationWorker

_apply_button_cursor = _forms.apply_button_cursor
_collapsible_form_section = _forms.collapsible_form_section
_disable_inner_scrollbars = _forms.disable_inner_scrollbars
_fit_table_columns_to_headers = _forms.fit_table_columns_to_headers
_fit_table_to_rows = _forms.fit_table_to_rows
_fit_tree_to_items = _forms.fit_tree_to_items
_form_layout = _forms.form_layout
_form_section = _forms.form_section
_form_section_with_actions = _forms.form_section_with_actions
_inline_table = _forms.inline_table
_read_only_line = _forms.read_only_line
_read_only_text = _forms.read_only_text
_table_section = _forms.table_section
_text_section = _forms.text_section
_untitled_form_section = _forms.untitled_form_section


def _configure_calculation_progress_dialog(
    progress: QProgressDialog,
    title: str,
) -> None:
    """Configure calculation progress without blocking workspace navigation."""
    progress.setWindowTitle(title)
    progress.setWindowModality(Qt.WindowModality.NonModal)
    progress.setAutoClose(False)
    progress.setAutoReset(False)
    progress.setMinimumDuration(0)


def _queued_calculation_label(position: int, calculation_name: str) -> str:
    """Return the queued calculation progress-dialog label."""
    return f"(Queue #{position}) {calculation_name}: calculation queued."


class CalculationRunnerMixin:
    """Calculation execution, phase selection, and result dialog methods."""

    def _calculation_verbose_logging_enabled(self) -> bool:
        """Return whether detailed calculation context should be logged."""
        explicit = bool(getattr(self, "calculation_verbose_logging", False))
        env_value = os.environ.get("EQUILIPY_VERBOSE_CALCULATION_LOG", "")
        return explicit or env_value.strip().lower() in {"1", "true", "yes", "on"}

    def _append_calculation_log_message(
        self,
        message: str,
        *,
        verbose: bool = False,
    ) -> None:
        """Append a calculation event to the GUI log and active log file."""
        if verbose and not self._calculation_verbose_logging_enabled():
            return
        if not message:
            return
        if self.calculation_log is not None:
            self.calculation_log.setPlainText(message)
            return
        if self.log_capture is not None:
            self._write_status_to_log_capture(message)

    def _log_calculation_start(
        self,
        *,
        database: CalculationDatabase,
        module: CalculationModule | None,
        module_kind: str,
        units: list[str],
        phases: list[str] | None,
        phase_names: list[str],
        condition: dict[str, Any] | None,
    ) -> None:
        """Write a concise start record for one calculation."""
        module_name = module.name if module is not None else module_kind.title()
        selected_phase_summary = _calculation_phase_selection_summary(
            phases,
            phase_names,
        )
        point_count = _calculation_condition_point_count(condition)
        self._append_calculation_log_message(
            f"Starting {module_kind} calculation: {module_name} "
            f"using {database.name}; {selected_phase_summary}; "
            f"{point_count} point(s)."
        )
        self._append_calculation_log_message(
            f"Units: T={units[0]}, P={units[1]}, amount={units[2]}.",
            verbose=True,
        )
        if phases:
            self._append_calculation_log_message(
                f"Selected phases: {', '.join(phases)}.",
                verbose=True,
            )

    def _log_calculation_finish(
        self,
        *,
        module: CalculationModule | None,
        result: Any,
    ) -> None:
        """Write a concise completion record for one calculation."""
        module_name = module.name if module is not None else "Calculation"
        summary = _calculation_result_summary(result)
        self._append_calculation_log_message(
            f"Finished {module_name}. {summary}"
        )

    def _log_calculation_canceled(self) -> None:
        """Write a cancellation record."""
        self._append_calculation_log_message("Calculation canceled.")

    def _set_calculation_running_status(self, status: QPlainTextEdit) -> None:
        """Reset per-run status metadata and show the running state."""
        status.setProperty("_equilipy_cancel_logged", False)
        status.setProperty("_equilipy_calculation_warnings", [])
        status.setPlainText("Calculation running...")

    def _append_calculation_status_warning(
        self,
        status: QPlainTextEdit,
        message: str,
    ) -> None:
        """Record a warning in both the calculation log and module status."""
        warnings = list(status.property("_equilipy_calculation_warnings") or [])
        warnings.append(message)
        status.setProperty("_equilipy_calculation_warnings", warnings)
        self._append_calculation_log_message(message)
        status.appendPlainText(message)

    def _append_result_warnings(self, status: QPlainTextEdit, result: Any) -> None:
        """Persist warnings carried by a calculation result object."""
        for message in list(getattr(result, "warnings", []) or []):
            self._append_calculation_status_warning(status, str(message))

    def _set_calculation_finished_status(self, status: QPlainTextEdit) -> None:
        """Set finished status while preserving visible per-run warnings."""
        warnings = list(status.property("_equilipy_calculation_warnings") or [])
        if warnings:
            status.setPlainText("\n".join([*warnings, "Calculation finished."]))
        else:
            status.setPlainText("Calculation finished.")

    def _set_calculation_canceled_status(self, status: QPlainTextEdit) -> None:
        """Set cancellation status and record it in the calculation log."""
        if not bool(status.property("_equilipy_cancel_logged")):
            self._log_calculation_canceled()
            status.setProperty("_equilipy_cancel_logged", True)
        status.setPlainText("Calculation canceled.")

    def _set_calculation_queued_status(
        self,
        status: QPlainTextEdit,
        position: int,
        module_name: str = "Calculation",
    ) -> None:
        """Set queued status for a calculation waiting on the global solver lane."""
        message = f"(Queue #{position}) {module_name}: calculation queued."
        if status is getattr(self, "calculation_log", None):
            status.setPlainText(message)
        else:
            self._append_calculation_log_message(message)
            status.setPlainText(message)

    def _set_calculation_queue_canceled_status(
        self,
        status: QPlainTextEdit,
        module_name: str = "Calculation",
    ) -> None:
        """Set canceled status for a calculation removed from the queue."""
        message = f"{module_name}: queued calculation canceled."
        if status is getattr(self, "calculation_log", None):
            status.setPlainText(message)
        else:
            self._append_calculation_log_message(message)
            status.setPlainText(message)

    def _enqueue_calculation(
        self,
        callback,
        queued_callback=None,
        *,
        calculation_name: str = "Calculation",
        canceled_callback=None,
    ) -> int:
        """Append one GUI calculation callback to the global solver queue."""
        queue = getattr(self, "_calculation_queue", None)
        if queue is None:
            self._calculation_queue = []
            queue = self._calculation_queue
        position = len(queue) + 1
        label = _queued_calculation_label(position, calculation_name)
        parent = self if isinstance(self, QWidget) else None
        progress = QProgressDialog(label, "Cancel", 0, 1, parent)
        _configure_calculation_progress_dialog(progress, "Queued Calculation")
        progress.setValue(0)
        cancel_button = QPushButton("Cancel")
        cancel_button.setObjectName("SmallActionButton")
        cancel_button.setCursor(Qt.CursorShape.PointingHandCursor)
        progress.setCancelButton(cancel_button)
        entry: dict[str, Any] = {
            "callback": callback,
            "name": calculation_name,
            "progress": progress,
            "canceled": False,
            "starting": False,
        }

        def request_cancel() -> None:
            if bool(entry.get("starting")):
                return
            if bool(entry.get("canceled")):
                return
            entry["canceled"] = True
            try:
                queue.remove(entry)
            except ValueError:
                pass
            progress.close()
            progress.deleteLater()
            if canceled_callback is not None:
                canceled_callback()
            self._refresh_queued_calculation_dialogs()
            self._set_status_bar_message(
                f"{calculation_name}: queued calculation canceled.",
                level="warning",
            )

        entry["cancel"] = request_cancel
        cancel_button.clicked.connect(request_cancel)
        progress.canceled.connect(request_cancel)
        queue.append(entry)
        progress.show()
        if queued_callback is not None:
            queued_callback(position)
        self._set_status_bar_message(
            f"Calculation queued. {position} waiting.",
            level="warning",
        )
        return position

    def _refresh_queued_calculation_dialogs(self) -> None:
        """Refresh queue labels after starting or canceling queued calculations."""
        queue = getattr(self, "_calculation_queue", [])
        for index, entry in enumerate(queue, start=1):
            progress = entry.get("progress")
            if isinstance(progress, QProgressDialog):
                progress.setLabelText(
                    _queued_calculation_label(
                        index,
                        str(entry.get("name") or "Calculation"),
                    )
                )

    def _start_next_queued_calculation(self) -> None:
        """Start the next queued calculation once the global solver lane is free."""
        if self._calculation_in_progress:
            return
        queue = getattr(self, "_calculation_queue", [])
        if not queue:
            return
        entry = queue.pop(0)
        if bool(entry.get("canceled")):
            self._start_next_queued_calculation()
            return
        entry["starting"] = True
        self._refresh_queued_calculation_dialogs()
        progress = entry.get("progress")
        if isinstance(progress, QProgressDialog):
            progress.close()
            progress.deleteLater()
        queued_calculation = entry["callback"]
        remaining = len(queue)
        if remaining:
            self._set_status_bar_message(
                f"Starting queued calculation. {remaining} waiting.",
            )
        else:
            self._set_status_bar_message("Starting queued calculation.")
        queued_calculation()
        if not self._calculation_in_progress:
            self._start_next_queued_calculation()

    def _calculate_all_modules(self) -> None:
        session, _module = self._current_calculation_context()
        if self._calculation_in_progress:
            QMessageBox.information(
                self,
                "Calculate All",
                "A calculation is already running.",
            )
            return
        if session is None:
            QMessageBox.information(
                self,
                "Calculate All",
                "Create or select a calculation session first.",
            )
            return
        if not session.modules:
            QMessageBox.information(
                self,
                "Calculate All",
                "No calculation modules are in this session yet.",
            )
            return
        if self._active_session_database(session) is None:
            QMessageBox.warning(
                self,
                "Calculate All",
                "Load and select a session database before recalculating modules.",
            )
            return

        runnable_modules = [
            module for module in session.modules.values() if module.kind != "nucleation"
        ]
        skipped_modules = len(session.modules) - len(runnable_modules)
        if not runnable_modules:
            QMessageBox.information(
                self,
                "Calculate All",
                "No implemented calculation modules are in this session yet.",
            )
            return

        failed_modules: list[str] = []
        module_units = 1000
        total_units = max(1, len(runnable_modules) * module_units)
        progress = QProgressDialog(
            "Preparing Calculate All...",
            "Cancel",
            0,
            total_units,
            self,
        )
        _configure_calculation_progress_dialog(progress, "Calculate All")
        progress.setValue(0)

        cancel_requested = False
        cancel_button = QPushButton("Cancel")
        cancel_button.setObjectName("SmallActionButton")
        cancel_button.setCursor(Qt.CursorShape.PointingHandCursor)
        progress.setCancelButton(cancel_button)

        def request_cancel() -> None:
            nonlocal cancel_requested
            if cancel_requested:
                return
            cancel_requested = True
            progress.setLabelText("Canceling after current solver step...")
            progress.setCancelButton(None)

        cancel_button.clicked.connect(request_cancel)
        progress.canceled.connect(request_cancel)
        progress.show()
        QApplication.processEvents()
        self._calculation_in_progress = True

        def module_progress_callback(module_index: int, module_name: str):
            def update(current: int, total: int, message: str = "") -> None:
                if cancel_requested:
                    raise CalculationCanceled("Calculation canceled.")
                base = module_index * module_units
                if total > 0:
                    current_value = max(0, min(int(current), int(total)))
                    module_value = int(module_units * current_value / int(total))
                    overall_value = min(total_units, base + module_value)
                else:
                    overall_value = min(total_units, base)
                percent = 100.0 * overall_value / total_units
                detail = message or "Running"
                progress.setLabelText(
                    f"{module_name} - {detail}: {percent:.2f}%"
                )
                progress.setValue(overall_value)
                QApplication.processEvents()

            return update

        canceled = False
        try:
            for module_index, module in enumerate(runnable_modules):
                if cancel_requested:
                    raise CalculationCanceled("Calculation canceled.")
                runtime = self._ensure_calculation_module_runtime(session, module)
                progress.setLabelText(f"{module.name} - running...")
                progress.setValue(module_index * module_units)
                QApplication.processEvents()
                if cancel_requested:
                    raise CalculationCanceled("Calculation canceled.")
                if not self._run_calculation(
                    session.id,
                    module.id,
                    module.kind,
                    runtime["composition_table"],
                    runtime["phase_tree"],
                    runtime["temperature"],
                    runtime["pressure"],
                    runtime["delta_t"],
                    runtime["liquid_phase"],
                    runtime["results_table"],
                    runtime["status"],
                    show_dialog=False,
                    progress_callback=module_progress_callback(
                        module_index,
                        module.name,
                    ),
                    lane_owned=True,
                ):
                    failed_modules.append(module.name)
                progress.setValue((module_index + 1) * module_units)
                QApplication.processEvents()
        except CalculationCanceled:
            canceled = True
        finally:
            progress.close()
            self._calculation_in_progress = False
            self._start_next_queued_calculation()

        if canceled:
            QMessageBox.information(
                self,
                "Calculate All",
                "Calculation canceled.",
            )
            return

        self._show_calculation_overview(session)
        skipped_text = (
            f" Skipped {skipped_modules} nucleation module(s)."
            if skipped_modules
            else ""
        )
        if failed_modules:
            QMessageBox.warning(
                self,
                "Calculate All",
                f"Calculated {len(runnable_modules) - len(failed_modules)} of "
                f"{len(runnable_modules)} implemented module(s). Check Log "
                f"for details.{skipped_text}",
            )
        else:
            QMessageBox.information(
                self,
                "Calculate All",
                f"Calculated {len(runnable_modules)} implemented module(s)."
                f"{skipped_text}",
            )

    def _confirm_calculation_phases(
        self,
        session_id: str,
        module_id: str,
        composition_table: QTableWidget,
        phase_tree: QTreeWidget,
        status: QPlainTextEdit,
        *,
        lane_owned: bool = False,
    ) -> bool:
        session = self.calculation_sessions[session_id]
        active_database = self._active_session_database(session)
        if active_database is None:
            status.setPlainText(
                "Load and select a session database before phase lookup."
            )
            return False
        module = session.modules.get(module_id)
        try:
            elements = _module_condition_elements(
                module,
                composition_table,
                active_database.database,
            )
        except ValueError as exc:
            status.setPlainText(f"Condition issue: {exc}")
            return False
        if (
            module is None
            or module.kind != "equilibrium"
            or module.calculation_type != "batch"
        ):
            if not self._validate_composition_table(
                session,
                composition_table,
                status,
            ):
                return False
        if not elements:
            status.setPlainText("Enter at least one element in the composition table.")
            return False
        if self._calculation_in_progress and not lane_owned:
            # list_phases drives the same process-global Fortran state the
            # calculation worker thread is using; a concurrent call corrupts
            # it and can crash the whole application.
            status.setPlainText(
                "Phase lookup is paused while a calculation is running; "
                "try again after it finishes."
            )
            return False
        try:
            import equilipy as eq

            phase_names = [
                phase
                for phase in eq.list_phases(active_database.database, elements)
                if not _is_reference_ser_phase(phase)
            ]
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            status.setPlainText(f"Phase lookup failed: {exc}")
            return False
        if not phase_names:
            status.setPlainText("No selectable phases were found for this condition.")
            return False
        if module is not None:
            module.phase_names = phase_names
            if module.kind == "solidification":
                liquid_widget = module.runtime.get("liquid_phase")
                if isinstance(liquid_widget, QLineEdit):
                    liquid_name = _resolve_liquid_phase_name(
                        liquid_widget.text(),
                        phase_names,
                    )
                    if liquid_name:
                        liquid_widget.setText(liquid_name)
        default_selection = _default_phase_selection_for_module(
            module,
            phase_names,
            session.default_phase_names,
            active_database.database,
        )
        _populate_phase_tree(phase_tree, phase_names, default_selection)
        if module is not None:
            self._refresh_nucleation_undercooling_runtime(module)
        self._update_current_calculation_heading(session)
        status.setPlainText(
            f"Loaded {len(phase_names)} phase(s) from {active_database.name}."
        )
        excluded_ordered = _default_scheil_excluded_ordered_phases_for_module(
            module,
            phase_names,
            session.default_phase_names,
            active_database.database,
        )
        if excluded_ordered:
            self._append_calculation_status_warning(
                status,
                scheil_ordered_phase_exclusion_notice(excluded_ordered),
            )
        return True

    def _ensure_calculation_phases_current(
        self,
        session_id: str,
        module_id: str,
        composition_table: QTableWidget,
        phase_tree: QTreeWidget,
        status: QPlainTextEdit,
        *,
        lane_owned: bool = False,
    ) -> bool:
        session = self.calculation_sessions[session_id]
        module = session.modules.get(module_id)
        active_database = self._active_session_database(session)
        if module is None or active_database is None:
            return self._confirm_calculation_phases(
                session_id,
                module_id,
                composition_table,
                phase_tree,
                status,
                lane_owned=lane_owned,
            )

        try:
            elements = _module_condition_elements(
                module,
                composition_table,
                active_database.database,
            )
        except ValueError:
            return self._confirm_calculation_phases(
                session_id,
                module_id,
                composition_table,
                phase_tree,
                status,
                lane_owned=lane_owned,
            )
        if not elements:
            return self._confirm_calculation_phases(
                session_id,
                module_id,
                composition_table,
                phase_tree,
                status,
                lane_owned=lane_owned,
            )

        if self._calculation_in_progress and not lane_owned:
            # Never touch the process-global Fortran state from the main
            # thread while a calculation worker owns it.
            return self._confirm_calculation_phases(
                session_id,
                module_id,
                composition_table,
                phase_tree,
                status,
            )

        try:
            import equilipy as eq

            phase_names = [
                phase
                for phase in eq.list_phases(active_database.database, elements)
                if not _is_reference_ser_phase(phase)
            ]
        except Exception:
            return self._confirm_calculation_phases(
                session_id,
                module_id,
                composition_table,
                phase_tree,
                status,
                lane_owned=lane_owned,
            )

        current_phase_names = _all_phase_names(phase_tree) or list(module.phase_names)
        if current_phase_names and set(current_phase_names) == set(phase_names):
            return True

        current_selection = _selected_phase_names(phase_tree)
        preferred_selection = (
            current_selection
            or _default_phase_selection_for_module(
                module,
                phase_names,
                session.default_phase_names,
                active_database.database,
            )
        )
        selected = [phase for phase in preferred_selection if phase in phase_names]
        module.phase_names = phase_names
        if module.kind == "solidification":
            liquid_widget = module.runtime.get("liquid_phase")
            if isinstance(liquid_widget, QLineEdit):
                liquid_name = _resolve_liquid_phase_name(
                    liquid_widget.text(),
                    phase_names,
                )
                if liquid_name:
                    liquid_widget.setText(liquid_name)
        _populate_phase_tree(phase_tree, phase_names, selected or None)
        self._refresh_nucleation_undercooling_runtime(module)
        self._update_current_calculation_heading(session)
        status.setPlainText(
            f"Confirmed {len(phase_names)} phase(s) from {active_database.name}."
        )
        if not current_selection:
            excluded_ordered = _default_scheil_excluded_ordered_phases_for_module(
                module,
                phase_names,
                session.default_phase_names,
                active_database.database,
            )
            if excluded_ordered:
                self._append_calculation_status_warning(
                    status,
                    scheil_ordered_phase_exclusion_notice(excluded_ordered),
                )
        return True

    def _open_module_phase_dialog(
        self,
        session_id: str,
        module_id: str,
        composition_table: QTableWidget,
        phase_tree: QTreeWidget,
        status: QPlainTextEdit,
    ) -> None:
        session = self.calculation_sessions.get(session_id)
        module = session.modules.get(module_id) if session is not None else None
        if session is None or module is None:
            return

        if not module.phase_names:
            self._confirm_calculation_phases(
                session_id,
                module_id,
                composition_table,
                phase_tree,
                status,
            )
        if not module.phase_names:
            return

        selected_phases = _selected_phase_names(phase_tree)
        if not selected_phases:
            selected_phases = _default_phase_selection_for_module(
                module,
                module.phase_names,
                session.default_phase_names,
            )

        dialog = QDialog(self)
        dialog.setWindowTitle(f"{module.name} Phase Selection")
        dialog.setModal(True)
        dialog.resize(620, 520)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(12)

        title = QLabel(module.name)
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        dialog_phase_tree = QTreeWidget()
        dialog_phase_tree.setObjectName("DatabaseTree")
        dialog_phase_tree.setHeaderHidden(True)
        dialog_phase_tree.setColumnCount(3)
        dialog_phase_tree.header().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        _populate_phase_tree(
            dialog_phase_tree,
            module.phase_names,
            selected_phases or None,
            fit_to_contents=False,
        )

        tool_row = QHBoxLayout()
        tool_row.setSpacing(8)
        all_solutions = QPushButton("All Solutions")
        all_solutions.setObjectName("SmallActionButton")
        all_solutions.clicked.connect(
            lambda _checked=False: _toggle_phase_group(
                dialog_phase_tree,
                "Solutions",
            )
        )
        all_compounds = QPushButton("All Compounds")
        all_compounds.setObjectName("SmallActionButton")
        all_compounds.clicked.connect(
            lambda _checked=False: _toggle_phase_group(
                dialog_phase_tree,
                "Compounds",
            )
        )
        tool_row.addWidget(all_solutions)
        tool_row.addWidget(all_compounds)
        tool_row.addStretch(1)
        layout.addLayout(tool_row)
        layout.addWidget(dialog_phase_tree, 1)

        actions = QHBoxLayout()
        actions.addStretch(1)
        apply = QPushButton("Apply")
        apply.setObjectName("PrimaryButton")
        cancel = QPushButton("Cancel")
        cancel.setObjectName("SmallActionButton")
        cancel.setDefault(True)
        cancel.setAutoDefault(True)
        apply.setAutoDefault(False)

        def apply_selection() -> None:
            selected = _selected_phase_names(dialog_phase_tree)
            if not selected:
                QMessageBox.warning(
                    dialog,
                    "Phase Selection",
                    "Select at least one phase.",
                )
                return
            _populate_phase_tree(phase_tree, module.phase_names, selected)
            if module.kind == "solidification":
                liquid_widget = module.runtime.get("liquid_phase")
                if isinstance(liquid_widget, QLineEdit):
                    liquid_name = _resolve_liquid_phase_name(
                        liquid_widget.text(),
                        module.phase_names,
                    )
                    if liquid_name:
                        liquid_widget.setText(liquid_name)
                self._refresh_nucleation_undercooling_runtime(module)
            self._update_current_calculation_heading(session)
            status.setPlainText(
                f"Selected {len(selected)} of {len(module.phase_names)} phase(s)."
            )
            dialog.accept()

        apply.clicked.connect(apply_selection)
        cancel.clicked.connect(dialog.reject)
        actions.addWidget(apply)
        actions.addWidget(cancel)
        layout.addLayout(actions)
        dialog.exec()

    def _open_result_columns_dialog(
        self,
        module: CalculationModule,
        table: QTableWidget,
    ) -> None:
        result = module.runtime.get("result") or module.result
        result_table = _result_table_for_module(module, result)
        if result_table is None:
            QMessageBox.information(
                self,
                "Result Columns",
                "Run or load a calculation result before selecting columns.",
            )
            return

        dialog = QDialog(self)
        dialog.setWindowTitle(f"{module.name} Result Columns")
        dialog.setModal(True)
        dialog.resize(560, 560)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(12)

        session_id = str(module.runtime.get("session_id") or "")
        session = self.calculation_sessions.get(session_id)
        system_default_columns = set(
            _default_result_column_keys(module, result_table, result)
        )
        saved_default_columns = set(
            _expand_result_column_defaults(
                list(session.result_column_defaults.get(module.kind, []))
                if session is not None
                else [],
                result_table,
                module,
            )
        )
        default_columns = (
            saved_default_columns & set(result_table.available_columns())
        ) or system_default_columns
        current_columns = (
            set(
                _expand_result_column_defaults(
                    list(module.result_columns),
                    result_table,
                    module,
                )
            )
            if module.result_columns
            else set(default_columns)
        )
        unit_filters = _result_unit_filters_from_columns(
            current_columns,
            module,
            result_table.available_columns(),
        )
        checkboxes: list[tuple[str, QCheckBox]] = []
        checkbox_by_key: dict[str, QCheckBox] = {}
        sync_callbacks: list[Any] = []
        syncing_section_checks = False

        scroll = QScrollArea(dialog)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        checkbox_host = QWidget()
        checkbox_layout = QVBoxLayout(checkbox_host)
        checkbox_layout.setContentsMargins(0, 0, 0, 0)
        checkbox_layout.setSpacing(6)

        columns_by_group: dict[str, list[Any]] = {}
        for column in result_table.columns:
            columns_by_group.setdefault(column.group, []).append(column)

        def set_checkboxes(
            target_checkboxes: list[QCheckBox],
            checked: bool,
        ) -> None:
            for checkbox in target_checkboxes:
                checkbox.setChecked(checked)

        def set_section_state(
            checkbox: QCheckBox,
            target_checkboxes: list[QCheckBox],
        ) -> None:
            nonlocal syncing_section_checks
            checked_count = sum(1 for target in target_checkboxes if target.isChecked())
            syncing_section_checks = True
            checkbox.setTristate(True)
            if checked_count == 0:
                checkbox.setCheckState(Qt.CheckState.Unchecked)
            elif checked_count == len(target_checkboxes):
                checkbox.setCheckState(Qt.CheckState.Checked)
            else:
                checkbox.setCheckState(Qt.CheckState.PartiallyChecked)
            syncing_section_checks = False

        scroll.setWidget(checkbox_host)
        layout.addWidget(scroll)

        def sync_all_headers() -> None:
            for sync in sync_callbacks:
                sync()

        def clear_checkbox_layout() -> None:
            while checkbox_layout.count():
                item = checkbox_layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.deleteLater()

        def visible_column_keys() -> set[str]:
            return {
                column.key
                for column in result_table.columns
                if _result_column_visible_for_unit_filters(column.key, unit_filters)
            }

        def rebuild_column_groups() -> None:
            nonlocal checkboxes, checkbox_by_key, sync_callbacks, current_columns
            if checkboxes:
                current_columns = {
                    key for key, checkbox in checkboxes if checkbox.isChecked()
                }
            visible_keys = visible_column_keys()
            current_columns &= visible_keys
            checkboxes = []
            checkbox_by_key = {}
            sync_callbacks = []
            clear_checkbox_layout()

            columns_by_group = {}
            for column in result_table.columns:
                if column.key not in visible_keys:
                    continue
                if _result_column_is_stability(column.key):
                    continue
                columns_by_group.setdefault(column.group, []).append(column)

            for group, columns in columns_by_group.items():
                group_frame = QFrame()
                group_frame.setObjectName("ResultColumnGroup")
                group_layout = QVBoxLayout(group_frame)
                group_layout.setContentsMargins(0, 0, 0, 0)
                group_layout.setSpacing(6)

                group_header = QHBoxLayout()
                group_header.setSpacing(6)
                group_toggle = QToolButton()
                group_toggle.setObjectName("ResultColumnGroupToggle")
                group_toggle.setText(group.replace("_", " ").title())
                group_toggle.setToolButtonStyle(
                    Qt.ToolButtonStyle.ToolButtonTextBesideIcon
                )
                group_toggle.setArrowType(Qt.ArrowType.DownArrow)
                group_toggle.setCheckable(True)
                group_toggle.setChecked(True)
                group_header.addWidget(group_toggle)

                section_box = QCheckBox("")
                section_box.setObjectName("ResultSectionCheck")
                section_box.setToolTip("Select all visible columns in this section")
                group_header.addWidget(section_box)
                group_header.addStretch(1)

                content = QWidget()
                content_layout = QVBoxLayout(content)
                content_layout.setContentsMargins(20, 0, 0, 0)
                content_layout.setSpacing(4)

                group_checkboxes: list[QCheckBox] = []
                for column in columns:
                    checkbox = QCheckBox(column.label or column.key)
                    checkbox.setObjectName("ResultColumnCheck")
                    checkbox.setToolTip(column.key)
                    checkbox.setChecked(column.key in current_columns)
                    content_layout.addWidget(checkbox)
                    checkboxes.append((column.key, checkbox))
                    checkbox_by_key[column.key] = checkbox
                    group_checkboxes.append(checkbox)

                def sync_header_checks(
                    all_checkbox: QCheckBox = section_box,
                    all_targets: list[QCheckBox] = group_checkboxes,
                ) -> None:
                    set_section_state(all_checkbox, all_targets)

                def on_section_state_changed(
                    state: int,
                    targets: list[QCheckBox] = group_checkboxes,
                    sync: Any = sync_header_checks,
                ) -> None:
                    if syncing_section_checks:
                        return
                    set_checkboxes(targets, state != 0)
                    sync()

                section_box.stateChanged.connect(on_section_state_changed)
                for checkbox in group_checkboxes:
                    checkbox.stateChanged.connect(
                        lambda _state, sync=sync_header_checks: sync()
                    )
                sync_callbacks.append(sync_header_checks)

                group_toggle.toggled.connect(
                    lambda checked, widget=content, toggle=group_toggle: (
                        widget.setVisible(checked),
                        toggle.setArrowType(
                            Qt.ArrowType.DownArrow
                            if checked
                            else Qt.ArrowType.RightArrow
                        ),
                    )
                )
                sync_header_checks()
                group_layout.addLayout(group_header)
                group_layout.addWidget(content)
                checkbox_layout.addWidget(group_frame)

            checkbox_layout.addStretch(1)
            for _key, checkbox in checkboxes:
                checkbox.stateChanged.connect(lambda _state: sync_all_headers())
            sync_all_headers()

        def set_selected_columns(keys: set[str]) -> None:
            for key, checkbox in checkboxes:
                checkbox.setChecked(key in keys)
            sync_all_headers()

        bulk_actions = QHBoxLayout()
        output_unit_button = QPushButton("Output Unit")
        output_unit_button.setObjectName("SmallActionButton")
        select_button = QToolButton()
        select_button.setText("Select")
        select_button.setObjectName("SmallMenuButton")
        select_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        select_menu = QMenu(select_button)
        select_all_action = QAction("Select All", select_menu)
        select_stable_action = QAction("Select Stable", select_menu)
        select_menu.addAction(select_all_action)
        select_menu.addAction(select_stable_action)
        select_button.setMenu(select_menu)
        clear_all = QPushButton("Clear")
        clear_all.setObjectName("SmallActionButton")
        save_default = QPushButton("Save as Default")
        save_default.setObjectName("SmallActionButton")
        reset_default = QPushButton("Default")
        reset_default.setObjectName("SmallActionButton")
        select_all_action.triggered.connect(
            lambda: set_selected_columns({key for key, _checkbox in checkboxes})
        )
        select_stable_action.triggered.connect(
            lambda: set_selected_columns(
                set(_stable_result_column_keys(result_table, result))
                & visible_column_keys()
            )
        )
        clear_all.clicked.connect(
            lambda: set_selected_columns(set())
        )
        output_unit_button.clicked.connect(
            lambda: self._open_output_unit_dialog(
                result_table.available_columns(),
                unit_filters,
                lambda filters: (
                    unit_filters.clear(),
                    unit_filters.update(filters),
                    rebuild_column_groups(),
                ),
            )
        )
        bulk_actions.addWidget(output_unit_button)
        bulk_actions.addWidget(select_button)
        bulk_actions.addWidget(clear_all)
        bulk_actions.addStretch(1)
        bulk_actions.addWidget(save_default)
        bulk_actions.addWidget(reset_default)
        layout.addLayout(bulk_actions)
        rebuild_column_groups()

        actions = QHBoxLayout()
        actions.addStretch(1)
        apply = QPushButton("Apply")
        apply.setObjectName("PrimaryButton")
        cancel = QPushButton("Cancel")
        cancel.setObjectName("SmallActionButton")
        actions.addWidget(apply)
        actions.addWidget(cancel)
        layout.addLayout(actions)

        def apply_selection() -> None:
            selected = [key for key, checkbox in checkboxes if checkbox.isChecked()]
            if not selected:
                QMessageBox.warning(
                    dialog,
                    "Result Columns",
                    "Select at least one result column.",
                )
                return
            module.result_columns = selected
            _populate_module_results(module, table, result)
            module.results_table_payload = _table_payload(table)
            dialog.accept()

        def save_selection_as_default() -> None:
            nonlocal default_columns
            selected = [key for key, checkbox in checkboxes if checkbox.isChecked()]
            if not selected:
                QMessageBox.warning(
                    dialog,
                    "Result Columns",
                    "Select at least one result column before saving a default.",
                )
                return
            if session is None:
                return
            saved_defaults = _generalized_result_column_defaults(
                selected,
                result_table,
                unit_filters,
            )
            session.result_column_defaults[module.kind] = saved_defaults
            default_columns = set(
                _expand_result_column_defaults(
                    saved_defaults,
                    result_table,
                    module,
                )
            )

        apply.clicked.connect(apply_selection)
        save_default.clicked.connect(save_selection_as_default)
        reset_default.clicked.connect(lambda: set_selected_columns(default_columns))
        cancel.clicked.connect(dialog.reject)
        dialog.exec()

    def _show_all_module_results(
        self,
        module: CalculationModule,
        table: QTableWidget,
    ) -> None:
        """Render every stored result row into the GUI table on request."""
        result = module.runtime.get("result") or module.result
        if result is None:
            QMessageBox.information(
                self,
                "Show All Results",
                "Run or load a calculation result before showing all rows.",
            )
            return

        table.setProperty("show_all_results", True)
        _populate_module_results(module, table, result)
        module.results_table_payload = _table_payload(table)

    def _open_output_unit_dialog(
        self,
        available_columns: list[str],
        current_filters: set[str],
        apply_callback: Any,
    ) -> None:
        options = _result_output_unit_filter_options(available_columns)
        if not options:
            QMessageBox.information(
                self,
                "Output Unit",
                "No alternate output units are available for this result.",
            )
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Output Unit")
        dialog.setModal(True)
        dialog.resize(420, 360)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(14)
        boxes: list[tuple[str, QCheckBox]] = []

        for section, section_options in options.items():
            if not section_options:
                continue
            section_label = QLabel(section)
            section_label.setObjectName("SmallSectionLabel")
            layout.addWidget(section_label)
            for label, filter_key in section_options:
                checkbox = QCheckBox(label)
                checkbox.setObjectName("ResultUnitCheck")
                checkbox.setChecked(filter_key in current_filters)
                layout.addWidget(checkbox)
                boxes.append((filter_key, checkbox))

        layout.addStretch(1)
        actions = QHBoxLayout()
        actions.addStretch(1)
        apply = QPushButton("Apply")
        apply.setObjectName("PrimaryButton")
        cancel = QPushButton("Cancel")
        cancel.setObjectName("SmallActionButton")
        actions.addWidget(apply)
        actions.addWidget(cancel)
        layout.addLayout(actions)

        def apply_units() -> None:
            selected = {key for key, checkbox in boxes if checkbox.isChecked()}
            apply_callback(selected)
            dialog.accept()

        apply.clicked.connect(apply_units)
        cancel.clicked.connect(dialog.reject)
        dialog.exec()

    def _run_with_progress_dialog(
        self,
        work,
        finished_callback,
        failure_callback,
        *,
        title: str,
        message: str,
        canceled_callback=None,
        lane_owned: bool = False,
    ) -> bool:
        """Run work in a thread and show progress only if it takes a moment."""
        if self._calculation_in_progress and not lane_owned:
            QMessageBox.information(
                self,
                title,
                "A calculation is already running.",
            )
            return False

        if not lane_owned:
            self._calculation_in_progress = True
        thread = QThread(self)
        worker = CalculationWorker(work)
        receiver = CalculationResultReceiver(self)
        worker.moveToThread(thread)
        self._calculation_thread_refs.append((thread, worker, receiver))

        progress = QProgressDialog(message, "Cancel", 0, 0, self)
        _configure_calculation_progress_dialog(progress, title)
        progress.reset()

        cancel_button = QPushButton("Cancel")
        cancel_button.setObjectName("SmallActionButton")
        cancel_button.setCursor(Qt.CursorShape.PointingHandCursor)
        progress.setCancelButton(cancel_button)

        cancel_requested = False

        def request_cancel() -> None:
            nonlocal cancel_requested
            if cancel_requested:
                return
            cancel_requested = True
            worker.cancel()
            progress.setLabelText("Canceling after current solver step...")
            self._set_status_bar_message(
                "Canceling after current solver step...",
                level="warning",
            )
            progress.setCancelButton(None)

        cancel_button.clicked.connect(request_cancel)
        progress.canceled.connect(request_cancel)

        timer = QTimer(self)
        timer.setSingleShot(True)
        timer.timeout.connect(progress.show)

        def handle_progress(current: int, total: int, progress_message: str) -> None:
            if total > 0:
                current_value = max(0, min(current, total))
                percent = 100.0 * current_value / total
                if progress_message:
                    status_message = f"{progress_message}: {percent:.2f}%"
                else:
                    status_message = f"{percent:.2f}%"
                progress.setLabelText(status_message)
                progress.setRange(0, total)
                progress.setValue(current_value)
            else:
                if progress_message:
                    status_message = progress_message
                    progress.setLabelText(status_message)
                else:
                    status_message = "Working..."
                progress.setRange(0, 0)
            self._set_status_bar_message(status_message)

        def cleanup_thread_refs() -> None:
            timer.stop()
            progress.close()
            try:
                self._calculation_thread_refs.remove((thread, worker, receiver))
            except ValueError:
                pass
            timer.deleteLater()
            receiver.deleteLater()

        def handle_finished(result, error) -> None:
            timer.stop()
            progress.close()
            try:
                if error is not None:
                    if isinstance(error, CalculationCanceled):
                        if canceled_callback is not None:
                            canceled_callback()
                        QMessageBox.information(
                            self,
                            title,
                            "Calculation canceled.",
                        )
                        return
                    failure_callback(error)
                    return
                finished_callback(result)
            finally:
                if not lane_owned:
                    self._calculation_in_progress = False
                    self._start_next_queued_calculation()

        receiver.finished.connect(handle_finished)
        receiver.progress.connect(handle_progress)
        worker.progress.connect(
            receiver.receive_progress,
            Qt.ConnectionType.QueuedConnection,
        )
        worker.finished.connect(
            receiver.receive,
            Qt.ConnectionType.QueuedConnection,
        )
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(cleanup_thread_refs)
        thread.finished.connect(thread.deleteLater)
        thread.started.connect(worker.run)
        timer.start(1000)
        thread.start()
        return True

    def _run_calculation(
        self,
        session_id: str,
        module_id: str,
        module_kind: str,
        composition_table: QTableWidget,
        phase_tree: QTreeWidget,
        temperature: QLineEdit,
        pressure: QLineEdit,
        delta_t: QLineEdit,
        liquid_phase: QLineEdit,
        results_table: QTableWidget,
        status: QPlainTextEdit,
        *,
        show_dialog: bool = True,
        progress_callback=None,
        lane_owned: bool = False,
    ) -> bool:
        session = self.calculation_sessions[session_id]
        active_database = self._active_session_database(session)
        if active_database is None:
            self._set_calculation_failure(
                status,
                "Load and select a session database before running.",
                show_dialog=show_dialog,
            )
            return False

        module = session.modules.get(module_id)
        units = [
            module.temperature_unit if module is not None else "K",
            module.pressure_unit if module is not None else "atm",
            module.amount_unit if module is not None else "moles",
        ]
        is_equilibrium_module = module is not None and module_kind == "equilibrium"
        is_solidification_module = (
            module is not None and module_kind == "solidification"
        )
        has_imported_batch_conditions = (
            module is not None
            and module.kind in {"equilibrium", "solidification"}
            and module.calculation_type == "batch"
            and bool(module.batch_condition)
        )
        has_manual_batch_conditions = (
            module is not None
            and module.kind in {"equilibrium", "solidification"}
            and module.calculation_type == "batch"
            and not module.batch_condition
        )
        has_imported_equilibrium_conditions = (
            is_equilibrium_module and has_imported_batch_conditions
        )
        has_imported_solidification_conditions = (
            is_solidification_module and has_imported_batch_conditions
        )
        condition: dict[str, Any] | None = None
        if has_imported_batch_conditions:
            if not module.batch_condition:
                self._set_calculation_failure(
                    status,
                    "Import batch conditions before running.",
                    show_dialog=show_dialog,
                )
                return False
            condition = module.batch_condition
        elif is_equilibrium_module:
            if not self._validate_composition_table(
                session,
                composition_table,
                status,
            ):
                self._show_calculation_failure_dialog(
                    "Calculation Warning",
                    _latest_status_text(status)
                    or "Fix composition issues before running.",
                    show_dialog=show_dialog,
                )
                return False
            if has_manual_batch_conditions:
                try:
                    condition = _manual_batch_condition(
                        active_database.database,
                        composition_table,
                        temperature.text(),
                        pressure.text(),
                        units,
                    )
                except ValueError as exc:
                    self._set_calculation_failure(
                        status,
                        f"Invalid batch input: {exc}",
                        show_dialog=show_dialog,
                    )
                    return False
        else:
            if not self._validate_composition_table(
                session,
                composition_table,
                status,
            ):
                self._show_calculation_failure_dialog(
                    "Calculation Warning",
                    _latest_status_text(status)
                    or "Fix composition issues before running.",
                    show_dialog=show_dialog,
                )
                return False
            try:
                if has_manual_batch_conditions:
                    condition = _manual_batch_condition(
                        active_database.database,
                        composition_table,
                        temperature.text(),
                        pressure.text(),
                        units,
                    )
                else:
                    condition = _calculation_condition(
                        composition_table,
                        float(temperature.text()),
                        float(pressure.text()),
                        units[2],
                        active_database.database,
                    )
            except ValueError as exc:
                self._set_calculation_failure(
                    status,
                    f"Invalid calculation input: {exc}",
                    show_dialog=show_dialog,
                )
                return False

        if module is not None:
            if not self._ensure_calculation_phases_current(
                session_id,
                module_id,
                composition_table,
                phase_tree,
                status,
                lane_owned=lane_owned,
            ):
                self._show_calculation_failure_dialog(
                    "Calculation Warning",
                    _latest_status_text(status)
                    or "Condition could not be confirmed before running.",
                    show_dialog=show_dialog,
                )
                return False

        phase_names = _all_phase_names(phase_tree)
        if not phase_names and module is not None:
            phase_names = list(module.phase_names)
        selected_phases = _selected_phase_names(phase_tree)
        if phase_names and not selected_phases and module is not None:
            selected_phases = _default_phase_selection_for_module(
                module,
                phase_names,
                session.default_phase_names,
                active_database.database,
            )
        phases = selected_phases or None
        self._log_calculation_start(
            database=active_database,
            module=module,
            module_kind=module_kind,
            units=units,
            phases=phases,
            phase_names=phase_names,
            condition=condition,
        )
        self._set_calculation_running_status(status)
        try:
            import equilipy as eq

            if module_kind == "solidification":
                liquid_name = _resolve_liquid_phase_name(
                    liquid_phase.text(),
                    phase_names,
                )
                if not liquid_name:
                    liquid_name = liquid_phase.text().strip() or "LIQUID"
                liquid_phase.setText(liquid_name)
                if phases is not None and liquid_name not in phases:
                    phases = [liquid_name, *phases]
                if phases is not None and not any(
                    phase != liquid_name for phase in phases
                ):
                    self._set_calculation_failure(
                        status,
                        "Solidification requires at least one non-liquid phase. "
                        "Open Phase Selection and select solid phases.",
                        show_dialog=show_dialog,
                    )
                    return False
                start_from_liquidus_widget = (
                    module.runtime.get("from_liquidus") if module else None
                )
                start_from_liquidus = True
                if isinstance(start_from_liquidus_widget, QCheckBox):
                    start_from_liquidus = start_from_liquidus_widget.isChecked()
                if module is not None:
                    module.start_from_liquidus = start_from_liquidus
                solidification_model = _solidification_model_value(
                    module.solidification_model if module is not None else "scheil"
                )
                if solidification_model == "scheil":
                    initial_solid_phases = _initial_scheil_solid_phase_names(
                        eq,
                        active_database.database,
                        condition,
                        units,
                        phases,
                        liquid_name,
                    )
                    if initial_solid_phases:
                        phase_list = ", ".join(initial_solid_phases)
                        self._set_calculation_failure(
                            status,
                            "The initial Scheil condition already contains "
                            f"solid phase(s): {phase_list}. Increase the "
                            "starting temperature above the liquidus and run "
                            "the calculation again.",
                            title="Scheil Initial Temperature",
                            show_dialog=show_dialog,
                        )
                        return False
                critical_undercooling = {}
                if solidification_model == "nucleoscheil":
                    stored_undercooling = (
                        dict(module.nucleation_undercooling)
                        if module is not None
                        else {}
                    )
                    undercooling_tree = (
                        module.runtime.get("nucleation_undercooling_tree")
                        if module
                        else None
                    )
                    try:
                        undercooling_values = (
                            _nucleation_undercooling_tree_values(
                                undercooling_tree,
                                strict=True,
                            )
                            if isinstance(undercooling_tree, QTreeWidget)
                            else dict(module.nucleation_undercooling)
                        )
                    except ValueError as exc:
                        self._set_calculation_failure(
                            status,
                            f"Invalid nucleation undercooling: {exc}",
                            show_dialog=show_dialog,
                        )
                        return False
                    if module is not None:
                        module.nucleation_undercooling = dict(undercooling_values)
                    missing_undercooling = _missing_nucleation_undercooling_phases(
                        liquid_name,
                        phases,
                        stored_undercooling,
                    )
                    if missing_undercooling:
                        preview = ", ".join(missing_undercooling[:6])
                        if len(missing_undercooling) > 6:
                            preview = f"{preview}, ..."
                        self._append_calculation_status_warning(
                            status,
                            "Warning: no critical undercooling set for "
                            f"{len(missing_undercooling)} phase(s) "
                            f"({preview}); using 0.0 K.",
                        )
                    critical_undercooling = _nucleation_undercooling_for_phases(
                        liquid_name,
                        phases,
                        undercooling_values,
                    )
                delta_t_value = float(delta_t.text())
                solidification_batch = (
                    has_imported_solidification_conditions
                    or (is_solidification_module and has_manual_batch_conditions)
                )
                n_cpu = (
                    _effective_batch_cpu_count(condition, module.batch_cpu_count)
                    if solidification_batch
                    else 1
                )

                def execute_solver(progress_callback) -> tuple[Any, str | None]:
                    if solidification_batch:
                        progress_callback(0, 0, "Preparing solidification batch...")
                        def run_batch(batch_cpu: int) -> Any:
                            if solidification_model == "equilibrium_cooling":
                                return _equilib_cooling_batch(
                                    eq,
                                    liquid_name,
                                    active_database.database,
                                    condition,
                                    delta_T=delta_t_value,
                                    unit=units,
                                    phases=phases,
                                    start_from_liquidus=start_from_liquidus,
                                    progress_callback=progress_callback,
                                )
                            if solidification_model == "nucleoscheil":
                                return eq.nucleoscheil_batch(
                                    liquid_name,
                                    active_database.database,
                                    condition,
                                    critical_undercooling or None,
                                    delta_T=delta_t_value,
                                    unit=units,
                                    n_cpu=batch_cpu,
                                    phases=phases,
                                    progress_callback=progress_callback,
                                )
                            return eq.scheil_batch(
                                liquid_name,
                                active_database.database,
                                condition,
                                delta_T=delta_t_value,
                                unit=units,
                                phases=phases,
                                start_from_liquidus=start_from_liquidus,
                                n_cpu=batch_cpu,
                                progress_callback=progress_callback,
                            )

                        if solidification_model == "equilibrium_cooling":
                            result = run_batch(1)
                        else:
                            result = _run_parallel_batch_with_serial_retry(
                                run_batch,
                                n_cpu,
                                status,
                                progress_callback,
                                calculation_name="solidification",
                            )
                    elif solidification_model == "equilibrium_cooling":
                        result = eq.equilib_cooling(
                            liquid_name,
                            active_database.database,
                            condition,
                            delta_T=delta_t_value,
                            unit=units,
                            phases=phases,
                            start_from_liquidus=start_from_liquidus,
                            progress_callback=progress_callback,
                        )
                    elif solidification_model == "nucleoscheil":
                        result = eq.nucleoscheil_cooling(
                            liquid_name,
                            active_database.database,
                            condition,
                            critical_undercooling or None,
                            delta_T=delta_t_value,
                            unit=units,
                            phases=phases,
                            progress_callback=progress_callback,
                        )
                    else:
                        result = eq.scheil_cooling(
                            liquid_name,
                            active_database.database,
                            condition,
                            delta_T=delta_t_value,
                            unit=units,
                            phases=phases,
                            start_from_liquidus=start_from_liquidus,
                            progress_callback=progress_callback,
                        )
                    progress_callback(1, 1, "Finalizing solidification result...")
                    return result, None

                def apply_solver_result(payload: tuple[Any, str | None]) -> bool:
                    result, _error_message = payload
                    if module is not None:
                        module.result = result
                        module.runtime["result"] = result
                        _populate_module_results(module, results_table, result)
                        module.result_payload = _result_payload(result)
                        module.results_table_payload = _table_payload(results_table)
                    else:
                        _populate_scheil_results(results_table, result)
                    self._log_calculation_finish(module=module, result=result)
                    self._append_result_warnings(status, result)
                    self._set_calculation_finished_status(status)
                    return True

                if show_dialog:
                    return self._run_with_progress_dialog(
                        execute_solver,
                        apply_solver_result,
                        lambda exc: self._set_calculation_failure(
                            status,
                            f"Calculation failed: {exc}",
                            title="Calculation Error",
                            critical=True,
                            show_dialog=show_dialog,
                        ),
                        title="Calculating",
                        message="Running solidification calculation...",
                        canceled_callback=lambda: self._set_calculation_canceled_status(status),
                        lane_owned=lane_owned,
                    )
                if progress_callback is not None:
                    return apply_solver_result(execute_solver(progress_callback))
                if (
                    has_imported_solidification_conditions
                    or (is_solidification_module and has_manual_batch_conditions)
                ):
                    n_cpu = _effective_batch_cpu_count(
                        condition,
                        module.batch_cpu_count,
                    )
                    def run_batch(batch_cpu: int) -> Any:
                        if solidification_model == "equilibrium_cooling":
                            return _equilib_cooling_batch(
                                eq,
                                liquid_name,
                                active_database.database,
                                condition,
                                delta_T=float(delta_t.text()),
                                unit=units,
                                phases=phases,
                                start_from_liquidus=start_from_liquidus,
                            )
                        if solidification_model == "nucleoscheil":
                            return eq.nucleoscheil_batch(
                                liquid_name,
                                active_database.database,
                                condition,
                                critical_undercooling or None,
                                delta_T=float(delta_t.text()),
                                unit=units,
                                n_cpu=batch_cpu,
                                phases=phases,
                            )
                        return eq.scheil_batch(
                            liquid_name,
                            active_database.database,
                            condition,
                            delta_T=float(delta_t.text()),
                            unit=units,
                            phases=phases,
                            start_from_liquidus=start_from_liquidus,
                            n_cpu=batch_cpu,
                        )

                    if solidification_model == "equilibrium_cooling":
                        result = run_batch(1)
                    else:
                        result = _run_parallel_batch_with_serial_retry(
                            run_batch,
                            n_cpu,
                            status,
                            None,
                            calculation_name="solidification",
                        )
                elif solidification_model == "equilibrium_cooling":
                    result = eq.equilib_cooling(
                        liquid_name,
                        active_database.database,
                        condition,
                        delta_T=float(delta_t.text()),
                        unit=units,
                        phases=phases,
                        start_from_liquidus=start_from_liquidus,
                    )
                elif solidification_model == "nucleoscheil":
                    result = eq.nucleoscheil_cooling(
                        liquid_name,
                        active_database.database,
                        condition,
                        critical_undercooling or None,
                        delta_T=float(delta_t.text()),
                        unit=units,
                        phases=phases,
                    )
                else:
                    result = eq.scheil_cooling(
                        liquid_name,
                        active_database.database,
                        condition,
                        delta_T=float(delta_t.text()),
                        unit=units,
                        phases=phases,
                        start_from_liquidus=start_from_liquidus,
                    )
                if module is not None:
                    module.result = result
                    module.runtime["result"] = result
                    _populate_module_results(module, results_table, result)
                else:
                    _populate_scheil_results(results_table, result)
                if module is not None:
                    module.result_payload = _result_payload(result)
                    module.results_table_payload = _table_payload(results_table)
                self._log_calculation_finish(module=module, result=result)
            else:
                if has_imported_equilibrium_conditions:
                    equilibrium_condition = module.batch_condition
                elif has_manual_batch_conditions and condition is not None:
                    equilibrium_condition = condition
                else:
                    transition_widget = module.runtime.get("transition_search")
                    transition_search = (
                        isinstance(transition_widget, QCheckBox)
                        and transition_widget.isVisible()
                        and transition_widget.isChecked()
                    )
                    equilibrium_condition = _manual_equilibrium_condition(
                        active_database.database,
                        composition_table,
                        temperature.text(),
                        pressure.text(),
                        units,
                        phases=phases,
                        include_transitions=transition_search,
                    )
                n_cpu = _effective_batch_cpu_count(
                    equilibrium_condition,
                    module.batch_cpu_count,
                )
                def execute_solver(progress_callback) -> tuple[Any, str | None]:
                    progress_callback(0, 0, "Running equilibrium calculation...")
                    result = _run_equilibrium_batch_with_serial_retry(
                        eq,
                        active_database.database,
                        equilibrium_condition,
                        units,
                        phases,
                        n_cpu,
                        None,
                        progress_callback,
                    )
                    error_message = _calculation_error_message(result)
                    if error_message is not None and n_cpu > 1:
                        retry_result = eq.equilib_batch(
                            active_database.database,
                            equilibrium_condition,
                            unit=units,
                            phases=phases,
                            n_cpu=1,
                            progress_callback=progress_callback,
                        )
                        retry_error = _calculation_error_message(retry_result)
                        if retry_error is None:
                            result = retry_result
                            error_message = None
                    progress_callback(1, 1, "Finalizing equilibrium result...")
                    return result, error_message

                def apply_solver_result(payload: tuple[Any, str | None]) -> bool:
                    result, error_message = payload
                    if module is not None:
                        module.result = result
                        module.runtime["result"] = result
                        _populate_module_results(module, results_table, result)
                        module.result_payload = _result_payload(result)
                        module.results_table_payload = _table_payload(results_table)
                    else:
                        _populate_equilibrium_results(results_table, result)
                    if error_message is not None:
                        self._set_calculation_failure(
                            status,
                            error_message,
                            show_dialog=show_dialog,
                        )
                        return False
                    self._log_calculation_finish(module=module, result=result)
                    self._set_calculation_finished_status(status)
                    return True

                if show_dialog:
                    return self._run_with_progress_dialog(
                        execute_solver,
                        apply_solver_result,
                        lambda exc: self._set_calculation_failure(
                            status,
                            f"Calculation failed: {exc}",
                            title="Calculation Error",
                            critical=True,
                            show_dialog=show_dialog,
                        ),
                        title="Calculating",
                        message="Running equilibrium calculation...",
                        canceled_callback=lambda: self._set_calculation_canceled_status(status),
                        lane_owned=lane_owned,
                    )
                if progress_callback is not None:
                    return apply_solver_result(execute_solver(progress_callback))
                result = _run_equilibrium_batch_with_serial_retry(
                    eq,
                    active_database.database,
                    equilibrium_condition,
                    units,
                    phases,
                    n_cpu,
                    status,
                )
                error_message = _calculation_error_message(result)
                if error_message is not None and n_cpu > 1:
                    retry_result = eq.equilib_batch(
                        active_database.database,
                        equilibrium_condition,
                        unit=units,
                        phases=phases,
                        n_cpu=1,
                    )
                    retry_error = _calculation_error_message(retry_result)
                    if retry_error is None:
                        result = retry_result
                        error_message = None
                if module is not None:
                    module.result = result
                    module.runtime["result"] = result
                    _populate_module_results(module, results_table, result)
                else:
                    _populate_equilibrium_results(results_table, result)
                if module is not None:
                    module.result_payload = _result_payload(result)
                    module.results_table_payload = _table_payload(results_table)
                if error_message is not None:
                    self._set_calculation_failure(
                        status,
                        error_message,
                        show_dialog=show_dialog,
                    )
                    return False
                self._log_calculation_finish(module=module, result=result)
        except CalculationCanceled:
            self._set_calculation_canceled_status(status)
            if show_dialog:
                QMessageBox.information(
                    self,
                    "Calculation Canceled",
                    "Calculation canceled.",
                )
                return False
            raise
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            self._set_calculation_failure(
                status,
                f"Calculation failed: {exc}",
                title="Calculation Error",
                critical=True,
                show_dialog=show_dialog,
            )
            return False

        self._set_calculation_finished_status(status)
        return True

    def _set_calculation_failure(
        self,
        status: QPlainTextEdit,
        message: str,
        *,
        title: str = "Calculation Warning",
        critical: bool = False,
        show_dialog: bool = True,
    ) -> None:
        self._append_calculation_log_message(message)
        status.setPlainText(message)
        self._show_calculation_failure_dialog(
            title,
            message,
            critical=critical,
            show_dialog=show_dialog,
        )

    def _show_calculation_failure_dialog(
        self,
        title: str,
        message: str,
        *,
        critical: bool = False,
        show_dialog: bool = True,
    ) -> None:
        if not show_dialog:
            return
        if critical:
            QMessageBox.critical(self, title, message)
        else:
            QMessageBox.warning(self, title, message)


def _calculation_condition_point_count(condition: dict[str, Any] | None) -> int:
    """Return the number of condition points represented by a condition dict."""
    if not condition:
        return 1
    point_count = 1
    for value in condition.values():
        if isinstance(value, str):
            continue
        try:
            length = len(value)
        except TypeError:
            continue
        point_count = max(point_count, int(length))
    return point_count


def _initial_scheil_solid_phase_names(
    eq_module: Any,
    database: Any,
    condition: dict[str, Any],
    units: list[str],
    phases: list[str] | None,
    liquid_phase_name: str,
) -> list[str]:
    """Return solid phases stable at the first Scheil starting condition."""
    try:
        result = eq_module.equilib_single(
            database,
            _first_condition_point(condition),
            unit=units,
            phases=phases,
            include_heat_capacity=False,
        )
    except Exception:
        return []

    stable_phases = getattr(result, "stable_phases", None)
    stable_names = list(getattr(stable_phases, "names", []) or [])
    liquid_key = _phase_name_key(liquid_phase_name)
    return [
        phase_name
        for phase_name in stable_names
        if _phase_name_key(phase_name) != liquid_key
    ]


def _first_condition_point(condition: dict[str, Any]) -> dict[str, Any]:
    """Return scalar values from the first row of a condition dictionary."""
    first_point: dict[str, Any] = {}
    for key, value in condition.items():
        if isinstance(value, str):
            first_point[key] = value
            continue
        try:
            if len(value) > 0:
                first_point[key] = value[0]
                continue
        except (IndexError, KeyError, TypeError):
            pass
        first_point[key] = value
    return first_point


def _phase_name_key(phase_name: str) -> str:
    """Return a comparison key that ignores instance suffixes."""
    return str(phase_name).split("#", 1)[0].strip().upper()


def _default_phase_selection_for_module(
    module: CalculationModule | None,
    phase_names: list[str],
    session_default_phase_names: list[str] | None = None,
    database: dict[str, Any] | None = None,
) -> list[str]:
    """Return the GUI's default phase selection for one calculation module."""
    if session_default_phase_names:
        return [
            phase
            for phase in session_default_phase_names
            if phase in phase_names and not _is_reference_ser_phase(phase)
        ]
    if (
        module is not None
        and module.kind == "solidification"
        and _solidification_model_value(module.solidification_model) == "scheil"
    ):
        return _default_scheil_phase_selection(phase_names, database)
    return [
        phase
        for phase in phase_names
        if not _is_reference_ser_phase(phase)
    ]


def _default_scheil_excluded_ordered_phases_for_module(
    module: CalculationModule | None,
    phase_names: list[str],
    session_default_phase_names: list[str] | None,
    database: dict[str, Any] | None,
) -> list[str]:
    """Return ordered phases excluded only by the Scheil default selection."""
    if session_default_phase_names:
        return []
    if (
        module is None
        or module.kind != "solidification"
        or _solidification_model_value(module.solidification_model) != "scheil"
    ):
        return []
    return _scheil_default_excluded_ordered_phases(phase_names, database)


def _calculation_phase_selection_summary(
    phases: list[str] | None,
    phase_names: list[str],
) -> str:
    """Return a compact phase-selection summary for the calculation log."""
    if phases:
        return f"{len(phases)} selected phase(s)"
    if phase_names:
        return f"all {len(phase_names)} listed phase(s)"
    return "automatic phase selection"


def _calculation_result_summary(result: Any) -> str:
    """Return a compact result summary for the calculation log."""
    stable_phases = getattr(result, "stable_phases", None)
    stable_names = getattr(stable_phases, "names", None)
    if stable_names is not None:
        names = [str(name) for name in stable_names if str(name)]
        if names:
            return f"Stable phases: {', '.join(names)}."
        return "No stable phases reported."
    if isinstance(result, ResultTable):
        try:
            row_count = len(result.to_rows())
        except Exception:
            row_count = 0
        return f"Result table rows: {row_count}."
    result_table = getattr(result, "table", None)
    if isinstance(result_table, ResultTable):
        try:
            row_count = len(result_table.to_rows())
        except Exception:
            row_count = 0
        return f"Result table rows: {row_count}."
    return f"Result type: {type(result).__name__}."
