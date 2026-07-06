"""Load and save calculation sessions and module files."""

from __future__ import annotations

# ruff: noqa: F401,F403,F405,I001,E501

import copy
import csv
import json
import os
import pickle
import re
import shutil
from pathlib import Path
from typing import Any

from PySide6.QtCore import QModelIndex, QRect, QSize, Qt, QThread, QTimer
from PySide6.QtGui import QAction, QIcon, QPainter, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
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
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)
try:
    from PySide6.QtSvg import QSvgGenerator
except ImportError:  # pragma: no cover - depends on optional QtSvg packaging
    QSvgGenerator = None

from equilipy.composition import expand_condition_species as _expand_condition_species
from equilipy.database_ir import (
    Diagnostic,
    read_tdb as _read_tdb_ir,
    remove_redundant_disordered_phase_aliases,
    to_chemsage_compounds,
)
from equilipy.results import ResultTable
from equilipy.results.equilib import EquilibResult
from equilipy.results.scheil import ScheilResult

from equilipy.gui.assets import calculation_icon_path
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
_RESULT_FIGURE_EXPORT_DPI = 96.0


class CalculationSessionFilesMixin:
    """Calculation save, load, export, and database attachment methods."""

    def _save_current_session(self) -> None:
        session, _module = self._current_calculation_context()
        if session is None:
            QMessageBox.information(
                self,
                "Save Session",
                "Create or select a calculation session first.",
            )
            return

        try:
            if self.active_calculation_directory is not None:
                saved_count = self._sync_active_calculation_directory_files()
                QMessageBox.information(
                    self,
                    "Save Session",
                    f"Synchronized {saved_count} saved session/module file(s).",
                )
                return
            path = self._save_target_path(session.name, "Save Session", "session")
            if path is None:
                return
            payload = _self_contained_eq_payload(
                self._session_eq_payload(session),
                path,
            )
            eq_path = _write_eq_payload(path, payload)
            self._remove_legacy_script_path(eq_path.with_suffix(".py"))
            session.source_path = _eq_source_path(eq_path)
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Save Session Failed", str(exc))
            return

        QMessageBox.information(
            self,
            "Save Session",
            f"Saved {eq_path.name}.",
        )

    def _save_current_module(self) -> None:
        session, module = self._current_calculation_context()
        if session is None or module is None:
            QMessageBox.information(
                self,
                "Save Module",
                "Select a calculation module first.",
            )
            return

        try:
            self._ensure_calculation_module_runtime(session, module)
            if self.active_calculation_directory is not None:
                saved_count = self._sync_active_calculation_directory_files()
                QMessageBox.information(
                    self,
                    "Save Module",
                    f"Synchronized {saved_count} saved session/module file(s).",
                )
                return
            path = self._save_target_path(module.name, "Save Module", "module")
            if path is None:
                return
            payload = _self_contained_eq_payload(
                self._module_eq_payload(session, module),
                path,
            )
            module_eq_path = _eq_path_for_save(path)
            script_path = self._module_script_path(module_eq_path, session)
            eq_path, script_path = _save_eq_and_script(
                path,
                payload,
                _module_script_text(payload["modules"][0], payload["session"]),
                script_path=script_path,
            )
            self._remove_legacy_script_path(eq_path.with_suffix(".py"), script_path)
            module.source_path = _eq_source_path(eq_path)
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Save Module Failed", str(exc))
            return

        QMessageBox.information(
            self,
            "Save Module",
            f"Saved {eq_path.name} and {script_path.name}.",
        )

    def _sync_active_calculation_directory_files(self) -> int:
        """Write current sessions/modules and remove stale active-directory files."""
        if self.active_calculation_directory is None:
            return 0

        directory = self.active_calculation_directory.expanduser().resolve()
        directory.mkdir(parents=True, exist_ok=True)
        self._remove_stale_calculation_source_paths(directory)

        used_paths: set[Path] = set()
        saved_count = 0

        for session in self.calculation_sessions.values():
            session_path = self._active_directory_save_path(
                directory,
                session.name,
                used_paths,
            )
            payload = _self_contained_eq_payload(
                self._session_eq_payload(session),
                session_path,
            )
            eq_path = _write_eq_payload(str(session_path), payload)
            self._remove_legacy_script_path(eq_path.with_suffix(".py"))
            session.source_path = _eq_source_path(eq_path)
            used_paths.add(eq_path.resolve())
            saved_count += 1

            for module in session.modules.values():
                module_path = self._active_directory_save_path(
                    directory,
                    module.name,
                    used_paths,
                    fallback_prefix=session.name,
                )
                payload = _self_contained_eq_payload(
                    self._module_eq_payload(session, module),
                    module_path,
                )
                script_path = self._module_script_path(module_path, session)
                eq_path, script_path = _save_eq_and_script(
                    str(module_path),
                    payload,
                    _module_script_text(payload["modules"][0], payload["session"]),
                    script_path=script_path,
                )
                self._remove_legacy_script_path(eq_path.with_suffix(".py"), script_path)
                module.source_path = _eq_source_path(eq_path)
                used_paths.add(eq_path.resolve())
                saved_count += 1

        return saved_count

    def _session_script_prefix(self, session: CalculationSession) -> str:
        match = re.search(r"(\d+)$", str(session.id))
        if match is not None:
            return f"Session{int(match.group(1))}"
        return _safe_file_stem(session.name).replace("_", "") or "Session"

    def _module_script_path(
        self,
        module_eq_path: str | Path,
        session: CalculationSession,
    ) -> Path:
        eq_path = _eq_path_for_save(module_eq_path)
        script_name = f"{self._session_script_prefix(session)}_{eq_path.stem}.py"
        return eq_path.with_name(script_name)

    def _remove_legacy_script_path(
        self,
        path: str | Path,
        keep_path: str | Path | None = None,
    ) -> None:
        legacy_path = Path(path).expanduser()
        keep_resolved = None
        if keep_path is not None:
            try:
                keep_resolved = Path(keep_path).expanduser().resolve()
            except OSError:
                keep_resolved = Path(keep_path).expanduser().absolute()
        try:
            legacy_resolved = legacy_path.resolve()
        except OSError:
            legacy_resolved = legacy_path.absolute()
        if keep_resolved is not None and legacy_resolved == keep_resolved:
            return
        try:
            if legacy_path.is_file():
                legacy_path.unlink()
        except OSError:
            return

    def _active_directory_save_path(
        self,
        directory: Path,
        name: str,
        used_paths: set[Path],
        *,
        fallback_prefix: str | None = None,
    ) -> Path:
        stem = _safe_file_stem(name)
        candidate = (directory / f"{stem}.eq").resolve()
        if candidate not in used_paths:
            return candidate

        if fallback_prefix:
            stem = f"{_safe_file_stem(fallback_prefix)}_{stem}"
            candidate = (directory / f"{stem}.eq").resolve()
            if candidate not in used_paths:
                return candidate

        index = 2
        while True:
            indexed = (directory / f"{stem}_{index}.eq").resolve()
            if indexed not in used_paths:
                return indexed
            index += 1

    def _remember_stale_calculation_source_path(
        self,
        path: str | Path | None,
    ) -> None:
        if not path or self.active_calculation_directory is None:
            return
        source_path = Path(path).expanduser()
        if source_path.suffix.lower() not in {".eq", ".py"}:
            return
        try:
            source_path = source_path.resolve()
            active_directory = self.active_calculation_directory.expanduser().resolve()
            source_path.relative_to(active_directory)
        except (OSError, ValueError):
            return
        stale_paths = getattr(self, "_stale_calculation_source_paths", None)
        if stale_paths is None:
            self._stale_calculation_source_paths = set()
            stale_paths = self._stale_calculation_source_paths
        stale_paths.add(str(source_path))

    def _remember_stale_module_source_path(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> None:
        if module.source_path and module.source_path != session.source_path:
            self._remember_stale_calculation_source_path(module.source_path)
            self._remember_stale_calculation_source_path(
                self._module_script_path(module.source_path, session)
            )

    def _remove_stale_calculation_source_paths(self, directory: Path) -> None:
        stale_paths = getattr(self, "_stale_calculation_source_paths", set())
        if not stale_paths:
            return

        active_directory = directory.resolve()
        for source_text in sorted(stale_paths):
            source_path = Path(source_text).expanduser()
            try:
                source_path = source_path.resolve()
                source_path.relative_to(active_directory)
            except (OSError, ValueError):
                continue
            for path in (source_path, source_path.with_suffix(".py")):
                try:
                    if path.is_file():
                        path.unlink()
                except OSError:
                    continue
        stale_paths.clear()

    def _save_target_path(
        self,
        name: str,
        title: str,
        scope: str,
    ) -> str | None:
        default_name = f"{_safe_file_stem(name)}.eq"
        if self.active_calculation_directory is not None:
            self.active_calculation_directory.mkdir(parents=True, exist_ok=True)
            return str(self.active_calculation_directory / default_name)

        selected_path, _selected_filter = QFileDialog.getSaveFileName(
            self,
            title,
            default_name,
            f"Equilipy {scope} (*.eq);;All files (*)",
        )
        return selected_path or None

    def _export_module_results(
        self,
        module: CalculationModule,
        table: QTableWidget,
    ) -> None:
        if _module_result_figure_is_active(module):
            self._export_module_figure(module)
            return

        if table.columnCount() == 0 or table.rowCount() == 0:
            QMessageBox.information(
                self,
                "Export Results",
                "Run or load a calculation result before exporting.",
            )
            return

        session = self.calculation_sessions.get(
            str(module.runtime.get("session_id") or "")
        )
        session_name = session.name if session is not None else "Session"
        default_name = (
            f"{_safe_file_stem(session_name)}_"
            f"{_safe_file_stem(module.name)}_result.csv"
        )
        if self.active_calculation_directory is not None:
            self.active_calculation_directory.mkdir(parents=True, exist_ok=True)
            target_path = self.active_calculation_directory / default_name
        else:
            selected_path, _selected_filter = QFileDialog.getSaveFileName(
                self,
                "Export Results",
                default_name,
                "CSV files (*.csv);;All files (*)",
            )
            if not selected_path:
                return
            target_path = Path(selected_path).expanduser()
            if target_path.suffix.lower() != ".csv":
                target_path = target_path.with_suffix(".csv")

        try:
            _write_module_result_csv(module, table, target_path)
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Export Results Failed", str(exc))
            return

        QMessageBox.information(
            self,
            "Export Results",
            f"Saved {target_path.name}.",
        )

    def _export_module_figure(self, module: CalculationModule) -> None:
        chart_view = module.runtime.get("result_figure_chart")
        if not _result_figure_chart_has_series(chart_view):
            QMessageBox.information(
                self,
                "Export Figure",
                "Generate a figure before exporting.",
            )
            return
        if QSvgGenerator is None:
            QMessageBox.critical(
                self,
                "Export Figure Failed",
                "This Qt build does not include SVG export support.",
            )
            return

        export_size = _result_figure_export_size(self, chart_view)
        if export_size is None:
            return

        session = self.calculation_sessions.get(
            str(module.runtime.get("session_id") or "")
        )
        session_name = session.name if session is not None else "Session"
        default_name = (
            f"{_safe_file_stem(session_name)}_"
            f"{_safe_file_stem(module.name)}_figure.svg"
        )
        if self.active_calculation_directory is not None:
            self.active_calculation_directory.mkdir(parents=True, exist_ok=True)
            target_path = self.active_calculation_directory / default_name
        else:
            selected_path, _selected_filter = QFileDialog.getSaveFileName(
                self,
                "Export Figure",
                default_name,
                "SVG files (*.svg);;All files (*)",
            )
            if not selected_path:
                return
            target_path = Path(selected_path).expanduser()
            if target_path.suffix.lower() != ".svg":
                target_path = target_path.with_suffix(".svg")

        try:
            _write_result_figure_svg(chart_view, target_path, export_size)
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Export Figure Failed", str(exc))
            return

        QMessageBox.information(
            self,
            "Export Figure",
            f"Saved {target_path.name}.",
        )

    def _set_active_calculation_directory(self) -> None:
        start_dir = (
            str(self.active_calculation_directory)
            if self.active_calculation_directory is not None
            else ""
        )
        directory = QFileDialog.getExistingDirectory(
            self,
            "Set Active Directory",
            start_dir,
        )
        if not directory:
            return

        self.active_calculation_directory = Path(directory).expanduser().resolve()
        self._update_active_directory_button()
        self._update_log_capture_file()
        self._offer_active_directory_load()

    def _show_active_calculation_directory_menu(self, position) -> None:
        if self.calculation_active_directory_button is None:
            return

        menu = self._active_calculation_directory_menu()
        menu.exec(self.calculation_active_directory_button.mapToGlobal(position))

    def _active_calculation_directory_menu(self) -> QMenu:
        menu = QMenu(self)
        set_directory = QAction("Set Active Directory", self)
        set_directory.triggered.connect(
            lambda _checked=False: self._set_active_calculation_directory()
        )
        menu.addAction(set_directory)

        if self.active_calculation_directory is not None:
            load_directory = QAction("Load Current Directory", self)
            close_directory = QAction("Close Current Directory", self)
            load_directory.triggered.connect(
                lambda _checked=False: self._offer_active_directory_load()
            )
            close_directory.triggered.connect(
                lambda _checked=False: self._close_active_calculation_directory()
            )
            menu.addAction(load_directory)
            menu.addSeparator()
            menu.addAction(close_directory)

        return menu

    def _close_active_calculation_directory(self) -> None:
        self.active_calculation_directory = None
        self._update_active_directory_button()
        self._update_log_capture_file()

    def _update_active_directory_button(self) -> None:
        if self.calculation_active_directory_button is None:
            return
        if self.active_calculation_directory is None:
            self.calculation_active_directory_button.setToolTip(
                "Set active calculation directory"
            )
            return
        self.calculation_active_directory_button.setToolTip(
            f"Active directory: {self.active_calculation_directory}"
        )

    def _offer_active_directory_load(self) -> None:
        if self.active_calculation_directory is None:
            return
        payloads = _eq_payloads_in_directory(self.active_calculation_directory)
        payloads = self._new_eq_payloads(payloads)
        if not payloads:
            return

        if not self._active_directory_load_choice(len(payloads)):
            return

        allow_pickle_results = self._pickled_result_load_choice_for_payloads(payloads)
        if allow_pickle_results is None:
            return

        try:
            loaded_count = self._load_eq_payloads_as_sessions(
                payloads,
                allow_pickle_results=allow_pickle_results,
            )
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Load Active Directory Failed", str(exc))
            return
        QMessageBox.information(
            self,
            "Load Active Directory",
            f"Loaded {loaded_count} session(s)/module group(s).",
        )

    def _active_directory_load_choice(self, payload_count: int) -> bool:
        dialog = QDialog(self)
        dialog.setWindowTitle("Load Active Directory")
        dialog.setModal(True)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(14)

        message = QLabel(
            f"Found {payload_count} saved session/module file(s) in:\n"
            f"{self.active_calculation_directory}\n\nLoad them now?"
        )
        message.setObjectName("DatabaseListLabel")
        message.setWordWrap(True)
        layout.addWidget(message)

        actions = QHBoxLayout()
        actions.setSpacing(8)
        yes = QPushButton("Yes")
        yes.setObjectName("PrimaryButton")
        no = QPushButton("No")
        no.setObjectName("SmallActionButton")
        yes.setDefault(True)
        yes.setAutoDefault(True)
        no.setAutoDefault(False)
        yes.clicked.connect(dialog.accept)
        no.clicked.connect(dialog.reject)

        actions.addWidget(yes)
        actions.addWidget(no)
        layout.addLayout(actions)
        return dialog.exec() == QDialog.DialogCode.Accepted

    def _pickled_result_load_choice_for_payloads(
        self,
        payloads: list[tuple[Path, dict[str, Any]]],
    ) -> bool:
        return True

    def _load_eq_payloads_as_sessions(
        self,
        payloads: list[tuple[Path, dict[str, Any]]],
        *,
        allow_pickle_results: bool = False,
    ) -> int:
        payloads = self._new_eq_payloads(payloads)
        session_payloads = [
            (path, payload)
            for path, payload in payloads
            if payload.get("scope") == "session"
        ]
        module_payloads = [
            (path, payload)
            for path, payload in payloads
            if payload.get("scope") != "session"
        ]
        loaded_sessions_by_name: dict[str, CalculationSession] = {}
        loaded_modules_by_key: dict[tuple[str, str, str], CalculationModule] = {}
        loaded_count = 0

        for path, payload in session_payloads:
            session_name = str(
                (payload.get("session") or {}).get("name") or "Session"
            ).strip() or "Session"
            session = self._create_calculation_session(session_name)
            self._add_session_tree_item(session)
            self._load_eq_payload_into_session(
                session,
                payload,
                replace=True,
                allow_pickle_results=allow_pickle_results,
                source_path=path,
            )
            loaded_sessions_by_name[session_name] = session
            loaded_sessions_by_name[session.name] = session
            for module in session.modules.values():
                loaded_modules_by_key[
                    (session_name, module.kind, module.name)
                ] = module
                loaded_modules_by_key[(session.name, module.kind, module.name)] = module
            loaded_count += 1

        grouped_modules: dict[str, list[tuple[Path, dict[str, Any]]]] = {}
        grouped_session_payloads: dict[str, dict[str, Any]] = {}
        grouped_session_sources: dict[str, Path] = {}
        for path, payload in module_payloads:
            session_payload = payload.get("session") or {}
            session_name = (
                str(session_payload.get("name") or "Session").strip() or "Session"
            )
            grouped_session_payloads.setdefault(session_name, session_payload)
            grouped_session_sources.setdefault(session_name, path)
            grouped_modules.setdefault(session_name, []).extend(
                (path, module_state) for module_state in payload.get("modules") or []
            )

        for session_name, module_entries in grouped_modules.items():
            session = loaded_sessions_by_name.get(session_name)
            if session is None:
                session = self._create_calculation_session(session_name)
                self._add_session_tree_item(session)
                self._load_session_state_only(
                    session,
                    grouped_session_payloads.get(session_name, {}),
                    source_path=grouped_session_sources.get(session_name),
                )
                loaded_sessions_by_name[session_name] = session
                loaded_sessions_by_name[session.name] = session
                loaded_count += 1

            for path, state in module_entries:
                kind = str(state.get("kind") or "equilibrium")
                name = str(state.get("name") or self._module_label(kind))
                saved_module_key = (session_name, kind, name)
                actual_module_key = (session.name, kind, name)
                existing_module = loaded_modules_by_key.get(
                    saved_module_key
                ) or loaded_modules_by_key.get(actual_module_key)
                if existing_module is not None:
                    existing_module.source_path = _eq_source_path(path)
                    current_state = {
                        "result": existing_module.result_payload,
                        "results_table": existing_module.results_table_payload,
                    }
                    if _module_state_result_score(state) > _module_state_result_score(
                        current_state
                    ):
                        self._apply_saved_module_state(
                            session,
                            existing_module,
                            state,
                            allow_pickle_results=allow_pickle_results,
                            source_path=path,
                        )
                        loaded_modules_by_key[saved_module_key] = existing_module
                        loaded_modules_by_key[
                            (
                                session.name,
                                existing_module.kind,
                                existing_module.name,
                            )
                        ] = existing_module
                    continue
                module = self._add_calculation_module(session, kind)
                self._apply_saved_module_state(
                    session,
                    module,
                    state,
                    allow_pickle_results=allow_pickle_results,
                    source_path=path,
                )
                loaded_modules_by_key[saved_module_key] = module
                loaded_modules_by_key[
                    (session.name, module.kind, module.name)
                ] = module

        active = self._active_calculation_session()
        if active is not None:
            self._show_calculation_overview(active)
        return loaded_count

    def _loaded_eq_source_paths(self) -> set[str]:
        paths: set[str] = set()
        for session in self.calculation_sessions.values():
            if session.source_path:
                paths.add(session.source_path)
            for module in session.modules.values():
                if module.source_path:
                    paths.add(module.source_path)
        return paths

    def _new_eq_payloads(
        self,
        payloads: list[tuple[Path, dict[str, Any]]],
    ) -> list[tuple[Path, dict[str, Any]]]:
        loaded_paths = self._loaded_eq_source_paths()
        seen_paths: set[str] = set()
        new_payloads: list[tuple[Path, dict[str, Any]]] = []
        for path, payload in payloads:
            source_path = _eq_source_path(path)
            if source_path in loaded_paths or source_path in seen_paths:
                continue
            seen_paths.add(source_path)
            new_payloads.append((path, payload))
        return new_payloads

    def _load_session_state_only(
        self,
        session: CalculationSession,
        session_payload: dict[str, Any],
        *,
        source_path: str | Path | None = None,
    ) -> None:
        saved_name = str(session_payload.get("name") or "").strip()
        if saved_name:
            session.name = self._unique_calculation_session_name(
                saved_name,
                exclude_session_id=session.id,
            )
        session.temperature_unit = str(
            session_payload.get("temperature_unit") or session.temperature_unit
        )
        session.pressure_unit = str(
            session_payload.get("pressure_unit") or session.pressure_unit
        )
        session.amount_unit = str(
            session_payload.get("amount_unit") or session.amount_unit
        )
        session.default_phase_names = list(
            session_payload.get("default_phase_names") or []
        )
        session.result_column_defaults = {
            str(kind): [str(column) for column in (columns or [])]
            for kind, columns in (
                session_payload.get("result_column_defaults") or {}
            ).items()
        }
        self._load_session_databases(
            session,
            session_payload.get("databases", []),
            source_path=source_path,
        )
        item = self._find_session_item(session.id)
        if item is not None:
            self._set_calculation_tree_text(item, session.name)

    def _load_session_from_file(self) -> None:
        session, _module = self._current_calculation_context()
        if session is None:
            session = self._active_calculation_session(create_if_missing=True)
        if session is None:
            return

        path, _selected_filter = QFileDialog.getOpenFileName(
            self,
            "Load Session",
            "",
            "Equilipy session/module (*.eq);;All files (*)",
        )
        if not path:
            return

        try:
            payload = _read_eq_payload(path)
            allow_pickle_results = self._pickled_result_load_choice(payload)
            if allow_pickle_results is None:
                return
            self._load_eq_payload_into_session(
                session,
                payload,
                replace=True,
                allow_pickle_results=allow_pickle_results,
                source_path=path,
            )
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Load Session Failed", str(exc))
            return

        self._show_calculation_overview(session)
        self._update_current_calculation_heading(session)

    def _load_module_from_file(self) -> None:
        session, module = self._current_calculation_context()
        if session is None or module is None:
            QMessageBox.information(
                self,
                "Load Module",
                "Select a calculation module first.",
            )
            return

        path, _selected_filter = QFileDialog.getOpenFileName(
            self,
            "Load Module",
            "",
            "Equilipy module/session (*.eq);;All files (*)",
        )
        if not path:
            return

        try:
            payload = _read_eq_payload(path)
            allow_pickle_results = self._pickled_result_load_choice(payload)
            if allow_pickle_results is None:
                return
            modules = payload.get("modules") or []
            if not modules:
                raise ValueError("No modules were found in the selected .eq file.")
            session_payload = payload.get("session", {})
            self._load_session_databases(
                session,
                session_payload.get("databases", []),
                source_path=path,
            )
            state = modules[0]
            if state.get("kind") == module.kind:
                loaded_module = module
            else:
                loaded_module = self._add_calculation_module(
                    session,
                    str(state.get("kind") or "equilibrium"),
                )
            self._apply_saved_module_state(
                session,
                loaded_module,
                state,
                allow_pickle_results=allow_pickle_results,
                source_path=path,
            )
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Load Module Failed", str(exc))
            return

        self._show_calculation_module(session, loaded_module)
        self._revalidate_session_composition_tables(session)

    def _session_eq_payload(self, session: CalculationSession) -> dict[str, Any]:
        modules = []
        for module in session.modules.values():
            self._ensure_calculation_module_runtime(session, module)
            modules.append(_module_state_payload(session, module))
        return {
            "format": "equilipy.eq",
            "version": 1,
            "scope": "session",
            "session": _session_state_payload(session),
            "modules": modules,
        }

    def _module_eq_payload(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> dict[str, Any]:
        return {
            "format": "equilipy.eq",
            "version": 1,
            "scope": "module",
            "session": _session_state_payload(session),
            "modules": [_module_state_payload(session, module)],
        }

    def _load_eq_payload_into_session(
        self,
        session: CalculationSession,
        payload: dict[str, Any],
        *,
        replace: bool,
        allow_pickle_results: bool = False,
        source_path: str | Path | None = None,
    ) -> None:
        session_payload = payload.get("session", {})
        modules = payload.get("modules") or []
        if replace:
            self._clear_session_modules(session)
        if source_path is not None:
            session.source_path = _eq_source_path(source_path)

        saved_name = str(session_payload.get("name") or "").strip()
        if saved_name:
            session.name = self._unique_calculation_session_name(
                saved_name,
                exclude_session_id=session.id,
            )
        session.temperature_unit = str(
            session_payload.get("temperature_unit") or session.temperature_unit
        )
        session.pressure_unit = str(
            session_payload.get("pressure_unit") or session.pressure_unit
        )
        session.amount_unit = str(
            session_payload.get("amount_unit") or session.amount_unit
        )
        session.default_phase_names = list(
            session_payload.get("default_phase_names") or []
        )
        session.result_column_defaults = {
            str(kind): [str(column) for column in (columns or [])]
            for kind, columns in (
                session_payload.get("result_column_defaults") or {}
            ).items()
        }
        self._load_session_databases(
            session,
            session_payload.get("databases", []),
            source_path=source_path,
        )

        session_item = self._find_session_item(session.id)
        if session_item is not None:
            self._set_calculation_tree_text(session_item, session.name)

        for state in modules:
            kind = str(state.get("kind") or "equilibrium")
            module = self._add_calculation_module(session, kind)
            self._apply_saved_module_state(
                session,
                module,
                state,
                allow_pickle_results=allow_pickle_results,
                source_path=source_path,
            )

    def _clear_session_modules(self, session: CalculationSession) -> None:
        for module_id in list(session.modules):
            self._remove_calculation_page(self._module_key(module_id))
        session.modules.clear()
        session.module_counts.clear()

        session_item = self._find_session_item(session.id)
        if session_item is not None:
            while session_item.childCount():
                session_item.removeChild(session_item.child(0))

    def _load_session_databases(
        self,
        session: CalculationSession,
        database_entries: list[dict[str, Any]],
        *,
        source_path: str | Path | None = None,
    ) -> None:
        session.databases.clear()
        session.database_counter = 0
        failures: list[str] = []
        base_dir = (
            Path(source_path).expanduser().resolve().parent
            if source_path
            else None
        )

        for entry in database_entries:
            path = str(entry.get("path") or "").strip()
            if not path:
                continue
            database_path_obj = _resolve_saved_database_path(path, base_dir)
            database_path = str(database_path_obj)
            try:
                database = self._load_runtime_database_for_calculation(
                    database_path_obj
                )
                if database is None:
                    continue
            except Exception as exc:  # pragma: no cover - exercised manually in GUI.
                failures.append(f"{Path(path).name}: {exc}")
                continue

            session.database_counter += 1
            database_id = f"{session.id}:database:{session.database_counter}"
            session.databases[database_id] = CalculationDatabase(
                id=database_id,
                name=str(entry.get("name") or Path(path).name),
                path=database_path,
                database=database,
                selected=bool(entry.get("selected", False)),
            )

        if session.databases and not any(
            db.selected for db in session.databases.values()
        ):
            next(iter(session.databases.values())).selected = True
        if failures:
            QMessageBox.warning(
                self,
                "Database Load Warning",
                "Some saved databases could not be loaded:\n" + "\n".join(failures),
            )

    def _apply_saved_module_state(
        self,
        session: CalculationSession,
        module: CalculationModule,
        state: dict[str, Any],
        *,
        allow_pickle_results: bool = False,
        source_path: str | Path | None = None,
    ) -> None:
        module.kind = str(state.get("kind") or module.kind)
        module.name = self._unique_calculation_module_name(
            session,
            str(state.get("name") or module.name),
            exclude_module_id=module.id,
        )
        if source_path is not None:
            module.source_path = _eq_source_path(source_path)
        module.temperature_unit = str(state.get("temperature_unit") or "K")
        module.pressure_unit = str(state.get("pressure_unit") or "atm")
        module.amount_unit = str(state.get("amount_unit") or "moles")
        module.calculation_type = str(state.get("calculation_type") or "single")
        module.solidification_model = str(
            state.get("solidification_model") or "scheil"
        )
        module.nucleation_undercooling = {
            str(phase): float(value)
            for phase, value in (
                state.get("nucleation_undercooling") or {}
            ).items()
        }
        module.batch_condition = {
            key: [float(value) for value in values]
            for key, values in (state.get("batch_condition") or {}).items()
        }
        module.batch_condition_path = str(state.get("batch_condition_path") or "")
        module.batch_cpu_count = _normalized_batch_cpu_count(
            state.get("batch_cpu_count", 0)
        )
        module.transition_search = bool(state.get("transition_search", True))
        module.start_from_liquidus = bool(state.get("start_from_liquidus", True))
        module.result_columns = [
            str(column) for column in (state.get("result_columns") or [])
        ]
        module.phase_names = list(state.get("phase_names") or [])

        item = self._find_module_item(session.id, module.id)
        if item is not None:
            self._updating_calculation_tree = True
            item.setIcon(0, QIcon(str(calculation_icon_path(module.kind))))
            item.setText(0, module.name)
            self._updating_calculation_tree = False

        runtime = self._ensure_calculation_module_runtime(session, module)
        _apply_module_state_to_runtime(
            module,
            runtime,
            state,
            allow_pickle_results=allow_pickle_results,
        )
        self._apply_calculation_type_to_runtime(module)
        self._set_header_unit_controls(module if module.kind != "nucleation" else None)

    def _pickled_result_load_choice(
        self,
        payload: dict[str, Any],
    ) -> bool:
        return True

    def _load_session_database(
        self,
        session_id: str,
    ) -> None:
        path, _selected_filter = QFileDialog.getOpenFileName(
            self,
            "Open Equilipy Database",
            "",
            "Equilipy databases (*.dat *.tdb);;DAT databases (*.dat);;TDB databases (*.tdb);;All files (*)",
        )
        if not path:
            return
        session = self.calculation_sessions[session_id]
        source_database_path = Path(path).expanduser().resolve()
        database_name = Path(path).name
        duplicate = next(
            (
                loaded
                for loaded in session.databases.values()
                if Path(loaded.path).expanduser().resolve() == source_database_path
                or loaded.name == database_name
            ),
            None,
        )
        if duplicate is not None:
            QMessageBox.warning(
                self,
                "Database Already Loaded",
                f"{duplicate.name} is already loaded in this session.",
            )
            return
        try:
            database_path_obj = source_database_path
            if self.active_calculation_directory is not None:
                database_path_obj = _copy_database_for_directory(
                    source_database_path,
                    self.active_calculation_directory,
                )
            database_path = str(database_path_obj)
            database = self._load_runtime_database_for_calculation(database_path_obj)
            if database is None:
                return
        except Exception as exc:  # pragma: no cover - exercised manually in GUI.
            QMessageBox.critical(self, "Database Load Failed", str(exc))
            return

        session.database_counter += 1
        database_id = f"{session.id}:database:{session.database_counter}"
        for loaded in session.databases.values():
            loaded.selected = False
        session.databases[database_id] = CalculationDatabase(
            id=database_id,
            name=database_name,
            path=database_path,
            database=database,
            selected=True,
        )
        self._show_calculation_overview(session)
        self._update_current_calculation_heading(session)
        self._revalidate_session_composition_tables(session)

    def _load_runtime_database_for_calculation(self, path: str | Path) -> Any | None:
        """Load a calculation runtime database from a DAT or TDB file."""
        import equilipy as eq

        source_path = Path(path).expanduser().resolve()
        suffix = source_path.suffix.lower()
        if suffix == ".dat":
            return eq.read_dat(str(source_path))
        if suffix == ".tdb":
            return self._load_tdb_runtime_database_for_calculation(source_path)
        raise ValueError(f"Unsupported database file type: {source_path.name}")

    def _load_tdb_runtime_database_for_calculation(
        self,
        path: Path,
    ) -> Any | None:
        """Validate, optionally correct, and load a TDB runtime database."""
        progress = self._tdb_load_progress_dialog(path)
        try:
            progress.setLabelText(f"Reading {path.name}...")
            progress.setValue(0)
            QApplication.processEvents()
            preview = _read_tdb_ir(
                path,
                strict=False,
                remove_redundant_phases=False,
            )

            progress.setLabelText(f"Validating {path.name}...")
            progress.setValue(1)
            QApplication.processEvents()
            report = preview.validate_tdb(auto_correct=False)

            progress.setLabelText(f"Checking corrections for {path.name}...")
            progress.setValue(2)
            QApplication.processEvents()
            corrected_preview = copy.deepcopy(preview)
            alias_count = remove_redundant_disordered_phase_aliases(
                corrected_preview
            )
            corrected_report = corrected_preview.validate_tdb(auto_correct=True)
            correction_count = alias_count + corrected_report.corrected_count
        finally:
            progress.close()

        if correction_count:
            choice = self._confirm_tdb_validation_action(
                path,
                report.diagnostics,
                corrected_report.diagnostics,
                correction_count,
                allow_use_as_is=not report.errors,
            )
            if choice == "cancel":
                return None
            if choice == "apply":
                if corrected_report.errors:
                    self._show_tdb_blocking_diagnostics(
                        path,
                        corrected_report.errors,
                    )
                    return None
                progress = self._tdb_load_progress_dialog(path)
                progress.setLabelText(
                    f"Applying corrections and loading {path.name}..."
                )
                progress.setValue(3)
                QApplication.processEvents()
                try:
                    return to_chemsage_compounds(corrected_preview)
                finally:
                    progress.close()
            runtime_preview = copy.deepcopy(preview)
            remove_redundant_disordered_phase_aliases(runtime_preview)
        else:
            runtime_preview = preview

        if report.errors:
            self._show_tdb_blocking_diagnostics(path, report.errors)
            return None
        progress = self._tdb_load_progress_dialog(path)
        progress.setLabelText(f"Loading runtime database from {path.name}...")
        progress.setValue(3)
        QApplication.processEvents()
        try:
            return to_chemsage_compounds(runtime_preview)
        finally:
            progress.close()

    def _tdb_load_progress_dialog(self, path: Path) -> QProgressDialog:
        """Return a modal progress dialog for TDB validation/loading."""
        progress = QProgressDialog(
            f"Validating database {path.name}...",
            "",
            0,
            4,
            self,
        )
        progress.setWindowTitle("Loading TDB Database")
        progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.setMinimumDuration(0)
        progress.setCancelButton(None)
        progress.show()
        QApplication.processEvents()
        return progress

    def _confirm_tdb_validation_action(
        self,
        path: Path,
        diagnostics: list[Diagnostic],
        corrected_diagnostics: list[Diagnostic],
        correction_count: int,
        *,
        allow_use_as_is: bool,
    ) -> str:
        """Ask the user whether to apply available TDB import corrections."""
        dialog = _ScrollableMessageDialog(
            parent=self,
            title="TDB Validation",
            icon=QMessageBox.Icon.Warning,
            headline=f"{path.name} has correctable TDB validation diagnostics.",
            details=_tdb_validation_prompt_text(
                diagnostics,
                corrected_diagnostics,
                correction_count,
                allow_use_as_is=allow_use_as_is,
            ),
        )
        apply_button = dialog.add_action("Apply Corrections", accept=True)
        use_as_is_button = dialog.add_action("Use As-Is")
        use_as_is_button.setEnabled(allow_use_as_is)
        dialog.add_action("Cancel")
        dialog.set_default_button(apply_button)
        dialog.exec()
        clicked = dialog.clicked_button
        if clicked is apply_button:
            return "apply"
        if clicked is use_as_is_button and allow_use_as_is:
            return "use_as_is"
        return "cancel"

    def _show_tdb_blocking_diagnostics(
        self,
        path: Path,
        diagnostics: list[Diagnostic],
    ) -> None:
        """Display blocking TDB diagnostics before calculation database load."""
        _ScrollableMessageDialog.information(
            parent=self,
            title="TDB Validation Failed",
            icon=QMessageBox.Icon.Critical,
            headline=(
                f"{path.name} has blocking validation errors and was not loaded."
            ),
            details=_format_tdb_diagnostics(diagnostics),
        )

    def _set_session_database_selected(
        self,
        session_id: str,
        database_id: str,
        selected: bool,
    ) -> None:
        session = self.calculation_sessions.get(session_id)
        if session is None or database_id not in session.databases:
            return
        if selected:
            for loaded_id, database in session.databases.items():
                database.selected = loaded_id == database_id
        else:
            session.databases[database_id].selected = False
        self._refresh_visible_session_database_page(session)
        self._update_current_calculation_heading(session)
        self._revalidate_session_composition_tables(session)

    def _remove_session_database(
        self,
        session_id: str,
        database_id: str,
    ) -> None:
        session = self.calculation_sessions.get(session_id)
        if session is None:
            return
        session.databases.pop(database_id, None)
        self._refresh_visible_session_database_page(session)
        self._update_current_calculation_heading(session)
        self._revalidate_session_composition_tables(session)

    def _refresh_visible_session_database_page(
        self,
        session: CalculationSession,
    ) -> None:
        if self.calculation_stack is None:
            return
        current = self.calculation_stack.currentWidget()
        if (
            current is not None
            and current.property("session_id") == session.id
            and current.property("module_id") is None
        ):
            self._show_calculation_overview(session)

    def _active_session_database(
        self,
        session: CalculationSession,
    ) -> CalculationDatabase | None:
        for database in session.databases.values():
            if database.selected:
                return database
        return None


class _ScrollableMessageDialog(QDialog):
    """Bounded message dialog with a scrollable plain-text details area."""

    def __init__(
        self,
        *,
        parent: QWidget,
        title: str,
        icon: QMessageBox.Icon,
        headline: str,
        details: str,
    ) -> None:
        super().__init__(parent)
        self.clicked_button: QPushButton | None = None
        self.setWindowTitle(title)
        self.setModal(True)
        self.resize(860, 560)

        layout = QVBoxLayout(self)
        header = QHBoxLayout()
        icon_label = QLabel()
        icon_label.setPixmap(
            self.style().standardIcon(_message_box_standard_icon(icon)).pixmap(
                QSize(48, 48)
            )
        )
        icon_label.setAlignment(
            Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignHCenter
        )
        header.addWidget(icon_label)

        headline_label = QLabel(headline)
        headline_label.setWordWrap(True)
        headline_label.setObjectName("DialogHeadline")
        headline_label.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Preferred,
        )
        header.addWidget(headline_label, 1)
        layout.addLayout(header)

        self.details = QPlainTextEdit()
        self.details.setReadOnly(True)
        self.details.setPlainText(details)
        self.details.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.details.setMinimumSize(720, 300)
        self.details.setMaximumHeight(430)
        self.details.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        layout.addWidget(self.details, 1)

        self._actions = QHBoxLayout()
        self._actions.addStretch(1)
        layout.addLayout(self._actions)

    def add_action(self, text: str, *, accept: bool = False) -> QPushButton:
        """Add an action button and return it for caller comparison."""
        button = QPushButton(text)
        button.setObjectName("SmallActionButton")
        button.clicked.connect(
            lambda _checked=False, clicked=button, should_accept=accept: (
                self._finish(clicked, should_accept)
            )
        )
        self._actions.addWidget(button)
        return button

    def set_default_button(self, button: QPushButton) -> None:
        """Mark one button as the dialog default."""
        button.setDefault(True)
        button.setAutoDefault(True)

    @classmethod
    def information(
        cls,
        *,
        parent: QWidget,
        title: str,
        icon: QMessageBox.Icon,
        headline: str,
        details: str,
    ) -> None:
        """Show a scrollable informational dialog with a single OK button."""
        dialog = cls(
            parent=parent,
            title=title,
            icon=icon,
            headline=headline,
            details=details,
        )
        ok_button = dialog.add_action("OK", accept=True)
        dialog.set_default_button(ok_button)
        dialog.exec()

    def _finish(self, button: QPushButton, accept: bool) -> None:
        self.clicked_button = button
        if accept:
            self.accept()
        else:
            self.reject()


def _message_box_standard_icon(icon: QMessageBox.Icon) -> QStyle.StandardPixmap:
    """Return a standard pixmap matching a QMessageBox icon enum."""
    if icon == QMessageBox.Icon.Critical:
        return QStyle.StandardPixmap.SP_MessageBoxCritical
    if icon == QMessageBox.Icon.Warning:
        return QStyle.StandardPixmap.SP_MessageBoxWarning
    if icon == QMessageBox.Icon.Question:
        return QStyle.StandardPixmap.SP_MessageBoxQuestion
    return QStyle.StandardPixmap.SP_MessageBoxInformation


def _tdb_validation_prompt_text(
    diagnostics: list[Diagnostic],
    corrected_diagnostics: list[Diagnostic],
    correction_count: int,
    *,
    allow_use_as_is: bool,
) -> str:
    """Return concise text for the TDB correction confirmation dialog."""
    lines = [
        f"Available corrections: {correction_count}.",
        "",
        "Current diagnostics:",
        _format_tdb_diagnostics(diagnostics),
    ]
    if corrected_diagnostics:
        lines.extend(
            [
                "",
                "Diagnostics after applying corrections:",
                _format_tdb_diagnostics(corrected_diagnostics),
            ]
        )
    else:
        lines.extend(["", "Applying corrections clears the current diagnostics."])
    if not allow_use_as_is:
        lines.extend(
            [
                "",
                "Use As-Is is disabled because the current file has blocking errors.",
            ]
        )
    return "\n".join(lines)


def _format_tdb_diagnostics(
    diagnostics: list[Diagnostic],
    *,
    limit: int = 8,
) -> str:
    """Format TDB diagnostics for compact GUI messages."""
    if not diagnostics:
        return "No diagnostics."
    lines: list[str] = []
    for diagnostic in diagnostics[:limit]:
        prefix = diagnostic.severity.upper()
        source = diagnostic.source
        location = ""
        if source.file and source.line:
            location = f" ({Path(source.file).name}:{source.line})"
        lines.append(f"- {prefix}{location}: {diagnostic.message}")
    if len(diagnostics) > limit:
        lines.append(f"... {len(diagnostics) - limit} additional diagnostic(s).")
    return "\n".join(lines)


def _write_module_result_csv(
    module: CalculationModule,
    table: QTableWidget,
    path: str | Path,
) -> Path:
    """Write the full stored module result when available, not just the preview."""
    result = module.runtime.get("result") or module.result
    result_table = _result_table_for_module(module, result)
    if result_table is None:
        return _write_result_table_csv(table, path)

    selected_columns = (
        list(module.result_columns)
        if module.result_columns
        else _default_result_column_keys(module, result_table, result)
    )
    expanded_columns = _expand_result_column_defaults(
        selected_columns,
        result_table,
        module,
    )
    available_columns = set(result_table.available_columns())
    columns = [column for column in expanded_columns if column in available_columns]
    if not columns:
        return _write_result_table_csv(table, path)

    selected_table = result_table.select(columns)
    selected_data = selected_table.to_dict()
    headers = [
        _result_column_display_header(module.kind, column.key, column.label)
        for column in selected_table.columns
    ]
    row_count = _result_table_row_count(selected_table)
    output_path = Path(path).expanduser()
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        for row_index in range(row_count):
            writer.writerow(
                [
                    _format_result_column_value(
                        module.kind,
                        column.key,
                        _sequence_value(selected_data.get(column.key), row_index),
                    )
                    for column in selected_table.columns
                ]
            )
    return output_path


def _module_result_figure_is_active(module: CalculationModule) -> bool:
    runtime = getattr(module, "runtime", {}) or {}
    stack = runtime.get("result_view_stack")
    return (
        isinstance(stack, QStackedWidget)
        and stack.currentIndex() == 1
        and _result_figure_chart_has_series(runtime.get("result_figure_chart"))
    )


def _result_figure_chart_has_series(chart_view: Any) -> bool:
    try:
        chart = chart_view.chart()
    except Exception:
        return False
    try:
        return bool(chart.series())
    except Exception:
        return False


def _result_figure_export_size(parent: QWidget, chart_view: Any) -> QSize | None:
    """Ask for SVG export dimensions and show a centered fit preview."""
    default_size = _result_figure_default_export_size(chart_view)
    state = {
        "unit": "inch",
        "width_px": default_size.width(),
        "height_px": default_size.height(),
    }
    dialog = QDialog(parent)
    dialog.setWindowTitle("Export Figure")
    dialog.setModal(True)
    dialog.resize(680, 500)

    layout = QVBoxLayout(dialog)
    layout.setContentsMargins(18, 16, 18, 16)
    layout.setSpacing(12)

    size_row = QHBoxLayout()
    size_row.setSpacing(10)
    unit_combo = QComboBox(dialog)
    unit_combo.setObjectName("FigureExportUnit")
    unit_combo.addItems(["px", "inch", "cm"])
    unit_combo.setCurrentText(str(state["unit"]))
    width_spin = QDoubleSpinBox(dialog)
    width_spin.setObjectName("FigureExportWidth")
    height_spin = QDoubleSpinBox(dialog)
    height_spin.setObjectName("FigureExportHeight")
    for spin in (width_spin, height_spin):
        spin.setMinimumWidth(120)
    size_row.addWidget(QLabel("Unit"))
    size_row.addWidget(unit_combo)
    size_row.addWidget(QLabel("Width"))
    size_row.addWidget(width_spin)
    size_row.addWidget(QLabel("Height"))
    size_row.addWidget(height_spin)
    size_row.addStretch(1)
    layout.addLayout(size_row)

    preview = _ResultFigurePreviewLabel(chart_view, dialog)
    preview.setMinimumSize(360, 260)
    layout.addWidget(preview, 1)

    def refresh_size_editor() -> None:
        unit = str(state["unit"])
        for spin in (width_spin, height_spin):
            spin.blockSignals(True)
        try:
            for spin in (width_spin, height_spin):
                spin.setDecimals(0 if unit == "px" else 2)
                spin.setRange(
                    _result_figure_px_to_unit(120, unit),
                    _result_figure_px_to_unit(10000, unit),
                )
                spin.setSingleStep(
                    max(_result_figure_px_to_unit(40, unit), 0.01)
                )
                spin.setSuffix(f" {_result_figure_unit_suffix(unit)}")
            width_spin.setValue(
                _result_figure_px_to_unit(int(state["width_px"]), unit)
            )
            height_spin.setValue(
                _result_figure_px_to_unit(int(state["height_px"]), unit)
            )
        finally:
            for spin in (width_spin, height_spin):
                spin.blockSignals(False)

    def update_preview() -> None:
        preview.set_export_size(QSize(int(state["width_px"]), int(state["height_px"])))

    def unit_changed(unit: str) -> None:
        state["unit"] = unit
        refresh_size_editor()
        update_preview()

    def width_changed(value: float) -> None:
        state["width_px"] = _result_figure_unit_to_px(
            value,
            str(state["unit"]),
        )
        update_preview()

    def height_changed(value: float) -> None:
        state["height_px"] = _result_figure_unit_to_px(
            value,
            str(state["unit"]),
        )
        update_preview()

    unit_combo.currentTextChanged.connect(unit_changed)
    width_spin.valueChanged.connect(width_changed)
    height_spin.valueChanged.connect(height_changed)
    refresh_size_editor()
    update_preview()

    button_row = QHBoxLayout()
    button_row.addStretch(1)
    cancel_button = QPushButton("Cancel")
    cancel_button.setObjectName("SmallActionButton")
    export_button = QPushButton("Export SVG")
    export_button.setObjectName("PrimaryButton")
    cancel_button.clicked.connect(dialog.reject)
    export_button.clicked.connect(dialog.accept)
    button_row.addWidget(cancel_button)
    button_row.addWidget(export_button)
    layout.addLayout(button_row)

    if dialog.exec() != QDialog.DialogCode.Accepted:
        return None
    return QSize(int(state["width_px"]), int(state["height_px"]))


def _result_figure_default_export_size(chart_view: Any) -> QSize:
    current_size = chart_view.size() if hasattr(chart_view, "size") else QSize()
    hint_size = chart_view.sizeHint() if hasattr(chart_view, "sizeHint") else QSize()
    width_source = current_size.width() if current_size.width() > 0 else hint_size.width()
    height_source = (
        current_size.height()
        if current_size.height() > 0
        else hint_size.height()
    )
    width = max(640, width_source)
    height = max(360, height_source if height_source > 0 else 480)
    return QSize(width, height)


class _ResultFigurePreviewLabel(QWidget):
    """Preview widget that centers a scaled figure in the available popup area."""

    def __init__(self, chart_view: Any, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._chart_view = chart_view
        self._export_size = QSize()
        self._preview_pixmap = QPixmap()
        self._preview_rect = QRect()
        self.setObjectName("FigureExportPreview")
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )

    def set_export_size(self, export_size: QSize) -> None:
        self._export_size = QSize(export_size)
        self._refresh_preview()

    def resizeEvent(self, event: Any) -> None:
        super().resizeEvent(event)
        self._refresh_preview()

    def _refresh_preview(self) -> None:
        available_rect = self.contentsRect()
        if (
            self._export_size.width() <= 0
            or self._export_size.height() <= 0
            or available_rect.width() <= 0
            or available_rect.height() <= 0
        ):
            self._preview_pixmap = QPixmap()
            self._preview_rect = QRect()
            self.update()
            return
        preview_rect = _result_figure_preview_rect(
            self._export_size,
            available_rect.size(),
        )
        preview_rect.translate(available_rect.topLeft())
        self._preview_pixmap = _result_figure_preview_pixmap(
            self._chart_view,
            self._export_size,
            target_size=preview_rect.size(),
        )
        self._preview_rect = preview_rect
        self.update()

    def paintEvent(self, event: Any) -> None:
        super().paintEvent(event)
        painter = QPainter(self)
        try:
            if self._preview_pixmap.isNull():
                painter.drawText(
                    self.contentsRect(),
                    int(Qt.AlignmentFlag.AlignCenter),
                    "Preview unavailable",
                )
            else:
                painter.drawPixmap(self._preview_rect, self._preview_pixmap)
        finally:
            painter.end()


def _result_figure_unit_to_px(value: float, unit: str) -> int:
    unit = _normalized_result_figure_unit(unit)
    if unit == "inch":
        pixels = float(value) * _RESULT_FIGURE_EXPORT_DPI
    elif unit == "cm":
        pixels = float(value) / 2.54 * _RESULT_FIGURE_EXPORT_DPI
    else:
        pixels = float(value)
    return max(1, int(round(pixels)))


def _result_figure_px_to_unit(pixels: int, unit: str) -> float:
    unit = _normalized_result_figure_unit(unit)
    if unit == "inch":
        return float(pixels) / _RESULT_FIGURE_EXPORT_DPI
    if unit == "cm":
        return float(pixels) / _RESULT_FIGURE_EXPORT_DPI * 2.54
    return float(pixels)


def _format_result_figure_size_value(pixels: int, unit: str) -> str:
    unit = _normalized_result_figure_unit(unit)
    value = _result_figure_px_to_unit(pixels, unit)
    if unit == "px":
        return f"{int(round(value))} px"
    return f"{value:.2f} {_result_figure_unit_suffix(unit)}"


def _result_figure_unit_suffix(unit: str) -> str:
    unit = _normalized_result_figure_unit(unit)
    if unit == "inch":
        return "inch"
    return unit


def _normalized_result_figure_unit(unit: str) -> str:
    unit = str(unit or "px").strip().lower()
    return unit if unit in {"px", "inch", "cm"} else "px"


def _result_figure_preview_rect(export_size: QSize, available_size: QSize) -> QRect:
    if (
        export_size.width() <= 0
        or export_size.height() <= 0
        or available_size.width() <= 0
        or available_size.height() <= 0
    ):
        return QRect()
    preview_size = export_size.scaled(
        available_size,
        Qt.AspectRatioMode.KeepAspectRatio,
    )
    x = (available_size.width() - preview_size.width()) // 2
    y = (available_size.height() - preview_size.height()) // 2
    return QRect(x, y, preview_size.width(), preview_size.height())


def _result_figure_preview_pixmap(
    chart_view: Any,
    export_size: QSize,
    *,
    target_size: QSize | None = None,
) -> QPixmap:
    preview_size = target_size or export_size.scaled(
        QSize(420, 260),
        Qt.AspectRatioMode.KeepAspectRatio,
    )
    if preview_size.width() <= 0 or preview_size.height() <= 0:
        return QPixmap()
    pixmap = QPixmap(preview_size)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    try:
        _paint_scaled_result_figure(
            chart_view,
            painter,
            preview_size,
            source_size=export_size,
        )
    finally:
        painter.end()
    return pixmap


def _write_result_figure_svg(
    chart_view: Any,
    path: str | Path,
    export_size: QSize,
) -> Path:
    """Write the active figure as SVG at the requested pixel dimensions."""
    if QSvgGenerator is None:
        raise RuntimeError("This Qt build does not include SVG export support.")
    if export_size.width() <= 0 or export_size.height() <= 0:
        raise ValueError("Figure export size must be positive.")

    output_path = Path(path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    generator = QSvgGenerator()
    generator.setFileName(str(output_path))
    generator.setSize(export_size)
    generator.setViewBox(QRect(0, 0, export_size.width(), export_size.height()))
    generator.setTitle("Equilipy figure")
    generator.setDescription("Exported from the Equilipy GUI.")

    old_size = chart_view.size() if hasattr(chart_view, "size") else QSize()
    restore_size = old_size.width() > 0 and old_size.height() > 0
    if hasattr(chart_view, "resize"):
        chart_view.resize(export_size)
        QApplication.processEvents()

    painter = QPainter(generator)
    try:
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        if (
            hasattr(chart_view, "size")
            and chart_view.size().width() == export_size.width()
            and chart_view.size().height() == export_size.height()
        ):
            chart_view.render(painter)
        else:
            _paint_scaled_result_figure(chart_view, painter, export_size)
    finally:
        painter.end()
        if restore_size and hasattr(chart_view, "resize"):
            chart_view.resize(old_size)
            QApplication.processEvents()
    return output_path


def _paint_scaled_result_figure(
    chart_view: Any,
    painter: QPainter,
    target_size: QSize,
    *,
    source_size: QSize | None = None,
) -> None:
    old_size = chart_view.size() if hasattr(chart_view, "size") else QSize()
    restore_size = old_size.width() > 0 and old_size.height() > 0
    if source_size is not None and hasattr(chart_view, "resize"):
        chart_view.resize(source_size)
        QApplication.processEvents()
    source_size = _result_figure_source_size(chart_view)
    painter.save()
    try:
        painter.scale(
            target_size.width() / source_size.width(),
            target_size.height() / source_size.height(),
        )
        chart_view.render(painter)
    finally:
        painter.restore()
        if restore_size and hasattr(chart_view, "resize"):
            chart_view.resize(old_size)
            QApplication.processEvents()


def _result_figure_source_size(chart_view: Any) -> QSize:
    size = chart_view.size() if hasattr(chart_view, "size") else QSize()
    if size.width() <= 0 or size.height() <= 0:
        size = chart_view.sizeHint() if hasattr(chart_view, "sizeHint") else QSize()
    return QSize(max(1, size.width()), max(1, size.height()))
