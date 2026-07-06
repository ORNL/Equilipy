"""Main PySide6 database editor window."""

from __future__ import annotations

# ruff: noqa: F401,F405,I001

import base64
import copy
import csv
import json
import os
import pickle
import re
import shutil
from dataclasses import asdict
from importlib.metadata import PackageNotFoundError, version
from itertools import product
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator

import numpy as np
from PySide6.QtCore import (
    QEvent,
    QItemSelectionModel,
    QModelIndex,
    QPersistentModelIndex,
    QSize,
    Qt,
    QThread,
    QTimer,
)
from PySide6.QtGui import (
    QAction,
    QBrush,
    QColor,
    QFontDatabase,
    QIcon,
    QKeySequence,
)
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
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
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from equilipy.composition import expand_condition_species as _expand_condition_species
from equilipy.database_ir import (
    ConstituentSet,
    DatabaseIR,
    Diagnostic,
    Element,
    FunctionDefinition,
    GibbsRange,
    Parameter,
    Phase,
    SourceRef,
    Species,
    load_database_ir,
    write_eqdb,
    write_tdb,
)
from equilipy.database_ir.tdb_syntax import parameter_command_prefix_match
from equilipy.results import ResultTable
from equilipy.results.equilib import EquilibResult
from equilipy.results.scheil import ScheilResult
from equilipy.utils import G2HSCp, HSCp2G, NeumanKoppHSCp

from .assets import (
    calculation_icon_path,
    equilipy_logo_path,
    math_font_path,
)
from .assistant import AssistantPanel
from .calculation.workspace import CalculationWorkspaceLogicMixin
from .calculation.editor import CalculationWorkspaceMixin
from .calculation.state import (
    CalculationDatabase,
    CalculationModule,
    CalculationSession,
)
from .calculation.workers import (
    CalculationCanceled,
    CalculationResultReceiver,
    CalculationWorker,
)
from .database.workspace import DatabaseWorkspaceLogicMixin
from .database.editor import DatabaseWorkspaceMixin
from .log_capture import GuiLogCapture
from .models import (
    CompoundPhasePayload,
    CompoundRecord,
    DatabaseTreeModel,
    DiagnosticsModel,
    ThermoRangePayload,
    function_error_icon,
    function_warning_icon,
    thermo_range_warning_key,
)
from .style import app_style_sheet as _style_sheet
from .units import AMOUNT_UNIT_OPTIONS as _AMOUNT_UNIT_OPTIONS
from .widgets import forms as _forms
from .widgets.composition import (
    COMPOSITION_INPUT_MIN_HEIGHT,
    COMPOSITION_ROW_HEIGHT,
    CompositionEditor,
)
from .widgets.log import CalculationLog
from .widgets.navigation import (
    ClearableTreeView,
    ClearableTreeWidget,
    LayoutToggleButton,
)
from .widgets.status_bar import AppStatusBar

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
_latest_status_text = _forms.latest_status_text
_read_only_line = _forms.read_only_line
_read_only_text = _forms.read_only_text
_table_section = _forms.table_section
_text_section = _forms.text_section
_untitled_form_section = _forms.untitled_form_section

_RESULT_OVERVIEW_ITEM_ID = "__result_overview__"
_THERMO_COEFFICIENT_BOX_WIDTH = 170
_THERMO_POWER_BOX_WIDTH = 80
_THERMO_TEMPERATURE_BOX_WIDTH = round(_THERMO_COEFFICIENT_BOX_WIDTH * 0.5)
_GIBBS_COEFFICIENT_SLOTS = (0, 1, 2, 3, 4, 5, 6, 8, 10, 12)
_FUNCTION_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)")
_PORTABLE_TDB_FUNCTION_NAME_WIDTH = 8
_MAX_TDB_FUNCTION_NAME_WIDTH = 25
_GIBBS_CONTINUITY_TOLERANCE = 1e-5
_MATH_FONT_FAMILY: str | None = None


def _application_version() -> str:
    """Return the installed package version, falling back to pyproject.toml."""
    try:
        return version("equilipy")
    except PackageNotFoundError:
        pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
        try:
            match = re.search(
                r'(?m)^version\s*=\s*"([^"]+)"',
                pyproject.read_text(encoding="utf-8"),
            )
        except OSError:
            match = None
        return match.group(1) if match is not None else "unknown"


class DatabaseEditorWindow(
    CalculationWorkspaceLogicMixin,
    DatabaseWorkspaceLogicMixin,
    DatabaseWorkspaceMixin,
    CalculationWorkspaceMixin,
    QMainWindow,
):
    """PySide6 GUI shell for DatabaseIR browsing and editing."""

    def __init__(self, database: DatabaseIR):
        super().__init__()
        self.databases = [database]
        self.database = database
        self.current_database = database
        self.setWindowTitle("Equilipy")
        self.setObjectName("DatabaseEditorWindow")
        self.setMinimumSize(1120, 720)
        self._font_size_delta = 0
        _load_math_font_family()
        self.setStyleSheet(_style_sheet(self._font_size_delta))

        self.function_continuity_warnings: dict[int, list[str]] = {}
        self.function_name_warnings: dict[int, list[str]] = {}
        self.function_range_warnings: dict[str, list[str]] = {}
        self.compound_range_warnings: dict[str, list[str]] = {}
        self.function_name_errors: dict[int, list[str]] = {}
        self.tree_model = self._new_database_tree_model()
        self.tree = ClearableTreeView()
        self.last_tree_index: QPersistentModelIndex | None = None
        self.editor_title = QLineEdit()
        self.editor_subtitle = QLabel()
        self.editor_subtitle_control = QWidget()
        self.editor_subtitle_control_layout: QHBoxLayout | None = None
        self.editor_badge = QLabel()
        self._editable_editor_title_record: FunctionDefinition | Phase | None = None
        self.editor_form = QWidget()
        self.editor_form_layout = QVBoxLayout(self.editor_form)
        self._active_database_thermo_editor: Any | None = None
        self.source_preview = QPlainTextEdit()
        self.workspace_stack = QStackedWidget()
        self.database_sidebar: QFrame | None = None
        self.calculation_sidebar: QFrame | None = None
        self.database_editor: QFrame | None = None
        self.calculation_editor: QFrame | None = None
        self.calculation_badge = QLineEdit()
        self.calculation_title = QLineEdit()
        self.calculation_subtitle = QLabel()
        self.calculation_elements_label = QLabel()
        self.calculation_unit_row: QFrame | None = None
        self.calculation_temperature_unit: QComboBox | None = None
        self.calculation_pressure_unit: QComboBox | None = None
        self.calculation_amount_unit: QComboBox | None = None
        self.calculation_load_session_button: QPushButton | None = None
        self.calculation_save_session_button: QPushButton | None = None
        self.calculation_load_module_button: QPushButton | None = None
        self.calculation_save_module_button: QPushButton | None = None
        self.calculation_calculate_button: QPushButton | None = None
        self.calculation_unit_button: QPushButton | None = None
        self.calculation_phases_button: QPushButton | None = None
        self.calculation_calculate_all_button: QPushButton | None = None
        self.calculation_active_directory_button: QPushButton | None = None
        self.active_calculation_directory: Path | None = None
        self.calculation_tree: ClearableTreeWidget | None = None
        self.calculation_stack: QStackedWidget | None = None
        self.calculation_pages: dict[str, QWidget] = {}
        self.result_overview_collapsed_modules: dict[str, set[str]] = {}
        self.last_calculation_item: QTreeWidgetItem | None = None
        self.calculation_session_counter = 0
        self._updating_calculation_tree = False
        self._updating_calculation_units = False
        self._stale_calculation_source_paths: set[str] = set()
        self.calculation_sessions: dict[str, CalculationSession] = {}
        self.database_chat_sidebar: QFrame | None = None
        self.calculation_chat_sidebar: QFrame | None = None
        self.database_inspector_tabs: QTabWidget | None = None
        self.database_diagnostics_table: QTableView | None = None
        self.calculation_log_tabs: QTabWidget | None = None
        self.calculation_log: CalculationLog | None = None
        self.log_capture: GuiLogCapture | None = None
        self.status_bar: AppStatusBar | None = None
        self.status_symbol: QLabel | None = None
        self.status_text: QLabel | None = None
        self._calculation_in_progress = False
        self._calculation_queue: list[dict[str, Any]] = []
        self._calculation_thread_refs: list[
            tuple[QThread, CalculationWorker, CalculationResultReceiver]
        ] = []
        self._transient_message_boxes: list[QMessageBox] = []
        self._detached_result_windows: dict[str, QWidget] = {}
        self.database_validation_diagnostics: dict[int, list[Diagnostic]] = {}

        self._build_actions()
        self._build_window()
        app = QApplication.instance()
        if app is not None:
            app.installEventFilter(self)
        self._show_database_summary()

    def closeEvent(self, event) -> None:
        """Shut down calculation threads before the window is destroyed.

        Destroying a running QThread aborts the whole process, so a close
        during a calculation must cancel the worker and wait for the thread
        to finish (or terminate it as a last resort) before accepting.
        """
        if self._calculation_in_progress and self._calculation_thread_refs:
            response = QMessageBox.question(
                self,
                "Quit Equilipy",
                "A calculation is still running. Quit anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if response != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
        self._calculation_queue.clear()
        for thread, worker, _receiver in list(self._calculation_thread_refs):
            worker.cancel()
            thread.quit()
        for thread, _worker, _receiver in list(self._calculation_thread_refs):
            if thread.isRunning() and not thread.wait(5000):
                thread.terminate()
                thread.wait(1000)
        super().closeEvent(event)

    def _build_window(self) -> None:
        root = QWidget()
        root.setObjectName("RootSurface")
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(0, 0, 0, 0)
        root_layout.setSpacing(0)
        root_layout.addWidget(self._build_header())

        self.workspace_stack.addWidget(self._build_calculation_workspace())
        self.workspace_stack.addWidget(self._build_database_workspace())
        root_layout.addWidget(self.workspace_stack, 1)
        root_layout.addWidget(self._build_status_bar())

        self.setCentralWidget(root)
        _apply_button_cursor(root)
        self._set_status_bar_message("Ready.")

    def _build_status_bar(self) -> QWidget:
        self.status_bar = AppStatusBar()
        self.status_symbol = self.status_bar.status_symbol
        self.status_text = self.status_bar.status_text
        return self.status_bar

    def _set_status_bar_message(self, message: str, level: str | None = None) -> None:
        if self.status_bar is None:
            return
        self.status_bar.set_message(message, level)

    def _build_header(self) -> QWidget:
        header = QFrame()
        header.setObjectName("AppHeader")
        layout = QGridLayout(header)
        layout.setContentsMargins(40, 0, 18, 0)
        layout.setHorizontalSpacing(12)
        layout.setVerticalSpacing(0)

        mode_selector = QFrame()
        mode_selector.setObjectName("ModeSelector")
        mode_layout = QHBoxLayout(mode_selector)
        mode_layout.setContentsMargins(4, 4, 4, 4)
        mode_layout.setSpacing(4)

        calculation_button = self._mode_button("Calculation")
        database_button = self._mode_button("Database")
        calculation_button.setChecked(True)

        self.mode_group = QButtonGroup(self)
        self.mode_group.setExclusive(True)
        self.mode_group.addButton(calculation_button, 0)
        self.mode_group.addButton(database_button, 1)
        self.mode_group.idClicked.connect(self._handle_workspace_mode_clicked)

        mode_layout.addWidget(calculation_button)
        mode_layout.addWidget(database_button)

        brand = self._build_brand()
        layout_toggles = self._build_layout_toggles()

        layout.addWidget(brand, 0, 0, Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(mode_selector, 0, 1, Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(layout_toggles, 0, 2, Qt.AlignmentFlag.AlignRight)
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(1, 0)
        layout.setColumnStretch(2, 1)
        return header

    def _build_brand(self) -> QWidget:
        brand = QFrame()
        brand.setObjectName("HeaderBrand")
        layout = QHBoxLayout(brand)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)

        logo = QSvgWidget(str(equilipy_logo_path()))
        logo.setObjectName("HeaderLogo")
        logo.setFixedSize(70, 70)

        brand_text = QFrame()
        brand_text.setObjectName("HeaderBrandText")
        brand_text_layout = QVBoxLayout(brand_text)
        brand_text_layout.setContentsMargins(0, 0, 0, 0)
        brand_text_layout.setSpacing(0)

        name = QLabel("Equilipy")
        name.setObjectName("HeaderBrandName")
        name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        version_label = QLabel(f"v{_application_version()}")
        version_label.setObjectName("HeaderBrandVersion")
        version_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        brand_text_layout.addWidget(name, alignment=Qt.AlignmentFlag.AlignHCenter)
        brand_text_layout.addWidget(
            version_label,
            alignment=Qt.AlignmentFlag.AlignHCenter,
        )
        layout.addWidget(logo)
        layout.addWidget(brand_text, alignment=Qt.AlignmentFlag.AlignVCenter)
        return brand

    def _build_layout_toggles(self) -> QWidget:
        toggles = QFrame()
        toggles.setObjectName("LayoutToggleGroup")
        layout = QHBoxLayout(toggles)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        left = LayoutToggleButton("left", "Toggle left sidebar")
        bottom = LayoutToggleButton("bottom", "Toggle bottom panel", checked=False)
        right = LayoutToggleButton("right", "Toggle assistant panel", checked=False)
        left.clicked.connect(self._toggle_left_sidebar)
        bottom.clicked.connect(self._toggle_bottom_panel)
        right.clicked.connect(self._toggle_right_chat)

        layout.addWidget(left)
        layout.addWidget(bottom)
        layout.addWidget(right)
        return toggles

    def _handle_workspace_mode_clicked(self, mode_index: int) -> None:
        """Switch workspace mode."""
        self.workspace_stack.setCurrentIndex(mode_index)
        mode_name = "Calculation" if mode_index == 0 else "Database"
        self._set_status_bar_message(f"{mode_name} workspace active.")

    def _show_future_feature_warning(self, feature: str, message: str) -> None:
        """Inform users that an unfinished feature is intentionally unavailable."""
        dialog = QMessageBox(
            QMessageBox.Icon.Information,
            f"{feature} Coming Soon",
            f"{message}\n\nThis feature will be implemented in a future release.",
            QMessageBox.StandardButton.Ok,
            self,
        )
        dialog.setModal(False)
        dialog.finished.connect(
            lambda _result, box=dialog: self._transient_message_boxes.remove(box)
            if box in self._transient_message_boxes
            else None
        )
        self._transient_message_boxes.append(dialog)
        dialog.open()

    def _build_chat_sidebar(self, workspace: str) -> QFrame:
        chat = AssistantPanel(self, workspace)
        if workspace == "database":
            self.database_chat_sidebar = chat
        else:
            self.calculation_chat_sidebar = chat
        return chat

    def _build_actions(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        tools_menu = self.menuBar().addMenu("&Tools")
        view_menu = self.menuBar().addMenu("&View")

        validate_action = QAction("Validate DatabaseIR", self)
        validate_action.triggered.connect(self._show_validation_summary)
        tools_menu.addAction(validate_action)

        zoom_in_action = QAction("Zoom In", self)
        zoom_in_action.setShortcuts(
            [QKeySequence("Ctrl+="), QKeySequence("Ctrl++")]
        )
        zoom_in_action.triggered.connect(lambda: self._change_font_zoom(1))
        view_menu.addAction(zoom_in_action)

        zoom_out_action = QAction("Zoom Out", self)
        zoom_out_action.setShortcut(QKeySequence("Ctrl+-"))
        zoom_out_action.triggered.connect(lambda: self._change_font_zoom(-1))
        view_menu.addAction(zoom_out_action)

        reset_zoom_action = QAction("Reset Font Size", self)
        reset_zoom_action.setShortcut(QKeySequence("Ctrl+0"))
        reset_zoom_action.triggered.connect(self._reset_font_zoom)
        view_menu.addAction(reset_zoom_action)

        open_tdb_action = QAction("Open TDB...", self)
        open_tdb_action.triggered.connect(self._open_database_ir_file)
        file_menu.addAction(open_tdb_action)

        save_tdb_action = QAction("Save Structured TDB...", self)
        save_tdb_action.triggered.connect(self._export_database_tdb)
        file_menu.addAction(save_tdb_action)

        save_eqdb_action = QAction("Save Equilipy DB...", self)
        save_eqdb_action.triggered.connect(self._save_database_eqdb)
        file_menu.addAction(save_eqdb_action)

        placeholder_actions = [
            ("Open DAT...", file_menu),
            ("Section by Elements...", tools_menu),
            ("Compare/Merge...", tools_menu),
        ]
        for label, menu in placeholder_actions:
            action = QAction(label, self)
            action.triggered.connect(
                lambda _checked=False, text=label: self._placeholder(text)
            )
            menu.addAction(action)

    def _change_font_zoom(self, steps: int) -> None:
        """Increase or decrease the GUI stylesheet font scale."""
        updated_delta = max(-4, min(10, self._font_size_delta + steps))
        if updated_delta == self._font_size_delta:
            return
        self._font_size_delta = updated_delta
        self._apply_font_zoom()

    def _reset_font_zoom(self) -> None:
        """Reset the GUI font scale to its default size."""
        if self._font_size_delta == 0:
            return
        self._font_size_delta = 0
        self._apply_font_zoom()

    def _apply_font_zoom(self) -> None:
        """Reapply the stylesheet with the current font-size delta."""
        self.setStyleSheet(_style_sheet(self._font_size_delta))

    def _mode_button(self, text: str) -> QPushButton:
        button = QPushButton(text)
        button.setObjectName("ModeButton")
        button.setCheckable(True)
        button.setCursor(Qt.CursorShape.PointingHandCursor)
        button.setMinimumWidth(118)
        return button

    def eventFilter(self, watched, event) -> bool:
        """Restore the tree highlight when the editor pane is clicked."""
        if event.type() == QEvent.Type.ChildAdded:
            child = event.child()
            if isinstance(child, QWidget):
                _apply_button_cursor(child)
        if event.type() == QEvent.Type.MouseButtonPress:
            if self._is_database_editor_widget(watched):
                self._restore_tree_selection()
            elif self._is_calculation_editor_widget(watched):
                self._restore_calculation_selection()
        return super().eventFilter(watched, event)

    def _is_database_editor_widget(self, watched) -> bool:
        if self.workspace_stack.currentIndex() != 1 or self.database_editor is None:
            return False
        return isinstance(watched, QWidget) and (
            watched is self.database_editor
            or self.database_editor.isAncestorOf(watched)
        )

    def _is_calculation_editor_widget(self, watched) -> bool:
        if self.workspace_stack.currentIndex() != 0 or self.calculation_editor is None:
            return False
        return isinstance(watched, QWidget) and (
            watched is self.calculation_editor
            or self.calculation_editor.isAncestorOf(watched)
        )

from .compat_exports import *  # noqa: E402,F403
