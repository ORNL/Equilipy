"""Calculation workspace construction for the main GUI window."""

from __future__ import annotations

from PySide6.QtCore import QSize, Qt, Slot
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QStackedWidget,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from equilipy.gui.assets import calculation_icon_path

from ..log_capture import GuiLogCapture
from ..units import AMOUNT_UNIT_OPTIONS
from ..widgets.log import CalculationLog
from ..widgets.navigation import ClearableTreeWidget


def _combo_box(values: list[str]) -> QComboBox:
    combo = QComboBox()
    combo.addItems(values)
    return combo


SIDEBAR_MAX_WIDTH = 900


class CalculationWorkspaceMixin:
    """Build the calculation workspace shell for ``DatabaseEditorWindow``."""

    def _build_calculation_workspace(self) -> QWidget:
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setObjectName("WorkspaceSplitter")
        splitter.setChildrenCollapsible(False)

        splitter.addWidget(self._build_calculation_sidebar())
        splitter.addWidget(self._build_calculation_center())
        splitter.addWidget(self._build_chat_sidebar("calculation"))
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 0)
        splitter.setSizes([330, 850, 320])
        return splitter

    def _build_calculation_center(self) -> QWidget:
        vertical_splitter = QSplitter(Qt.Orientation.Vertical)
        vertical_splitter.setChildrenCollapsible(False)
        vertical_splitter.addWidget(self._build_calculation_editor())
        vertical_splitter.addWidget(self._build_calculation_log_panel())
        vertical_splitter.setStretchFactor(0, 4)
        vertical_splitter.setStretchFactor(1, 1)
        vertical_splitter.setSizes([620, 180])
        return vertical_splitter

    def _build_calculation_sidebar(self) -> QWidget:
        self.calculation_sidebar = QFrame()
        self.calculation_sidebar.setObjectName("SidebarPanel")
        self.calculation_sidebar.setMinimumWidth(290)
        self.calculation_sidebar.setMaximumWidth(SIDEBAR_MAX_WIDTH)
        sidebar_layout = QVBoxLayout(self.calculation_sidebar)
        sidebar_layout.setContentsMargins(18, 18, 18, 18)
        sidebar_layout.setSpacing(12)

        session_action_row = QHBoxLayout()
        session_action_row.setSpacing(8)

        self.calculation_active_directory_button = QPushButton("Open")
        self.calculation_active_directory_button.setObjectName("SmallActionButton")
        self.calculation_active_directory_button.setToolTip(
            "Set active calculation directory"
        )
        self.calculation_active_directory_button.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self.calculation_active_directory_button.clicked.connect(
            self._set_active_calculation_directory
        )
        self.calculation_active_directory_button.setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu
        )
        self.calculation_active_directory_button.customContextMenuRequested.connect(
            self._show_active_calculation_directory_menu
        )
        session_action_row.addWidget(self.calculation_active_directory_button)

        new_session_button = QPushButton("New Session")
        new_session_button.setObjectName("SmallActionButton")
        new_session_button.setToolTip("Create new calculation session")
        new_session_button.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        new_session_button.clicked.connect(self._add_calculation_session)
        session_action_row.addWidget(new_session_button)

        save_session_button = QPushButton("Save")
        save_session_button.setObjectName("SmallActionButton")
        save_session_button.setToolTip("Save current calculation session")
        save_session_button.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        save_session_button.clicked.connect(
            lambda _checked=False: self._save_current_session()
        )
        session_action_row.addWidget(save_session_button)
        sidebar_layout.addLayout(session_action_row)

        module_row = QHBoxLayout()
        module_row.setSpacing(10)
        module_row.addWidget(
            self._calculation_module_button(
                "equilibrium",
                "Equilibrium",
            )
        )
        module_row.addWidget(
            self._calculation_module_button(
                "solidification",
                "Solidification",
            )
        )
        module_row.addWidget(
            self._calculation_module_button(
                "nucleation",
                "Nucleation",
            )
        )
        module_row.addStretch(1)
        sidebar_layout.addLayout(module_row)

        self.calculation_tree = ClearableTreeWidget()
        self.calculation_tree.setObjectName("DatabaseTree")
        self.calculation_tree.setHeaderHidden(True)
        self.calculation_tree.setIconSize(QSize(18, 18))
        self.calculation_tree.setIndentation(12)
        self.calculation_tree.itemSelectionChanged.connect(
            self._remember_calculation_selection
        )
        self.calculation_tree.itemChanged.connect(self._rename_calculation_item)
        self.calculation_tree.setContextMenuPolicy(
            Qt.ContextMenuPolicy.CustomContextMenu
        )
        self.calculation_tree.customContextMenuRequested.connect(
            self._show_calculation_tree_menu
        )
        sidebar_layout.addWidget(self.calculation_tree, 1)
        return self.calculation_sidebar

    def _build_calculation_log_panel(self) -> QWidget:
        self.calculation_log_tabs = QTabWidget()
        self.calculation_log_tabs.setObjectName("InspectorTabs")
        self.calculation_log = CalculationLog()
        self.calculation_log.status_sink = self._handle_status_message
        self.calculation_log.setReadOnly(True)
        self.calculation_log.setObjectName("SourcePreview")
        self.calculation_log.setPlaceholderText("Calculation log will appear here.")
        self.calculation_log_tabs.addTab(self.calculation_log, "Log")
        self.calculation_log_tabs.setVisible(False)
        return self.calculation_log_tabs

    def attach_log_capture(self, log_capture: GuiLogCapture) -> None:
        """Attach process stdout/stderr capture to the bottom Log panel."""
        self.log_capture = log_capture
        log_capture.text_received.connect(self._append_captured_log_text)
        self._update_log_capture_file()

    @Slot(str)
    def _append_captured_log_text(self, text: str) -> None:
        if self.calculation_log is None:
            return
        self.calculation_log.append_stream_text(text)

    def _write_status_to_log_capture(self, text: str) -> None:
        if self.log_capture is not None:
            self.log_capture.write_message(text)

    def _handle_status_message(self, text: str) -> None:
        self._write_status_to_log_capture(text)
        self._set_status_bar_message(text)

    def _update_log_capture_file(self) -> None:
        if self.log_capture is None:
            return
        if self.active_calculation_directory is None:
            self.log_capture.set_log_file(None)
            return
        self.log_capture.set_log_file(
            self.active_calculation_directory / "equilipy_gui_log.txt"
        )

    def _build_calculation_editor(self) -> QWidget:
        self.calculation_editor = QFrame()
        self.calculation_editor.setObjectName("EditorPanel")
        panel_layout = QVBoxLayout(self.calculation_editor)
        panel_layout.setContentsMargins(22, 20, 22, 20)
        panel_layout.setSpacing(14)

        header = QHBoxLayout()
        header.setSpacing(16)
        title_block = QVBoxLayout()
        title_block.setSpacing(4)

        self.calculation_badge.setText("Calculation")
        self.calculation_badge.setObjectName("EditorBadge")
        self.calculation_badge.setReadOnly(True)
        self.calculation_badge.editingFinished.connect(
            self._rename_current_session_from_heading
        )
        title_block.addWidget(self.calculation_badge)

        self.calculation_title.setText("Welcome to Equilipy")
        self.calculation_title.setObjectName("EditorTitle")
        self.calculation_title.setReadOnly(True)
        self.calculation_title.editingFinished.connect(
            self._rename_current_module_from_heading
        )
        title_block.addWidget(self.calculation_title)

        self.calculation_subtitle.setText(
            "Open an active directory, create a session, select a module "
            "for calculations"
        )
        self.calculation_subtitle.setObjectName("EditorSubtitle")
        self.calculation_subtitle.setWordWrap(True)
        self.calculation_subtitle.setMinimumWidth(0)
        self.calculation_subtitle.setSizePolicy(
            QSizePolicy.Policy.Ignored,
            QSizePolicy.Policy.Preferred,
        )
        title_block.addWidget(self.calculation_subtitle)
        self.calculation_elements_label.setObjectName("EditorSubtitle")
        self.calculation_elements_label.setWordWrap(True)
        self.calculation_elements_label.setMinimumWidth(0)
        self.calculation_elements_label.setSizePolicy(
            QSizePolicy.Policy.Ignored,
            QSizePolicy.Policy.Preferred,
        )
        title_block.addWidget(self.calculation_elements_label)
        title_block.addWidget(self._build_calculation_unit_row())
        header.addLayout(title_block, 1)

        header_actions = QHBoxLayout()
        header_actions.setSpacing(8)
        self.calculation_unit_button = QPushButton("Unit")
        self.calculation_unit_button.setObjectName("SmallActionButton")
        self.calculation_unit_button.clicked.connect(
            lambda _checked=False: self._open_session_unit_dialog()
        )
        header_actions.addWidget(self.calculation_unit_button)

        self.calculation_phases_button = QPushButton("Phases")
        self.calculation_phases_button.setObjectName("SmallActionButton")
        self.calculation_phases_button.clicked.connect(
            lambda _checked=False: self._open_session_phase_dialog()
        )
        header_actions.addWidget(self.calculation_phases_button)

        self.calculation_load_session_button = QPushButton("Load")
        self.calculation_load_session_button.setObjectName("PrimaryButton")
        self.calculation_load_session_button.clicked.connect(
            lambda _checked=False: self._load_session_from_file()
        )
        header_actions.addWidget(self.calculation_load_session_button)

        self.calculation_save_session_button = QPushButton("Save")
        self.calculation_save_session_button.setObjectName("PrimaryButton")
        self.calculation_save_session_button.clicked.connect(
            lambda _checked=False: self._save_current_session()
        )
        header_actions.addWidget(self.calculation_save_session_button)

        self.calculation_calculate_all_button = QPushButton("Calculate All")
        self.calculation_calculate_all_button.setObjectName("PrimaryButton")
        self.calculation_calculate_all_button.clicked.connect(
            lambda _checked=False: self._calculate_all_modules()
        )
        header_actions.addWidget(self.calculation_calculate_all_button)

        self.calculation_load_module_button = QPushButton("Load")
        self.calculation_load_module_button.setObjectName("SmallActionButton")
        self.calculation_load_module_button.clicked.connect(
            lambda _checked=False: self._load_module_from_file()
        )
        header_actions.addWidget(self.calculation_load_module_button)

        self.calculation_save_module_button = QPushButton("Save")
        self.calculation_save_module_button.setObjectName("SmallActionButton")
        self.calculation_save_module_button.clicked.connect(
            lambda _checked=False: self._save_current_module()
        )
        header_actions.addWidget(self.calculation_save_module_button)

        self.calculation_calculate_button = QPushButton("Calculate")
        self.calculation_calculate_button.setObjectName("PrimaryButton")
        self.calculation_calculate_button.clicked.connect(
            lambda _checked=False: self._calculate_current_module()
        )
        header_actions.addWidget(self.calculation_calculate_button)
        header.addLayout(header_actions)
        panel_layout.addLayout(header)

        self.calculation_stack = QStackedWidget()
        self.calculation_stack.setObjectName("CalculationStack")
        self._set_calculation_page(
            "empty",
            self._calculation_overview_page(None),
            make_current=True,
        )
        self._set_calculation_header_actions("empty")
        panel_layout.addWidget(self.calculation_stack, 1)
        return self.calculation_editor

    def _build_calculation_unit_row(self) -> QFrame:
        row = QFrame()
        row.setObjectName("HeaderUnitRow")
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 4, 0, 0)
        layout.setSpacing(10)

        self.calculation_temperature_unit = _combo_box(["K", "C", "F", "R"])
        self.calculation_pressure_unit = _combo_box(["atm", "psi", "bar", "Pa", "kPa"])
        self.calculation_amount_unit = _combo_box(AMOUNT_UNIT_OPTIONS)
        for combo in (
            self.calculation_temperature_unit,
            self.calculation_pressure_unit,
            self.calculation_amount_unit,
        ):
            combo.currentTextChanged.connect(self._store_current_module_units)

        layout.addWidget(QLabel("Temperature"))
        layout.addWidget(self.calculation_temperature_unit)
        layout.addWidget(QLabel("Pressure"))
        layout.addWidget(self.calculation_pressure_unit)
        layout.addWidget(QLabel("Amount"))
        layout.addWidget(self.calculation_amount_unit)
        layout.addStretch(1)
        row.setVisible(False)
        self.calculation_unit_row = row
        return row

    def _calculation_module_button(
        self,
        module_id: str,
        tooltip: str,
    ) -> QToolButton:
        button = QToolButton()
        button.setObjectName("CalculationModuleButton")
        button.setIcon(QIcon(str(calculation_icon_path(module_id))))
        button.setIconSize(QSize(34, 34))
        button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        button.setFixedSize(46, 46)
        button.setCursor(Qt.CursorShape.PointingHandCursor)
        button.setToolTip(tooltip)
        button.clicked.connect(
            lambda _checked=False, module=module_id: self._open_calculation_module(
                module
            )
        )
        return button
