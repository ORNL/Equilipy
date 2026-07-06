"""Calculation session tree and overview pages."""

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
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
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
from equilipy.results import ResultTable
from equilipy.results.equilib import EquilibResult
from equilipy.results.scheil import ScheilResult

from equilipy.gui.assets import calculation_icon_path, ornl_logo_path
from equilipy.gui.compat_exports import *
from equilipy.gui.units import AMOUNT_UNIT_OPTIONS as _AMOUNT_UNIT_OPTIONS
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


def _session_unit_combo_box(values: list[str]) -> QComboBox:
    combo = QComboBox()
    combo.addItems(values)
    return combo


def _developers_section() -> QFrame:
    """Return the welcome-page developer table."""
    section = _form_section("Lead Developer")
    lead_name = QLabel("<b>Sun Yong Kwon</b> (kwons@ornl.gov)")
    lead_name.setObjectName("LeadDeveloperName")
    lead_name.setTextFormat(Qt.TextFormat.RichText)
    lead_name.setWordWrap(True)

    affiliation_row = QWidget()
    affiliation_row.setObjectName("LeadDeveloperAffiliationRow")
    affiliation_layout = QHBoxLayout(affiliation_row)
    affiliation_layout.setContentsMargins(0, 0, 0, 0)
    affiliation_layout.setSpacing(8)
    logo = QLabel()
    logo.setObjectName("LeadDeveloperAffiliationLogo")
    logo.setPixmap(QIcon(str(ornl_logo_path())).pixmap(QSize(22, 22)))
    logo.setFixedSize(24, 24)
    affiliation = QLabel("Oak Ridge National Laboratory")
    affiliation.setObjectName("LeadDeveloperAffiliation")
    affiliation.setWordWrap(True)
    affiliation_layout.addWidget(logo)
    affiliation_layout.addWidget(affiliation, 1)

    developers_heading = QLabel("Developers")
    developers_heading.setObjectName("DevelopersTableTitle")
    heading_font = developers_heading.font()
    heading_font.setBold(True)
    developers_heading.setFont(heading_font)

    headers = ("Name", "Role", "Affiliation")
    rows = [
        ("Eric Thibodeau", "Developer", "Independent"),
        ("Ying Yang", "Developer", "Oak Ridge National Laboratory"),
        ("Alex Plotkowski", "Collaborator", "Oak Ridge National Laboratory"),
        ("Sam Reeve", "Contributor", "Oak Ridge National Laboratory"),
    ]
    table = QTableWidget(len(rows), len(headers))
    table.setObjectName("DevelopersTable")
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setAlternatingRowColors(False)
    table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    for row_index, row in enumerate(rows):
        for column_index, value in enumerate(row):
            item = QTableWidgetItem(str(value))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row_index, column_index, item)

    table.horizontalHeader().setStretchLastSection(True)
    table.resizeColumnsToContents()
    _fit_table_columns_to_headers(table)
    _disable_inner_scrollbars(table)
    _fit_table_to_rows(table)
    section.layout().addWidget(lead_name)
    section.layout().addWidget(affiliation_row)
    section.layout().addWidget(developers_heading)
    section.layout().addWidget(table)
    return section


def _sponsoring_organization_section() -> QFrame:
    """Return the compact welcome-page sponsoring organization box."""
    section = _form_section("Sponsoring Organization")
    sponsor_text = _read_only_text(
        "US Department of Energy (DOE), Office of Critical Minerals "
        "and Energy Innovation (CMEI), Transportation Technologies "
        "Office (TTO)",
        80,
    )
    sponsor_text.setMaximumHeight(80)
    section.layout().addWidget(sponsor_text)
    return section


class CalculationSessionOverviewMixin:
    """Overview, tree, selection, and rename controller methods."""

    def _calculation_overview_page(
        self,
        session: CalculationSession | None,
    ) -> QWidget:
        page = QScrollArea()
        page.setObjectName("EditorScroll")
        page.setWidgetResizable(True)
        page.setFrameShape(QFrame.Shape.NoFrame)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(14)

        if session is None:
            layout.addWidget(_developers_section())
            layout.addWidget(_sponsoring_organization_section())
            layout.addStretch(1)
            page.setWidget(content)
            return page

        layout.addWidget(self._session_database_section(session))

        selected_dbs = [db for db in session.databases.values() if db.selected]
        if not selected_dbs:
            layout.addWidget(
                _text_section(
                    "Database Overview",
                    "No database selected.\n"
                    "Load a .dat database and check its box to see\n"
                    "available elements and phases here.",
                )
            )
        else:
            for db in selected_dbs:
                layout.addWidget(self._database_overview_tables(db))

        layout.addStretch(1)
        page.setWidget(content)
        return page

    def _session_result_overview_page(
        self,
        session: CalculationSession,
    ) -> QWidget:
        page = QScrollArea()
        page.setObjectName("EditorScroll")
        page.setWidgetResizable(True)
        page.setFrameShape(QFrame.Shape.NoFrame)

        content = QWidget()
        content.setObjectName("SessionResultOverview")
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(14)

        result_modules = [
            module
            for module in session.modules.values()
            if _module_has_displayable_result(module)
        ]

        if not result_modules:
            layout.addWidget(
                _text_section(
                    "Results Overview",
                    "No calculation results are available yet.\n"
                    "Run Equilib or Solidification modules to populate this page.",
                )
            )
        else:
            for module in result_modules:
                layout.addWidget(
                    self._module_result_overview_section(session, module)
                )

        layout.addStretch(1)
        page.setWidget(content)
        return page

    def _module_result_overview_section(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> QFrame:
        collapsed = module.id in self.result_overview_collapsed_modules.get(
            session.id,
            set(),
        )

        table = QTableWidget(0, 0)
        table.setObjectName("FormTable")
        table.verticalHeader().setVisible(False)
        table.horizontalHeader().setStretchLastSection(True)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        _disable_inner_scrollbars(table)

        columns_button = QPushButton("Columns")
        columns_button.setObjectName("SmallActionButton")
        columns_button.clicked.connect(
            lambda: self._open_result_columns_dialog(module, table)
        )
        show_all_button = QPushButton("Show All")
        show_all_button.setObjectName("SmallActionButton")
        show_all_button.clicked.connect(
            lambda: self._show_all_module_results(module, table)
        )
        export_button = QPushButton("Export")
        export_button.setObjectName("PrimaryButton")
        export_button.clicked.connect(
            lambda: self._export_module_results(module, table)
        )
        section, body_layout = _collapsible_form_section(
            f"Result {module.name}",
            [columns_button, show_all_button, export_button],
            collapsed=collapsed,
        )
        section.setProperty("module_id", module.id)
        section.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        section.customContextMenuRequested.connect(
            lambda position, selected_module=module, widget=section: (
                self._show_result_overview_subwindow_menu(
                    session,
                    selected_module,
                    widget,
                    position,
                )
            )
        )

        result = module.runtime.get("result") or module.result
        if result is not None:
            try:
                _populate_module_results(module, table, result)
            except Exception:
                _restore_table_payload(table, module.results_table_payload)
        else:
            _restore_table_payload(table, module.results_table_payload)

        body_layout.addWidget(table)
        return section

    def _session_database_section(self, session: CalculationSession) -> QFrame:
        section = QFrame()
        section.setObjectName("FormSection")
        layout = QVBoxLayout(section)
        layout.setContentsMargins(16, 14, 16, 16)
        layout.setSpacing(12)

        header = QHBoxLayout()
        header.setSpacing(8)
        title = QLabel("Database")
        title.setObjectName("SectionTitle")
        header.addWidget(title)
        open_button = QToolButton()
        open_button.setObjectName("IconActionButton")
        open_button.setIcon(
            self.style().standardIcon(QStyle.StandardPixmap.SP_DirOpenIcon)
        )
        open_button.setToolTip("Load database")
        open_button.clicked.connect(lambda: self._load_session_database(session.id))
        header.addWidget(open_button)
        header.addStretch(1)
        layout.addLayout(header)

        database_list = QWidget()
        database_list.setObjectName("DatabaseListViewport")
        database_list.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed
        )
        grid = QGridLayout(database_list)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(8)
        grid.setSizeConstraint(QGridLayout.SizeConstraint.SetMinimumSize)

        if not session.databases:
            empty_label = QLabel("No database loaded")
            empty_label.setObjectName("DatabaseListLabel")
            empty_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            grid.addWidget(empty_label, 0, 0, 1, 2)
        else:
            for index, database in enumerate(session.databases.values()):
                grid.addWidget(
                    self._session_database_row(session.id, database),
                    index // 2,
                    index % 2,
                )
        grid.setColumnStretch(0, 1)
        grid.setColumnStretch(1, 1)
        row_count = max(1, (len(session.databases) + 1) // 2)
        database_list.setMinimumHeight(row_count * 48 + max(0, row_count - 1) * 8)

        layout.addWidget(database_list)
        return section

    def _session_database_row(
        self,
        session_id: str,
        database: CalculationDatabase,
    ) -> QWidget:
        row = QFrame()
        row.setObjectName("DatabaseListItem")
        layout = QHBoxLayout(row)
        layout.setContentsMargins(8, 6, 8, 6)
        layout.setSpacing(8)

        selected = QCheckBox()
        selected.setFixedSize(22, 22)
        selected.setChecked(database.selected)
        selected.toggled.connect(
            lambda checked, db_id=database.id: self._set_session_database_selected(
                session_id,
                db_id,
                checked,
            )
        )
        layout.addWidget(selected, 0, Qt.AlignmentFlag.AlignVCenter)
        layout.addSpacing(4)

        remove = QToolButton()
        remove.setObjectName("DatabaseRemoveButton")
        remove.setFixedSize(22, 22)
        remove.setIconSize(QSize(18, 18))
        remove.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_TrashIcon))
        remove.setToolTip("Remove database")
        remove.clicked.connect(
            lambda _checked=False, db_id=database.id: self._remove_session_database(
                session_id,
                db_id,
            )
        )
        layout.addWidget(remove, 0, Qt.AlignmentFlag.AlignVCenter)

        name = QLabel(database.name)
        name.setObjectName("DatabaseListLabel")
        name.setToolTip(database.path)
        name.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        layout.addWidget(name, 1)
        return row

    def _set_calculation_page(
        self,
        key: str,
        page: QWidget | None = None,
        make_current: bool = True,
    ) -> None:
        if self.calculation_stack is None:
            return
        if page is not None and key in self.calculation_pages:
            self._close_detached_result_window(key)
            old = self.calculation_pages.pop(key)
            self.calculation_stack.removeWidget(old)
            old.deleteLater()
        if page is not None:
            _apply_button_cursor(page)
            self.calculation_pages[key] = page
            self.calculation_stack.addWidget(page)
        if make_current and key in self.calculation_pages:
            self.calculation_stack.setCurrentWidget(self.calculation_pages[key])

    def _remove_calculation_page(self, key: str) -> None:
        if self.calculation_stack is None or key not in self.calculation_pages:
            return
        self._close_detached_result_window(key)
        page = self.calculation_pages.pop(key)
        self.calculation_stack.removeWidget(page)
        page.deleteLater()

    def _overview_key(self, session_id: str) -> str:
        return f"{session_id}:overview"

    def _result_overview_key(self, session_id: str) -> str:
        return f"{session_id}:result-overview"

    def _module_key(self, module_id: str) -> str:
        return module_id

    def _set_calculation_header_actions(
        self,
        mode: str,
        session: CalculationSession | None = None,
    ) -> None:
        show_module_actions = mode == "module"
        show_session_io = mode == "overview" and session is not None
        show_session_actions = (
            mode == "overview"
            and session is not None
            and bool(session.databases)
        )
        show_calculate_all = (
            show_session_actions
            and sum(
                module.kind != "nucleation"
                for module in session.modules.values()
            )
            > 1
        )
        for button in (
            self.calculation_load_module_button,
            self.calculation_save_module_button,
            self.calculation_calculate_button,
        ):
            if button is not None:
                button.setVisible(show_module_actions)
        for button in (
            self.calculation_unit_button,
            self.calculation_phases_button,
        ):
            if button is not None:
                button.setVisible(show_session_actions)
        for button in (
            self.calculation_load_session_button,
            self.calculation_save_session_button,
        ):
            if button is not None:
                button.setVisible(show_session_io)
        if self.calculation_calculate_all_button is not None:
            self.calculation_calculate_all_button.setVisible(show_calculate_all)

    def _update_calculation_heading(
        self,
        session: CalculationSession | None,
        module_name: str = "",
    ) -> None:
        if session is None:
            self._set_header_unit_controls(None)
            self._set_calculation_header_actions("empty")
            self.calculation_badge.setReadOnly(True)
            self.calculation_title.setReadOnly(True)
            self.calculation_badge.setText("Calculation")
            self.calculation_title.setText("Welcome to Equilipy")
            self.calculation_subtitle.setText(
                "Open an active directory, create a session, select a module "
                "for calculations"
            )
            self.calculation_elements_label.setText("")
            return
        self.calculation_badge.setReadOnly(False)
        self.calculation_badge.setText(session.name)
        if module_name:
            self._set_calculation_header_actions("module", session)
            self.calculation_title.setReadOnly(False)
            self.calculation_title.setText(module_name)
            self.calculation_subtitle.setText(
                f"{self._session_database_summary(session)} | "
                f"{self._current_module_phase_count(session)} phase(s)"
            )
            self.calculation_elements_label.setText(
                self._session_elements_summary(session)
            )
        else:
            self._set_header_unit_controls(None)
            self._set_calculation_header_actions("overview", session)
            self.calculation_title.setReadOnly(True)
            self.calculation_title.setText("Session Overview")
            self.calculation_subtitle.setText(
                self._session_database_summary(session)
            )
            self.calculation_elements_label.setText(
                self._session_elements_summary(session)
            )

    def _update_current_calculation_heading(
        self,
        session: CalculationSession,
    ) -> None:
        module_name = ""
        module: CalculationModule | None = None
        if self.calculation_stack is not None:
            current = self.calculation_stack.currentWidget()
            if current is not None and current.property("session_id") == session.id:
                if bool(current.property("result_overview")):
                    self._set_result_overview_heading(session)
                    return
                module_id = current.property("module_id")
                if isinstance(module_id, str):
                    module = session.modules.get(module_id)
                    if module is not None:
                        module_name = module.name
        self._update_calculation_heading(session, module_name)
        self._set_header_unit_controls(module)

    def _session_database_summary(self, session: CalculationSession) -> str:
        if not session.databases:
            return "No database loaded"
        active_database = self._active_session_database(session)
        if active_database is None:
            return "No active database"
        return f"{active_database.name} active"

    def _session_elements_summary(self, session: CalculationSession) -> str:
        all_elements: set[str] = set()
        for db in session.databases.values():
            if not db.selected or db.database is None:
                continue
            elements = db.database.get("cElementNameCS")
            if elements is not None:
                for e in elements:
                    if isinstance(e, str) and e.strip():
                        all_elements.add(e.strip())
        if not all_elements:
            return ""
        return f"Available elements: {', '.join(sorted(all_elements))}"

    def _database_overview_text(self, session: CalculationSession) -> str:
        selected_dbs = [
            db for db in session.databases.values() if db.selected
        ]
        if not selected_dbs:
            return (
                "No database selected.\n\n"
                "Load a .dat database and check its box to see\n"
                "available elements and phases here."
            )
        lines: list[str] = []
        for db in selected_dbs:
            data = db.database
            if data is None:
                continue
            lines.append(f"── {db.name} ──")

            elements = data.get("cElementNameCS")
            if elements is not None:
                names = sorted(
                    e.strip() for e in elements if isinstance(e, str) and e.strip()
                )
                lines.append(f"  Elements ({len(names)}): {', '.join(names)}")

            phase_names = data.get("cSolnPhaseNameCS")
            species_per_phase = data.get("nSolnPhaseCS")
            phase_types = data.get("cSolnPhaseTypeCS")
            if phase_names is not None:
                soln_lines: list[str] = []
                for i, name in enumerate(phase_names):
                    name_str = name.strip() if isinstance(name, str) else ""
                    if not name_str:
                        continue
                    ptype = ""
                    if phase_types is not None and i < len(phase_types):
                        pt = phase_types[i]
                        ptype = pt.strip() if isinstance(pt, str) else ""
                    nsp = ""
                    if species_per_phase is not None and i < len(species_per_phase):
                        nsp = f", {species_per_phase[i]} species"
                    tag = f" ({ptype}{nsp})" if ptype or nsp else ""
                    soln_lines.append(f"    {name_str}{tag}")
                lines.append(f"  Solution phases ({len(soln_lines)}):")
                lines.extend(soln_lines)

            pure_names = _pure_compound_names(data)
            lines.append(f"  Pure compounds: {len(pure_names)}")
            if pure_names:
                lines.append(f"    {', '.join(pure_names)}")

            lines.append("")
        return "\n".join(lines).rstrip()

    def _database_overview_tables(self, db: CalculationDatabase) -> QFrame:
        data = db.database or {}
        section = _form_section(db.name)

        # Elements line
        elements = data.get("cElementNameCS")
        if elements is not None:
            names = [e.strip() for e in elements if isinstance(e, str) and e.strip()]
            el_label = QLabel(f"Elements ({len(names)}): {', '.join(names)}")
            el_label.setObjectName("DatabaseListLabel")
            el_label.setWordWrap(True)
            el_label.setMinimumWidth(0)
            el_label.setSizePolicy(
                QSizePolicy.Policy.Ignored,
                QSizePolicy.Policy.Preferred,
            )
            section.layout().addWidget(el_label)

        # Solution phases table
        soln_rows = _solution_phase_overview_rows(data)
        if soln_rows:
            section.layout().addWidget(
                _inline_table(
                    f"Solution Phases ({len(soln_rows)})",
                    ("#", "Phase", "Model", "Structure", "Species (count)"),
                    soln_rows,
                )
            )

        # Pure compounds table
        pure_names = _pure_compound_names(data)
        if pure_names:
            pure_rows = [(i + 1, name) for i, name in enumerate(pure_names)]
            section.layout().addWidget(
                _inline_table(
                    f"Pure Compounds ({len(pure_names)})",
                    ("#", "Compound"),
                    pure_rows,
                )
            )

        return section

    def _current_module_phase_count(self, session: CalculationSession) -> int:
        if self.calculation_stack is None:
            return 0
        current = self.calculation_stack.currentWidget()
        if current is None or current.property("session_id") != session.id:
            return 0
        module_id = current.property("module_id")
        module = session.modules.get(module_id) if isinstance(module_id, str) else None
        return len(module.phase_names) if module else 0

    def _show_calculation_overview(
        self,
        session: CalculationSession | None,
    ) -> None:
        if session is not None:
            self._select_calculation_session_item(session.id)
        self._update_calculation_heading(session)
        key = self._overview_key(session.id) if session else "empty"
        page = self._calculation_overview_page(session)
        if session is not None:
            page.setProperty("session_id", session.id)
        self._set_calculation_page(
            key,
            page,
            make_current=True,
        )

    def _show_session_result_overview(
        self,
        session: CalculationSession,
        *,
        reset_collapsed: bool = False,
    ) -> None:
        if reset_collapsed:
            self.result_overview_collapsed_modules[session.id] = set()
        overview_item = self._ensure_result_overview_tree_item(session.id)
        if self.calculation_tree is not None and overview_item is not None:
            previous_blocked = self.calculation_tree.blockSignals(True)
            self.calculation_tree.setCurrentItem(overview_item)
            self.calculation_tree.blockSignals(previous_blocked)
            self.last_calculation_item = overview_item
        self._set_result_overview_heading(session)
        key = self._result_overview_key(session.id)
        page = self._session_result_overview_page(session)
        page.setProperty("session_id", session.id)
        page.setProperty("result_overview", True)
        self._set_calculation_page(key, page, make_current=True)

    def _set_result_overview_heading(self, session: CalculationSession) -> None:
        self._update_calculation_heading(session)
        self._set_header_unit_controls(None)
        self._set_calculation_header_actions("empty")
        self.calculation_title.setReadOnly(True)
        self.calculation_title.setText("Overview")
        result_count = sum(
            1
            for module in session.modules.values()
            if _module_has_displayable_result(module)
        )
        self.calculation_subtitle.setText(
            f"{result_count} result(s) | {self._session_database_summary(session)}"
        )
        self.calculation_elements_label.setText(
            self._session_elements_summary(session)
        )

    def _add_calculation_session(
        self,
        _checked: bool = False,
    ) -> CalculationSession:
        session = self._create_calculation_session()
        self._add_session_tree_item(session)
        self._show_calculation_overview(session)
        return session

    def _create_calculation_session(
        self,
        name: str | None = None,
    ) -> CalculationSession:
        self.calculation_session_counter += 1
        session_id = f"session-{self.calculation_session_counter}"
        session_name = self._unique_calculation_session_name(
            name or f"Session#{self.calculation_session_counter}"
        )
        session = CalculationSession(id=session_id, name=session_name)
        self.calculation_sessions[session.id] = session
        return session

    def _add_session_tree_item(
        self,
        session: CalculationSession,
    ) -> QTreeWidgetItem:
        item = QTreeWidgetItem([session.name])
        item.setData(0, Qt.ItemDataRole.UserRole, session.id)
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
        item.setExpanded(True)
        if self.calculation_tree is not None:
            self.calculation_tree.addTopLevelItem(item)
            self.calculation_tree.setCurrentItem(item)
        self.last_calculation_item = item
        return item

    def _ensure_result_overview_tree_item(
        self,
        session_id: str,
    ) -> QTreeWidgetItem | None:
        item = self._find_session_item(session_id)
        if item is None:
            return None
        for row in range(item.childCount()):
            child = item.child(row)
            if child.data(0, Qt.ItemDataRole.UserRole + 2) == "result_overview":
                return child
        child = QTreeWidgetItem(["Overview"])
        child.setIcon(0, QIcon(str(calculation_icon_path("result_overview"))))
        child.setData(0, Qt.ItemDataRole.UserRole, session_id)
        child.setData(0, Qt.ItemDataRole.UserRole + 1, _RESULT_OVERVIEW_ITEM_ID)
        child.setData(0, Qt.ItemDataRole.UserRole + 2, "result_overview")
        child.setFlags(child.flags() & ~Qt.ItemFlag.ItemIsEditable)
        item.insertChild(0, child)
        item.setExpanded(True)
        return child

    def _show_calculation_tree_menu(self, position) -> None:
        if self.calculation_tree is None:
            return
        item = self.calculation_tree.itemAt(position)
        if item is None:
            return

        self.calculation_tree.setCurrentItem(item)
        menu = self._calculation_tree_context_menu(item)
        if menu is not None:
            menu.exec(self.calculation_tree.viewport().mapToGlobal(position))

    def _calculation_tree_context_menu(
        self,
        item: QTreeWidgetItem,
    ) -> QMenu | None:
        menu = QMenu(self)
        session_id = self._session_id_for_item(item)
        session = self.calculation_sessions.get(session_id)
        module_id = item.data(0, Qt.ItemDataRole.UserRole + 1)
        item_kind = item.data(0, Qt.ItemDataRole.UserRole + 2)
        if session is not None and item_kind == "result_overview":
            refresh = QAction(
                QIcon(str(calculation_icon_path("result_overview"))),
                "Refresh Overview",
                self,
            )
            delete = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete",
                self,
            )
            refresh.triggered.connect(
                lambda _checked=False, target_session=session: (
                    self._show_session_result_overview(
                        target_session,
                        reset_collapsed=True,
                    )
                )
            )
            delete.triggered.connect(
                lambda _checked=False, target_session=session, target_item=item: (
                    self._delete_result_overview_item(target_session, target_item)
                )
            )
            menu.addAction(refresh)
            menu.addAction(delete)
            return menu
        if session is not None and not isinstance(module_id, str):
            add_equilib = QAction(
                QIcon(str(calculation_icon_path("equilibrium"))),
                "Add Equilib",
                self,
            )
            add_solidification = QAction(
                QIcon(str(calculation_icon_path("solidification"))),
                "Add Solidification",
                self,
            )
            overview = QAction(
                QIcon(str(calculation_icon_path("result_overview"))),
                "Add Result Overview",
                self,
            )
            add_equilib.triggered.connect(
                lambda _checked=False, target_session=session: (
                    self._add_calculation_module_from_menu(
                        target_session,
                        "equilibrium",
                    )
                )
            )
            add_solidification.triggered.connect(
                lambda _checked=False, target_session=session: (
                    self._add_calculation_module_from_menu(
                        target_session,
                        "solidification",
                    )
                )
            )
            overview.triggered.connect(
                lambda _checked=False, target_session=session: (
                    self._show_session_result_overview(
                        target_session,
                        reset_collapsed=True,
                    )
                )
            )
            menu.addAction(add_equilib)
            menu.addAction(add_solidification)
            menu.addAction(overview)
            menu.addSeparator()
        duplicate = QAction(
            self._theme_icon("edit-copy", QStyle.StandardPixmap.SP_FileIcon),
            "Duplicate",
            self,
        )
        delete = QAction(
            self._theme_icon("user-trash", QStyle.StandardPixmap.SP_TrashIcon),
            "Delete",
            self,
        )
        duplicate.triggered.connect(
            lambda _checked=False, tree_item=item: self._duplicate_calculation_item(
                tree_item
            )
        )
        delete.triggered.connect(
            lambda _checked=False, tree_item=item: self._delete_calculation_item(
                tree_item
            )
        )
        menu.addAction(duplicate)
        menu.addAction(delete)
        return menu

    def _add_calculation_module_from_menu(
        self,
        session: CalculationSession,
        module_kind: str,
    ) -> CalculationModule:
        module = self._add_calculation_module(session, module_kind)
        self._show_calculation_module(session, module)
        return module

    def _show_result_overview_subwindow_menu(
        self,
        session: CalculationSession,
        module: CalculationModule,
        widget: QWidget,
        position,
    ) -> None:
        menu = QMenu(self)
        collapsed = module.id in self.result_overview_collapsed_modules.get(
            session.id,
            set(),
        )
        collapse = QAction(
            "Expand" if collapsed else "Collapse",
            self,
        )
        collapse.triggered.connect(
            lambda _checked=False: self._toggle_result_overview_module(
                session,
                module.id,
            )
        )
        menu.addAction(collapse)
        menu.exec(widget.mapToGlobal(position))

    def _toggle_result_overview_module(
        self,
        session: CalculationSession,
        module_id: str,
    ) -> None:
        collapsed = self.result_overview_collapsed_modules.setdefault(session.id, set())
        if module_id in collapsed:
            collapsed.remove(module_id)
        else:
            collapsed.add(module_id)
        self._show_session_result_overview(session)

    def _theme_icon(
        self,
        theme_name: str,
        fallback: QStyle.StandardPixmap,
    ) -> QIcon:
        icon = QIcon.fromTheme(theme_name)
        if icon.isNull():
            icon = self.style().standardIcon(fallback)
        return icon

    def _duplicate_calculation_item(
        self,
        item: QTreeWidgetItem,
    ) -> CalculationSession | CalculationModule | None:
        session_id = self._session_id_for_item(item)
        session = self.calculation_sessions.get(session_id)
        if session is None:
            return None
        if item.data(0, Qt.ItemDataRole.UserRole + 2) == "result_overview":
            self._show_session_result_overview(session, reset_collapsed=True)
            return None

        module_id = item.data(0, Qt.ItemDataRole.UserRole + 1)
        if isinstance(module_id, str):
            module = session.modules.get(module_id)
            if module is None:
                return None
            return self._duplicate_calculation_module(session, module)
        return self._duplicate_calculation_session(session)

    def _duplicate_calculation_session(
        self,
        session: CalculationSession,
    ) -> CalculationSession:
        existing_names = [item.name for item in self.calculation_sessions.values()]
        new_session = self._create_calculation_session(
            self._copy_name(session.name, existing_names)
        )
        new_session.temperature_unit = session.temperature_unit
        new_session.pressure_unit = session.pressure_unit
        new_session.amount_unit = session.amount_unit
        new_session.default_phase_names = list(session.default_phase_names)
        new_session.result_column_defaults = {
            str(kind): list(columns)
            for kind, columns in session.result_column_defaults.items()
        }

        for database in session.databases.values():
            new_session.database_counter += 1
            database_id = f"{new_session.id}:database:{new_session.database_counter}"
            new_session.databases[database_id] = CalculationDatabase(
                id=database_id,
                name=database.name,
                path=database.path,
                database=database.database,
                selected=database.selected,
            )

        session_item = self._add_session_tree_item(new_session)
        for module in session.modules.values():
            new_module = self._copy_calculation_module_state(module, new_session)
            self._add_session_module_item(new_session.id, new_module)

        if self.calculation_tree is not None:
            self.calculation_tree.setCurrentItem(session_item)
        self._show_calculation_overview(new_session)
        return new_session

    def _duplicate_calculation_module(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> CalculationModule:
        existing_names = [item.name for item in session.modules.values()]
        new_module = self._copy_calculation_module_state(
            module,
            session,
            name=self._copy_name(module.name, existing_names),
        )
        self._add_session_module_item(session.id, new_module)
        item = self._find_module_item(session.id, new_module.id)
        if item is not None:
            self._set_calculation_tree_text(item, new_module.name)
        self._show_calculation_module(session, new_module)
        return new_module

    def _copy_calculation_module_state(
        self,
        module: CalculationModule,
        session: CalculationSession,
        *,
        name: str | None = None,
    ) -> CalculationModule:
        count = session.module_counts.get(module.kind, 0) + 1
        session.module_counts[module.kind] = count
        new_module = CalculationModule(
            id=f"{session.id}:{module.kind}:{count}",
            kind=module.kind,
            name=self._unique_calculation_module_name(session, name or module.name),
            phase_names=list(module.phase_names),
            temperature_unit=module.temperature_unit,
            pressure_unit=module.pressure_unit,
            amount_unit=module.amount_unit,
            calculation_type=module.calculation_type,
            solidification_model=module.solidification_model,
            nucleation_undercooling=dict(module.nucleation_undercooling),
            batch_condition={
                key: list(values) for key, values in module.batch_condition.items()
            },
            batch_condition_path=module.batch_condition_path,
            batch_cpu_count=module.batch_cpu_count,
            transition_search=module.transition_search,
            start_from_liquidus=module.start_from_liquidus,
            result_columns=list(module.result_columns),
        )
        session.modules[new_module.id] = new_module
        return new_module

    def _delete_calculation_item(
        self,
        item: QTreeWidgetItem,
        *,
        confirm: bool = True,
    ) -> bool:
        session_id = self._session_id_for_item(item)
        session = self.calculation_sessions.get(session_id)
        if session is None:
            return False
        if item.data(0, Qt.ItemDataRole.UserRole + 2) == "result_overview":
            self._delete_result_overview_item(session, item)
            return False

        module_id = item.data(0, Qt.ItemDataRole.UserRole + 1)
        if isinstance(module_id, str):
            return self._delete_calculation_module_item(
                session,
                module_id,
                item,
                confirm=confirm,
            )
        return self._delete_calculation_session_item(session, item, confirm=confirm)

    def _delete_result_overview_item(
        self,
        session: CalculationSession,
        item: QTreeWidgetItem,
    ) -> None:
        parent = item.parent()
        if parent is not None:
            parent.removeChild(item)
        self._remove_calculation_page(self._result_overview_key(session.id))
        self.result_overview_collapsed_modules.pop(session.id, None)
        if self.last_calculation_item is item:
            self.last_calculation_item = parent
        if self.calculation_tree is not None and parent is not None:
            self.calculation_tree.setCurrentItem(parent)
        self._show_calculation_overview(session)

    def _delete_calculation_module_item(
        self,
        session: CalculationSession,
        module_id: str,
        item: QTreeWidgetItem,
        *,
        confirm: bool,
    ) -> bool:
        module = session.modules.get(module_id)
        if module is None:
            return False
        if confirm and not self._confirm_delete("Delete Module", module.name):
            return False

        self._remember_stale_module_source_path(session, module)
        session.modules.pop(module_id, None)
        self.result_overview_collapsed_modules.setdefault(session.id, set()).discard(
            module_id
        )
        self._remove_calculation_page(self._module_key(module_id))
        self._remove_calculation_page(self._result_overview_key(session.id))
        parent = item.parent()
        if parent is not None:
            parent.removeChild(item)
            if self.last_calculation_item is item:
                self.last_calculation_item = parent
            if self.calculation_tree is not None:
                self.calculation_tree.setCurrentItem(parent)
        self._show_calculation_overview(session)
        return True

    def _delete_calculation_session_item(
        self,
        session: CalculationSession,
        item: QTreeWidgetItem,
        *,
        confirm: bool,
    ) -> bool:
        if confirm and not self._confirm_delete("Delete Session", session.name):
            return False

        self._remember_stale_calculation_source_path(session.source_path)
        for module in session.modules.values():
            self._remember_stale_module_source_path(session, module)

        self._remove_calculation_page(self._overview_key(session.id))
        self._remove_calculation_page(self._result_overview_key(session.id))
        self.result_overview_collapsed_modules.pop(session.id, None)
        for module_id in list(session.modules):
            self._remove_calculation_page(self._module_key(module_id))
        self.calculation_sessions.pop(session.id, None)

        next_session = None
        if self.calculation_tree is not None:
            index = self.calculation_tree.indexOfTopLevelItem(item)
            if index >= 0:
                self.calculation_tree.takeTopLevelItem(index)
                if self.calculation_tree.topLevelItemCount() > 0:
                    next_index = min(
                        index,
                        self.calculation_tree.topLevelItemCount() - 1,
                    )
                    next_item = self.calculation_tree.topLevelItem(next_index)
                    self.calculation_tree.setCurrentItem(next_item)
                    next_id = self._session_id_for_item(next_item)
                    next_session = self.calculation_sessions.get(next_id)

        self.last_calculation_item = (
            self.calculation_tree.currentItem()
            if self.calculation_tree is not None
            else None
        )
        self._show_calculation_overview(next_session)
        return True

    def _confirm_delete(self, title: str, name: str) -> bool:
        dialog = QDialog(self)
        dialog.setWindowTitle(title)
        dialog.setModal(True)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(14)

        message = QLabel(f"Delete {name}?")
        message.setObjectName("DatabaseListLabel")
        layout.addWidget(message)

        actions = QHBoxLayout()
        actions.setSpacing(8)
        yes = QPushButton("Yes")
        yes.setObjectName("SmallActionButton")
        cancel = QPushButton("Cancel")
        cancel.setObjectName("SmallActionButton")
        cancel.setDefault(True)
        cancel.setAutoDefault(True)
        yes.setAutoDefault(False)
        yes.clicked.connect(dialog.accept)
        cancel.clicked.connect(dialog.reject)

        actions.addWidget(yes)
        actions.addWidget(cancel)
        layout.addLayout(actions)

        return dialog.exec() == QDialog.DialogCode.Accepted

    def _copy_name(
        self,
        name: str,
        existing_names: list[str],
    ) -> str:
        base_name = f"{name} Copy"
        if base_name not in existing_names:
            return base_name
        index = 2
        while f"{base_name} {index}" in existing_names:
            index += 1
        return f"{base_name} {index}"

    def _unique_calculation_session_name(
        self,
        name: str,
        *,
        exclude_session_id: str | None = None,
    ) -> str:
        existing_names = [
            session.name
            for session in self.calculation_sessions.values()
            if session.id != exclude_session_id
        ]
        return _unique_numbered_name(name, existing_names, fallback="Session")

    def _unique_calculation_module_name(
        self,
        session: CalculationSession,
        name: str,
        *,
        exclude_module_id: str | None = None,
    ) -> str:
        existing_names = [
            module.name
            for module in session.modules.values()
            if module.id != exclude_module_id
        ]
        return _unique_numbered_name(name, existing_names, fallback="Module")

    def _open_calculation_module(self, module_id: str) -> None:
        if self._is_unavailable_calculation_module(module_id):
            self._show_unavailable_calculation_module(module_id)
            return

        session = self._selected_session_for_new_module()
        if session is None:
            return

        module = self._add_calculation_module(session, module_id)
        self._show_calculation_module(session, module)

    def _selected_session_for_new_module(self) -> CalculationSession | None:
        if self.calculation_tree is None:
            return None
        selected = self.calculation_tree.selectedItems()
        if not selected:
            QMessageBox.information(
                self,
                "Select Session",
                "Select a session in the left tree before adding a module.",
            )
            return None
        item = selected[0]
        session_id = self._session_id_for_item(item)
        return self.calculation_sessions.get(session_id)

    def _add_calculation_module(
        self,
        session: CalculationSession,
        module_kind: str,
    ) -> CalculationModule:
        count = session.module_counts.get(module_kind, 0) + 1
        session.module_counts[module_kind] = count
        module_name = self._unique_calculation_module_name(
            session,
            f"{self._module_label(module_kind)}#{count}",
        )
        module = CalculationModule(
            id=f"{session.id}:{module_kind}:{count}",
            kind=module_kind,
            name=module_name,
            temperature_unit=session.temperature_unit,
            pressure_unit=session.pressure_unit,
            amount_unit=session.amount_unit,
            result_columns=list(
                session.result_column_defaults.get(module_kind, [])
            ),
        )
        session.modules[module.id] = module
        self._add_session_module_item(session.id, module)
        return module

    def _show_calculation_module(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> None:
        if self._is_unavailable_calculation_module(module.kind):
            self._show_unavailable_calculation_module(module.kind)
            self._show_calculation_overview(session)
            return

        key = self._module_key(module.id)
        if key not in self.calculation_pages:
            page = self._calculation_module_page(session, module)
            page.setProperty("session_id", session.id)
            page.setProperty("module_id", module.id)
            page.setProperty("module_kind", module.kind)
            self._set_calculation_page(key, page, make_current=False)
        self._update_calculation_heading(session, module.name)
        self._set_header_unit_controls(None if module.kind == "nucleation" else module)
        self._set_calculation_page(key)
        results_table = module.runtime.get("results_table")
        if isinstance(results_table, QTableWidget):
            _restore_cached_module_result(module, results_table)

    def _is_unavailable_calculation_module(self, module_kind: str) -> bool:
        """Return whether a calculation module is intentionally unavailable."""
        return module_kind == "nucleation"

    def _show_unavailable_calculation_module(self, module_kind: str) -> None:
        """Show a release warning for unfinished calculation modules."""
        if module_kind == "nucleation":
            self._show_future_feature_warning(
                "Nucleation Module",
                "The nucleation calculation module is not fully implemented "
                "in this release. Please use Equilibrium or Solidification "
                "modules for the current calculation workflow.",
            )

    def _set_header_unit_controls(
        self,
        module: CalculationModule | None,
    ) -> None:
        if self.calculation_unit_row is None:
            return
        self.calculation_unit_row.setVisible(module is not None)
        if module is None:
            return
        if (
            self.calculation_temperature_unit is None
            or self.calculation_pressure_unit is None
            or self.calculation_amount_unit is None
        ):
            return
        self._updating_calculation_units = True
        self.calculation_temperature_unit.setCurrentText(module.temperature_unit)
        self.calculation_pressure_unit.setCurrentText(module.pressure_unit)
        self.calculation_amount_unit.setCurrentText(module.amount_unit)
        self._updating_calculation_units = False

    def _store_current_module_units(self) -> None:
        if self._updating_calculation_units:
            return
        if (
            self.calculation_temperature_unit is None
            or self.calculation_pressure_unit is None
            or self.calculation_amount_unit is None
        ):
            return
        _session, module = self._current_calculation_context()
        if module is None:
            return
        module.temperature_unit = self.calculation_temperature_unit.currentText()
        module.pressure_unit = self.calculation_pressure_unit.currentText()
        module.amount_unit = self.calculation_amount_unit.currentText()
        if module.runtime:
            self._apply_calculation_type_to_runtime(module)

    def _open_session_unit_dialog(self) -> None:
        session = self._active_calculation_session()
        if session is None:
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Session Units")
        dialog.setModal(True)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(14)

        form = QFormLayout()
        temperature = _session_unit_combo_box(["K", "C", "F", "R"])
        temperature.setObjectName("SessionTemperatureUnit")
        pressure = _session_unit_combo_box(["atm", "psi", "bar", "Pa", "kPa"])
        pressure.setObjectName("SessionPressureUnit")
        amount = _session_unit_combo_box(_AMOUNT_UNIT_OPTIONS)
        amount.setObjectName("SessionAmountUnit")
        temperature.setCurrentText(session.temperature_unit)
        pressure.setCurrentText(session.pressure_unit)
        amount.setCurrentText(session.amount_unit)
        form.addRow("Temperature", temperature)
        form.addRow("Pressure", pressure)
        form.addRow("Amount", amount)
        layout.addLayout(form)

        actions = QHBoxLayout()
        actions.addStretch(1)
        apply = QPushButton("Apply")
        apply.setObjectName("PrimaryButton")
        cancel = QPushButton("Cancel")
        cancel.setObjectName("SmallActionButton")
        cancel.setDefault(True)
        cancel.setAutoDefault(True)
        apply.setAutoDefault(False)
        apply.clicked.connect(
            lambda _checked=False: (
                self._set_session_default_units(
                    session,
                    temperature.currentText(),
                    pressure.currentText(),
                    amount.currentText(),
                ),
                dialog.accept(),
            )
        )
        cancel.clicked.connect(dialog.reject)
        actions.addWidget(apply)
        actions.addWidget(cancel)
        layout.addLayout(actions)
        dialog.exec()

    def _set_session_default_units(
        self,
        session: CalculationSession,
        temperature_unit: str,
        pressure_unit: str,
        amount_unit: str,
    ) -> None:
        session.temperature_unit = temperature_unit
        session.pressure_unit = pressure_unit
        session.amount_unit = amount_unit
        for module in session.modules.values():
            module.temperature_unit = temperature_unit
            module.pressure_unit = pressure_unit
            module.amount_unit = amount_unit
        self._update_current_calculation_heading(session)

    def _open_session_phase_dialog(self) -> None:
        session = self._active_calculation_session()
        if session is None:
            return
        active_database = self._active_session_database(session)
        if active_database is None:
            QMessageBox.information(
                self,
                "Phase Selection",
                "Load and select a session database before choosing phases.",
            )
            return

        phase_names = _database_phase_names(active_database.database)
        if not phase_names:
            QMessageBox.warning(
                self,
                "Phase Selection",
                f"No phases were found in {active_database.name}.",
            )
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Session Phase Selection")
        dialog.setModal(True)
        dialog.resize(620, 520)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(12)

        title = QLabel(active_database.name)
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        phase_tree = QTreeWidget()
        phase_tree.setObjectName("DatabaseTree")
        phase_tree.setHeaderHidden(True)
        phase_tree.setColumnCount(3)
        phase_tree.header().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        _populate_phase_tree(
            phase_tree,
            phase_names,
            session.default_phase_names or None,
            fit_to_contents=False,
        )

        tool_row = QHBoxLayout()
        tool_row.setSpacing(8)
        all_solutions = QPushButton("All Solutions")
        all_solutions.setObjectName("SmallActionButton")
        all_solutions.clicked.connect(
            lambda _checked=False: _toggle_phase_group(phase_tree, "Solutions")
        )
        all_compounds = QPushButton("All Compounds")
        all_compounds.setObjectName("SmallActionButton")
        all_compounds.clicked.connect(
            lambda _checked=False: _toggle_phase_group(phase_tree, "Compounds")
        )
        tool_row.addWidget(all_solutions)
        tool_row.addWidget(all_compounds)
        tool_row.addStretch(1)
        layout.addLayout(tool_row)
        layout.addWidget(phase_tree, 1)

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
            selected = _selected_phase_names(phase_tree)
            if not selected:
                QMessageBox.warning(
                    dialog,
                    "Phase Selection",
                    "Select at least one phase.",
                )
                return
            self._set_session_default_phases(session, selected)
            dialog.accept()

        apply.clicked.connect(apply_selection)
        cancel.clicked.connect(dialog.reject)
        actions.addWidget(apply)
        actions.addWidget(cancel)
        layout.addLayout(actions)
        dialog.exec()

    def _set_session_default_phases(
        self,
        session: CalculationSession,
        phase_names: list[str],
    ) -> None:
        session.default_phase_names = list(phase_names)
        for module in session.modules.values():
            phase_tree = module.runtime.get("phase_tree")
            if isinstance(phase_tree, QTreeWidget) and module.phase_names:
                _set_phase_tree_selection(phase_tree, session.default_phase_names)
                self._refresh_nucleation_undercooling_runtime(module)
        self._update_current_calculation_heading(session)

    def _calculate_current_module(self) -> None:
        session, module = self._current_calculation_context()
        if session is None or module is None:
            QMessageBox.information(
                self,
                "Calculate",
                "Select a calculation module first.",
            )
            return
        if module.kind == "nucleation":
            self._show_unavailable_calculation_module(module.kind)
            return

        if self._calculation_in_progress:
            runtime = self._ensure_calculation_module_runtime(session, module)
            status = runtime["status"]
            self._enqueue_calculation(
                lambda session_id=session.id, module_id=module.id: (
                    self._calculate_module_by_id(session_id, module_id)
                ),
                queued_callback=lambda position: self._set_calculation_queued_status(
                    status,
                    position,
                    module.name,
                ),
                calculation_name=module.name,
                canceled_callback=lambda: self._set_calculation_queue_canceled_status(
                    status,
                    module.name,
                ),
            )
            return

        self._calculate_module(session, module)

    def _calculate_module_by_id(self, session_id: str, module_id: str) -> None:
        """Run a queued module calculation if it still exists."""
        session = self.calculation_sessions.get(session_id)
        module = session.modules.get(module_id) if session is not None else None
        if session is None or module is None:
            self._append_calculation_log_message(
                "Queued calculation skipped because the module no longer exists."
            )
            return
        if module.kind == "nucleation":
            self._show_unavailable_calculation_module(module.kind)
            return
        self._calculate_module(session, module)

    def _calculate_module(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> None:
        """Run one calculation module through the serialized solver lane."""
        runtime = self._ensure_calculation_module_runtime(session, module)
        self._run_calculation(
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
        )

    def _active_calculation_session(
        self,
        create_if_missing: bool = False,
    ) -> CalculationSession | None:
        item = None
        if self.calculation_tree is not None:
            selected = self.calculation_tree.selectedItems()
            if selected:
                item = selected[0]
        item = item or self.last_calculation_item

        if item is not None:
            session_id = self._session_id_for_item(item)
            if session_id in self.calculation_sessions:
                return self.calculation_sessions[session_id]

        if create_if_missing:
            return self._add_calculation_session()
        return None

    def _session_id_for_item(self, item: QTreeWidgetItem | None) -> str:
        while item is not None:
            value = item.data(0, Qt.ItemDataRole.UserRole)
            if isinstance(value, str):
                return value
            item = item.parent()
        return ""

    def _find_session_item(self, session_id: str) -> QTreeWidgetItem | None:
        if self.calculation_tree is None:
            return None
        for row in range(self.calculation_tree.topLevelItemCount()):
            item = self.calculation_tree.topLevelItem(row)
            if item.data(0, Qt.ItemDataRole.UserRole) == session_id:
                return item
        return None

    def _select_calculation_session_item(self, session_id: str) -> None:
        """Select a session tree item without re-entering tree navigation."""
        if self.calculation_tree is None:
            return
        item = self._find_session_item(session_id)
        if item is None:
            return
        previous_blocked = self.calculation_tree.blockSignals(True)
        self.calculation_tree.setCurrentItem(item)
        self.calculation_tree.blockSignals(previous_blocked)
        self.last_calculation_item = item

    def _add_session_module_item(
        self,
        session_id: str,
        module: CalculationModule,
    ) -> None:
        item = self._find_session_item(session_id)
        if item is None:
            return
        child = QTreeWidgetItem([module.name])
        child.setIcon(0, QIcon(str(calculation_icon_path(module.kind))))
        child.setData(0, Qt.ItemDataRole.UserRole, session_id)
        child.setData(0, Qt.ItemDataRole.UserRole + 1, module.id)
        child.setData(0, Qt.ItemDataRole.UserRole + 2, module.kind)
        child.setFlags(child.flags() | Qt.ItemFlag.ItemIsEditable)
        item.addChild(child)
        item.setExpanded(True)
        if self.calculation_tree is not None:
            self.calculation_tree.setCurrentItem(child)

    def _module_label(self, module_id: str) -> str:
        if module_id == "solidification":
            return "Solidification"
        if module_id == "nucleation":
            return "Nucleation"
        return "Equilibrium"

    def _remember_calculation_selection(self) -> None:
        if self.calculation_tree is None:
            return
        selected = self.calculation_tree.selectedItems()
        if selected:
            self.last_calculation_item = selected[0]
            item = selected[0]
            item_kind = item.data(0, Qt.ItemDataRole.UserRole + 2)
            module_id = item.data(0, Qt.ItemDataRole.UserRole + 1)
            session_id = self._session_id_for_item(item)
            if item_kind == "result_overview":
                session = self.calculation_sessions.get(session_id)
                if session is not None:
                    self._show_session_result_overview(session)
                return
            if isinstance(module_id, str):
                session = self.calculation_sessions.get(session_id)
                module = session.modules.get(module_id) if session else None
                if session is not None and module is not None:
                    self._show_calculation_module(session, module)
                return
            self._show_calculation_overview(self.calculation_sessions.get(session_id))

    def _rename_calculation_item(
        self,
        item: QTreeWidgetItem,
        column: int,
    ) -> None:
        if self._updating_calculation_tree or column != 0:
            return
        session_id = self._session_id_for_item(item)
        session = self.calculation_sessions.get(session_id)
        if session is None:
            return
        if item.data(0, Qt.ItemDataRole.UserRole + 2) == "result_overview":
            self._set_calculation_tree_text(item, "Overview")
            return

        new_name = item.text(0).strip()
        if not new_name:
            self._updating_calculation_tree = True
            module_id = item.data(0, Qt.ItemDataRole.UserRole + 1)
            item.setText(
                0,
                session.modules[module_id].name
                if isinstance(module_id, str) and module_id in session.modules
                else session.name,
            )
            self._updating_calculation_tree = False
            return

        module_id = item.data(0, Qt.ItemDataRole.UserRole + 1)
        if isinstance(module_id, str):
            module = session.modules.get(module_id)
            if module is None:
                return
            old_source_path = module.source_path
            module.name = self._unique_calculation_module_name(
                session,
                new_name,
                exclude_module_id=module.id,
            )
            if module.source_path == old_source_path:
                self._remember_stale_module_source_path(session, module)
            self._set_calculation_tree_text(item, module.name)
            self._update_current_calculation_heading(session)
        else:
            old_source_path = session.source_path
            session.name = self._unique_calculation_session_name(
                new_name,
                exclude_session_id=session.id,
            )
            if session.source_path == old_source_path:
                self._remember_stale_calculation_source_path(session.source_path)
            self._set_calculation_tree_text(item, session.name)
            current = (
                self.calculation_stack.currentWidget()
                if self.calculation_stack is not None
                else None
            )
            if (
                current is not None
                and current.property("session_id") == session.id
                and current.property("module_id") is None
            ):
                self._show_calculation_overview(session)
            else:
                self._update_current_calculation_heading(session)

    def _current_calculation_context(
        self,
    ) -> tuple[CalculationSession | None, CalculationModule | None]:
        current = (
            self.calculation_stack.currentWidget()
            if self.calculation_stack is not None
            else None
        )
        session_id = current.property("session_id") if current is not None else None
        module_id = current.property("module_id") if current is not None else None
        if not isinstance(session_id, str) and self.last_calculation_item is not None:
            session_id = self._session_id_for_item(self.last_calculation_item)
            module_id = self.last_calculation_item.data(0, Qt.ItemDataRole.UserRole + 1)
        if module_id == _RESULT_OVERVIEW_ITEM_ID:
            module_id = None
        session = (
            self.calculation_sessions.get(session_id)
            if isinstance(session_id, str)
            else None
        )
        module = (
            session.modules.get(module_id)
            if session is not None and isinstance(module_id, str)
            else None
        )
        return session, module

    def _rename_current_session_from_heading(self) -> None:
        if self.calculation_badge.isReadOnly():
            return
        session, module = self._current_calculation_context()
        if session is None:
            return
        new_name = self.calculation_badge.text().strip()
        if not new_name:
            self.calculation_badge.setText(session.name)
            return
        old_source_path = session.source_path
        session.name = self._unique_calculation_session_name(
            new_name,
            exclude_session_id=session.id,
        )
        if session.source_path == old_source_path:
            self._remember_stale_calculation_source_path(session.source_path)
        self.calculation_badge.setText(session.name)
        item = self._find_session_item(session.id)
        if item is not None:
            self._set_calculation_tree_text(item, session.name)
        current = (
            self.calculation_stack.currentWidget()
            if self.calculation_stack is not None
            else None
        )
        if (
            current is not None
            and current.property("session_id") == session.id
            and current.property("module_id") is None
        ):
            self._show_calculation_overview(session)
        else:
            self._update_calculation_heading(session, module.name if module else "")

    def _rename_current_module_from_heading(self) -> None:
        if self.calculation_title.isReadOnly():
            return
        session, module = self._current_calculation_context()
        if session is None or module is None:
            return
        new_name = self.calculation_title.text().strip()
        if not new_name:
            self.calculation_title.setText(module.name)
            return
        old_source_path = module.source_path
        module.name = self._unique_calculation_module_name(
            session,
            new_name,
            exclude_module_id=module.id,
        )
        if module.source_path == old_source_path:
            self._remember_stale_module_source_path(session, module)
        self.calculation_title.setText(module.name)
        item = self._find_module_item(session.id, module.id)
        if item is not None:
            self._set_calculation_tree_text(item, module.name)
        self._update_calculation_heading(session, module.name)

    def _find_module_item(
        self,
        session_id: str,
        module_id: str,
    ) -> QTreeWidgetItem | None:
        session_item = self._find_session_item(session_id)
        if session_item is None:
            return None
        for row in range(session_item.childCount()):
            item = session_item.child(row)
            if item.data(0, Qt.ItemDataRole.UserRole + 1) == module_id:
                return item
        return None

    def _set_calculation_tree_text(
        self,
        item: QTreeWidgetItem,
        text: str,
    ) -> None:
        self._updating_calculation_tree = True
        item.setText(0, text)
        self._updating_calculation_tree = False

    def _restore_calculation_selection(self) -> None:
        if self.calculation_tree is None or self.last_calculation_item is None:
            return
        self.calculation_tree.setCurrentItem(self.last_calculation_item)
        self.last_calculation_item.setSelected(True)
