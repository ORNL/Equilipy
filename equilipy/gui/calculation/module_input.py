"""Build and update calculation module input forms."""

from __future__ import annotations

# ruff: noqa: F401,F403,F405,I001,E501

import copy
import json
import os
import pickle
import shutil
from pathlib import Path
from typing import Any

from PySide6.QtCore import QModelIndex, QRect, QSize, Qt, QThread, QTimer
from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap, QTransform
from PySide6.QtWidgets import (
    QAbstractItemView,
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


# Glyphs for the Results detach toggle: open-in-window / dock-back.
_DETACH_GLYPH = "⧉"  # ⧉
_ATTACH_SOURCE_GLYPH = "⇱"  # mirrored to point top-right (Unicode has no NE)
_DETACH_TOOLTIP = "Move the Results panel into its own window"
_ATTACH_TOOLTIP = "Return the Results panel to the main window"

_GLYPH_ICON_CACHE: dict[tuple[str, bool], QIcon] = {}


def _glyph_icon(reference: QWidget, glyph: str, mirror: bool = False) -> QIcon:
    """Tight-cropped icon rendered from a text glyph.

    The glyph is drawn large, cropped to its alpha bounds, and optionally
    mirrored horizontally; scaling down through the button's icon size then
    fills the icon area with no built-in whitespace bezel.
    """
    key = (glyph, mirror)
    cached = _GLYPH_ICON_CACHE.get(key)
    if cached is not None:
        return cached
    canvas = 96
    pixmap = QPixmap(canvas, canvas)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    font = reference.font()
    font.setPointSizeF(60.0)
    painter.setFont(font)
    painter.setPen(QColor("#e7ebec"))
    painter.drawText(
        QRect(0, 0, canvas, canvas),
        Qt.AlignmentFlag.AlignCenter,
        glyph,
    )
    painter.end()

    image = pixmap.toImage()
    left, right, top, bottom = canvas, -1, canvas, -1
    for y in range(canvas):
        for x in range(canvas):
            if image.pixelColor(x, y).alpha() > 0:
                left = min(left, x)
                right = max(right, x)
                top = min(top, y)
                bottom = max(bottom, y)
    if right >= left and bottom >= top:
        pixmap = pixmap.copy(
            QRect(left, top, right - left + 1, bottom - top + 1)
        )
    if mirror:
        pixmap = pixmap.transformed(QTransform().scale(-1, 1))
    icon = QIcon(pixmap)
    _GLYPH_ICON_CACHE[key] = icon
    return icon


class _DetachedResultWindow(QWidget):
    """Top-level window hosting a detached module Results section."""

    def __init__(self, parent: QWidget, title: str) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle(title)
        self.setObjectName("DetachedResultsWindow")
        self.resize(940, 660)
        self._section: QWidget | None = None
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(18, 16, 18, 16)
        self.placeholder: QWidget | None = None
        self.page_layout = None
        self.detach_button: QPushButton | None = None
        self.reattach_callback = None
        self.force_closing = False

    def host(self, section: QWidget) -> None:
        """Adopt the Results section into this window."""
        self._section = section
        self._layout.addWidget(section)

    def take_section(self) -> QWidget | None:
        """Release the hosted Results section for reattachment."""
        section = self._section
        if section is not None:
            self._layout.removeWidget(section)
            section.setParent(None)
        self._section = None
        return section

    def closeEvent(self, event) -> None:
        """Reattach the Results section when the user closes the window."""
        if not self.force_closing and self.reattach_callback is not None:
            callback = self.reattach_callback
            self.reattach_callback = None
            callback()
        super().closeEvent(event)


class CalculationModuleInputMixin:
    """Module editor form and option dialog controller methods."""

    def _calculation_module_page(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> QWidget:
        if module.kind == "nucleation":
            return self._nucleation_placeholder_page(module)

        page = QScrollArea()
        page.setObjectName("EditorScroll")
        page.setWidgetResizable(True)
        page.setFrameShape(QFrame.Shape.NoFrame)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(14)

        composition_table = _composition_table()
        composition_table.setProperty("amount_unit", module.amount_unit)
        temperature = QLineEdit("2000" if module.kind == "solidification" else "1600")
        composition_table.set_boundary_focus_widgets(next_widget=temperature)
        pressure = QLineEdit("1")
        delta_t = QLineEdit(_default_solidification_delta_t(module))
        liquid_phase = QLineEdit("LIQUID")
        from_liquidus = QCheckBox("From Liquidus")
        from_liquidus.setChecked(module.start_from_liquidus)
        from_liquidus.toggled.connect(
            lambda checked: setattr(module, "start_from_liquidus", checked)
        )
        status = self._module_status_widget()
        status.setPlainText("Ready.")
        def sync_condition_state() -> None:
            self._validate_composition_table(
                session,
                composition_table,
                status,
                quiet=True,
            )
            self._update_manual_batch_point_count(module)
            self._refresh_nucleation_undercooling_runtime(module)

        composition_table.itemChanged.connect(lambda _item: sync_condition_state())
        liquid_phase.textChanged.connect(
            lambda _text: self._refresh_nucleation_undercooling_runtime(module)
        )

        add_row = QPushButton("Add Species")
        add_row.setObjectName("SmallActionButton")
        add_row.clicked.connect(
            lambda: self._add_composition_table_row(
                session,
                composition_table,
                status,
            )
        )
        remove_row = QPushButton("Remove Species")
        remove_row.setObjectName("SmallActionButton")
        remove_row.clicked.connect(
            lambda: self._remove_composition_table_rows(
                session,
                composition_table,
                status,
            )
        )
        phase_selection = QPushButton("Phase Selection")
        phase_selection.setObjectName("SmallActionButton")
        type_button = QPushButton("Type")
        type_button.setObjectName("PrimaryButton")
        type_button.setVisible(module.kind in {"equilibrium", "solidification"})
        condition_section, condition_layout = _collapsible_form_section(
            "Set Condition",
            [add_row, remove_row, phase_selection, type_button],
        )
        condition_layout.addWidget(composition_table)
        batch_conditions_table = QTableWidget(0, 0)
        batch_conditions_table.setObjectName("FormTable")
        batch_conditions_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        batch_conditions_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        batch_conditions_table.verticalHeader().setVisible(False)
        batch_conditions_table.horizontalHeader().setStretchLastSection(True)
        _disable_inner_scrollbars(batch_conditions_table)
        batch_conditions_table.setVisible(False)
        condition_layout.addWidget(batch_conditions_table)
        batch_point_count = QLabel("")
        batch_point_count.setObjectName("DatabaseListLabel")
        batch_point_count.setVisible(False)
        condition_layout.addWidget(batch_point_count)

        def field_column(label_text: str, editor: QWidget) -> QWidget:
            container = QWidget()
            column_layout = QVBoxLayout(container)
            column_layout.setContentsMargins(0, 0, 0, 0)
            column_layout.setSpacing(4)
            label = QLabel(label_text)
            column_layout.addWidget(label)
            column_layout.addWidget(editor)
            return container

        condition_grid = QGridLayout()
        condition_grid.setHorizontalSpacing(12)
        condition_grid.setVerticalSpacing(10)
        temperature_label = QLabel("Temperature:")
        pressure_label = QLabel("Pressure:")
        transition_search = QCheckBox("Transitions")
        transition_search.setChecked(module.transition_search)
        transition_search.setVisible(
            module.kind == "equilibrium"
            and _temperature_text_is_range(temperature.text())
        )
        transition_search.toggled.connect(
            lambda checked: setattr(module, "transition_search", checked)
        )
        temperature.textChanged.connect(
            lambda text: transition_search.setVisible(
                module.calculation_type != "batch"
                and module.kind == "equilibrium"
                and _temperature_text_is_range(text)
            )
        )
        temperature.textChanged.connect(lambda _text: sync_condition_state())
        pressure.textChanged.connect(lambda _text: sync_condition_state())
        if module.kind == "solidification":
            temperature_field = field_column("Temperature:", temperature)
            delta_t_field = field_column("Delta T:", delta_t)
            pressure_field = field_column("Pressure:", pressure)
            liquid_phase_field = field_column("Liquid Phase:", liquid_phase)
            from_liquidus_container = QWidget()
            from_liquidus_layout = QVBoxLayout(from_liquidus_container)
            from_liquidus_layout.setContentsMargins(0, 0, 0, 0)
            from_liquidus_layout.setSpacing(4)
            from_liquidus_label_spacer = QLabel("")
            from_liquidus_label_spacer.setFixedHeight(
                QLabel("Temperature:").sizeHint().height()
            )
            from_liquidus_row = QWidget()
            from_liquidus_row.setMinimumHeight(temperature.sizeHint().height())
            from_liquidus_row_layout = QHBoxLayout(from_liquidus_row)
            from_liquidus_row_layout.setContentsMargins(0, 0, 0, 0)
            from_liquidus_row_layout.addWidget(
                from_liquidus,
                0,
                Qt.AlignmentFlag.AlignVCenter,
            )
            from_liquidus_row_layout.addStretch(1)
            from_liquidus_layout.addWidget(from_liquidus_label_spacer)
            from_liquidus_layout.addWidget(from_liquidus_row)
            condition_grid.addWidget(temperature_field, 0, 0)
            condition_grid.addWidget(delta_t_field, 0, 1)
            condition_grid.addWidget(from_liquidus_container, 0, 2)
            condition_grid.addWidget(pressure_field, 1, 0)
            condition_grid.addWidget(liquid_phase_field, 1, 1)
            condition_grid.setColumnStretch(0, 1)
            condition_grid.setColumnStretch(1, 1)
            condition_grid.setColumnStretch(2, 1)
        else:
            condition_grid.addWidget(temperature_label, 0, 0)
            condition_grid.addWidget(temperature, 0, 1)
            condition_grid.addWidget(pressure_label, 1, 0)
            condition_grid.addWidget(pressure, 1, 1)
            condition_grid.addWidget(transition_search, 0, 2, 1, 2)
            condition_grid.setColumnStretch(1, 1)
            condition_grid.setColumnStretch(3, 1)
        condition_layout.addLayout(condition_grid)
        layout.addWidget(condition_section)

        nucleation_undercooling_tree = QTreeWidget(content)
        nucleation_undercooling_tree.setObjectName("DatabaseTree")
        nucleation_undercooling_tree.setColumnCount(2)
        nucleation_undercooling_tree.setHeaderLabels(
            ["Phase", f"Undercooling [{module.temperature_unit}]"]
        )
        nucleation_undercooling_tree.header().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        _disable_inner_scrollbars(nucleation_undercooling_tree)
        (
            nucleation_undercooling_section,
            nucleation_undercooling_layout,
        ) = _collapsible_form_section("Nucleation Undercooling")
        nucleation_undercooling_layout.addWidget(
            nucleation_undercooling_tree
        )
        layout.addWidget(nucleation_undercooling_section)
        nucleation_undercooling_tree.itemChanged.connect(
            lambda _item, _column: _store_nucleation_undercooling_tree(
                module,
                nucleation_undercooling_tree,
            )
        )

        phase_tree = QTreeWidget(content)
        phase_tree.setObjectName("DatabaseTree")
        phase_tree.setHeaderHidden(True)
        phase_tree.setColumnCount(3)
        phase_tree.header().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        _disable_inner_scrollbars(phase_tree)
        phase_tree.setVisible(False)

        results_table = QTableWidget(0, 0)
        results_table.setObjectName("FormTable")
        results_table.setVisible(False)
        results_table.verticalHeader().setVisible(False)
        results_table.horizontalHeader().setStretchLastSection(True)
        results_table.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        results_table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        result_table_button = QPushButton("Table")
        result_table_button.setObjectName("ResultViewButton")
        result_table_button.setCheckable(True)
        result_table_button.setChecked(True)
        result_figure_button = QPushButton("Figure")
        result_figure_button.setObjectName("ResultViewButton")
        result_figure_button.setCheckable(True)
        result_figure_button.setEnabled(False)
        result_columns_button = QPushButton("Columns")
        result_columns_button.setObjectName("SmallActionButton")
        result_columns_button.clicked.connect(
            lambda: self._open_result_columns_dialog(module, results_table)
        )
        result_show_all_button = QPushButton("Show All")
        result_show_all_button.setObjectName("SmallActionButton")
        result_show_all_button.clicked.connect(
            lambda: self._show_all_module_results(module, results_table)
        )
        result_export_button = QPushButton("Export")
        result_export_button.setObjectName("PrimaryButton")
        result_export_button.clicked.connect(
            lambda: self._export_module_results(module, results_table)
        )
        result_detach_button = QPushButton("")
        result_detach_button.setObjectName("SmallActionButton")
        result_detach_button.setToolTip(_DETACH_TOOLTIP)
        # Fixed square: the shared SmallActionButton style pads 12px per
        # side, which clips content inside the narrow button; keep only the
        # vertical padding and pin the size so the icon toggle is stable.
        result_detach_button.setFixedSize(34, 33)
        result_detach_button.setStyleSheet("padding: 5px 0px;")
        result_detach_button.setIcon(
            _glyph_icon(result_detach_button, _DETACH_GLYPH)
        )
        result_detach_button.setIconSize(QSize(17, 17))
        result_view_stack = QStackedWidget()
        result_view_stack.setMinimumHeight(552)
        result_view_stack.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        result_table_page = QWidget()
        result_table_layout = QVBoxLayout(result_table_page)
        result_table_layout.setContentsMargins(0, 0, 0, 0)
        result_table_layout.setSpacing(8)
        result_table_controls = QWidget()
        result_table_controls.setObjectName("ResultTableControls")
        result_table_controls_layout = QHBoxLayout(result_table_controls)
        result_table_controls_layout.setContentsMargins(0, 0, 0, 0)
        result_table_controls_layout.setSpacing(8)
        result_table_controls_layout.addWidget(result_columns_button)
        result_table_controls_layout.addWidget(result_show_all_button)
        result_table_controls_layout.addStretch(1)
        result_table_layout.addWidget(result_table_controls)
        result_table_layout.addWidget(results_table)
        result_table_layout.addStretch(1)

        result_figure_page = QWidget()
        result_figure_layout = QVBoxLayout(result_figure_page)
        result_figure_layout.setContentsMargins(0, 0, 0, 0)
        result_figure_layout.setSpacing(8)
        result_figure_controls = QWidget()
        result_figure_controls.setObjectName("ResultFigureControls")
        result_figure_controls.setMaximumHeight(40)
        result_figure_controls_layout = QHBoxLayout(result_figure_controls)
        result_figure_controls_layout.setContentsMargins(0, 0, 0, 0)
        result_figure_controls_layout.setSpacing(8)
        result_figure_x_axis = QComboBox()
        result_figure_x_axis.setObjectName("ResultFigureXAxis")
        result_figure_y_axis = QComboBox()
        result_figure_y_axis.setObjectName("ResultFigureYAxis")
        result_figure_y_columns_button = QPushButton("select")
        result_figure_y_columns_button.setObjectName("SmallActionButton")
        result_figure_palette = QComboBox()
        result_figure_palette.setObjectName("ResultFigurePalette")
        result_figure_qdot_label = QLabel("q&#775; [J/s]")
        result_figure_qdot_label.setTextFormat(Qt.TextFormat.RichText)
        result_figure_qdot = QLineEdit("1")
        result_figure_qdot.setObjectName("ResultFigureQdot")
        result_figure_qdot.setMaximumWidth(96)
        for combo, width in (
            (result_figure_x_axis, 132),
            (result_figure_y_axis, 150),
            (result_figure_palette, 140),
        ):
            combo.setMinimumContentsLength(4)
            combo.setSizeAdjustPolicy(
                QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
            )
            combo.setMinimumWidth(width)
            combo.setMaximumWidth(width)
            combo.setSizePolicy(
                QSizePolicy.Policy.Fixed,
                QSizePolicy.Policy.Fixed,
            )
        result_figure_y_columns_button.setMinimumWidth(150)
        result_figure_y_columns_button.setMaximumWidth(150)
        result_figure_y_columns_button.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Fixed,
        )
        result_figure_y_columns_button.setFont(result_figure_x_axis.font())
        result_figure_y_columns_button.setVisible(False)

        result_figure_qdot_field = QWidget()
        result_figure_qdot_layout = QHBoxLayout(result_figure_qdot_field)
        result_figure_qdot_layout.setContentsMargins(0, 0, 0, 0)
        result_figure_qdot_layout.setSpacing(8)
        result_figure_qdot_layout.addWidget(result_figure_qdot_label)
        result_figure_qdot_layout.addWidget(result_figure_qdot)
        result_figure_qdot_field.setVisible(False)
        result_figure_controls_layout.addWidget(QLabel("x-axis:"))
        result_figure_controls_layout.addWidget(result_figure_x_axis)
        result_figure_controls_layout.addWidget(QLabel("y-axis:"))
        result_figure_controls_layout.addWidget(result_figure_y_axis)
        result_figure_controls_layout.addWidget(result_figure_y_columns_button)
        result_figure_controls_layout.addWidget(QLabel("Palette:"))
        result_figure_controls_layout.addWidget(result_figure_palette)
        result_figure_controls_layout.addWidget(result_figure_qdot_field)
        result_figure_controls_layout.addStretch(1)
        result_figure_status = QLabel("")
        result_figure_status.setObjectName("ResultFigureStatus")
        result_figure_status.setVisible(False)
        result_figure_chart = _new_scheil_figure_chart_view()
        result_figure_legend = _new_result_figure_legend_widget()
        result_figure_canvas = QWidget()
        result_figure_canvas.setObjectName("ResultFigureCanvas")
        result_figure_canvas_layout = QVBoxLayout(result_figure_canvas)
        result_figure_canvas_layout.setContentsMargins(0, 0, 0, 0)
        result_figure_canvas_layout.setSpacing(6)
        result_figure_canvas_layout.addWidget(result_figure_chart, 1)
        result_figure_canvas_layout.addWidget(result_figure_legend)
        result_figure_layout.addWidget(result_figure_controls)
        result_figure_layout.addWidget(result_figure_status)
        result_figure_layout.addWidget(result_figure_canvas, 1)
        result_view_stack.addWidget(result_table_page)
        result_view_stack.addWidget(result_figure_page)

        def show_result_table() -> None:
            result_view_stack.setCurrentIndex(0)
            _refresh_module_result_view(module)

        def show_result_figure() -> None:
            if not result_figure_button.isEnabled():
                return
            result_view_stack.setCurrentIndex(1)
            _refresh_module_result_view(module)

        result_table_button.clicked.connect(show_result_table)
        result_figure_button.clicked.connect(show_result_figure)
        result_figure_x_axis.currentTextChanged.connect(
            lambda _text: _refresh_module_result_view(module)
        )
        result_figure_y_axis.currentTextChanged.connect(
            lambda _text: _refresh_module_result_view(module)
        )
        result_figure_y_columns_button.clicked.connect(
            lambda: _open_result_figure_y_columns_dialog(module)
        )
        result_figure_palette.currentTextChanged.connect(
            lambda _text: _populate_result_figure(module)
            if result_view_stack.currentIndex() == 1
            else None
        )
        result_figure_qdot.textChanged.connect(
            lambda _text: _populate_result_figure(module)
            if result_view_stack.currentIndex() == 1
            else None
        )

        phase_selection.clicked.connect(
            lambda: self._open_module_phase_dialog(
                session.id,
                module.id,
                composition_table,
                phase_tree,
                status,
            )
        )
        type_button.clicked.connect(
            lambda: self._open_calculation_type_dialog(module)
        )
        module.runtime = {
            "session_id": session.id,
            "composition_table": composition_table,
            "batch_conditions_table": batch_conditions_table,
            "batch_point_count": batch_point_count,
            "phase_tree": phase_tree,
            "temperature": temperature,
            "temperature_label": temperature_label,
            "pressure": pressure,
            "pressure_label": pressure_label,
            "transition_search": transition_search,
            "delta_t": delta_t,
            "liquid_phase": liquid_phase,
            "from_liquidus": from_liquidus,
            "nucleation_undercooling_section": nucleation_undercooling_section,
            "nucleation_undercooling_tree": nucleation_undercooling_tree,
            "results_table": results_table,
            "result_table_button": result_table_button,
            "result_figure_button": result_figure_button,
            "result_view_stack": result_view_stack,
            "result_figure_x_axis": result_figure_x_axis,
            "result_figure_y_axis": result_figure_y_axis,
            "result_figure_y_columns_button": result_figure_y_columns_button,
            "result_figure_palette": result_figure_palette,
            "result_figure_qdot_label": result_figure_qdot_label,
            "result_figure_qdot": result_figure_qdot,
            "result_figure_qdot_field": result_figure_qdot_field,
            "result_figure_controls": result_figure_controls,
            "result_figure_status": result_figure_status,
            "result_figure_chart": result_figure_chart,
            "result_figure_legend": result_figure_legend,
            "result_figure_canvas": result_figure_canvas,
            "result_columns_button": result_columns_button,
            "result_show_all_button": result_show_all_button,
            "result_export_button": result_export_button,
            "status": status,
            "single_condition_widgets": [
                composition_table,
                temperature_label,
                temperature,
                pressure_label,
                pressure,
            ],
            "single_condition_buttons": [add_row, remove_row],
        }
        if module.kind == "solidification":
            module.runtime["single_condition_widgets"] = [
                composition_table,
                temperature_field,
                pressure_field,
            ]
            module.runtime["batch_visible_widgets"] = [
                delta_t_field,
                liquid_phase_field,
                from_liquidus_container,
            ]

        results_section, results_layout = _collapsible_form_section(
            "Results",
            [
                result_table_button,
                result_figure_button,
                result_export_button,
                result_detach_button,
            ],
        )
        results_section.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        results_layout.addWidget(result_view_stack)
        layout.addWidget(results_section, 1)

        def toggle_results_detached() -> None:
            if module.id in self._detached_result_windows:
                self._reattach_module_results(module.id)
            else:
                self._detach_module_results(
                    module,
                    results_section,
                    layout,
                    result_detach_button,
                )

        result_detach_button.clicked.connect(
            lambda _checked=False: toggle_results_detached()
        )

        page.setWidget(content)
        _restore_cached_module_result(module, results_table)
        _refresh_module_result_view(module)
        self._apply_calculation_type_to_runtime(module)
        self._validate_composition_table(session, composition_table, status, quiet=True)
        return page

    def _detach_module_results(
        self,
        module,
        results_section: QWidget,
        page_layout,
        detach_button: QPushButton,
    ) -> None:
        """Move a module's Results section into its own top-level window."""
        index = page_layout.indexOf(results_section)
        if index < 0:
            return
        placeholder = QFrame()
        placeholder.setObjectName("DetachedResultsPlaceholder")
        placeholder_layout = QVBoxLayout(placeholder)
        placeholder_layout.setContentsMargins(0, 24, 0, 24)
        note = QLabel("Results are open in a separate window.")
        note.setObjectName("EditorSubtitle")
        note.setAlignment(Qt.AlignmentFlag.AlignCenter)
        placeholder_layout.addWidget(note)

        page_layout.removeWidget(results_section)
        page_layout.insertWidget(index, placeholder, 1)

        window = _DetachedResultWindow(self, f"{module.name} - Results")
        window.host(results_section)
        window.placeholder = placeholder
        window.page_layout = page_layout
        window.detach_button = detach_button
        window.reattach_callback = (
            lambda module_key=module.id: self._reattach_module_results(module_key)
        )
        self._detached_result_windows[module.id] = window
        detach_button.setIcon(
            _glyph_icon(detach_button, _ATTACH_SOURCE_GLYPH, mirror=True)
        )
        detach_button.setToolTip(_ATTACH_TOOLTIP)
        window.show()
        window.raise_()

    def _reattach_module_results(self, module_key: str) -> None:
        """Return a detached Results section to its module page."""
        window = self._detached_result_windows.pop(module_key, None)
        if window is None:
            return
        window.force_closing = True
        section = window.take_section()
        placeholder = window.placeholder
        page_layout = window.page_layout
        if placeholder is not None and page_layout is not None:
            index = page_layout.indexOf(placeholder)
            page_layout.removeWidget(placeholder)
            placeholder.deleteLater()
            if section is not None and index >= 0:
                page_layout.insertWidget(index, section, 1)
                section.show()
        if window.detach_button is not None:
            window.detach_button.setIcon(
                _glyph_icon(window.detach_button, _DETACH_GLYPH)
            )
            window.detach_button.setToolTip(_DETACH_TOOLTIP)
        window.close()
        window.deleteLater()

    def _close_detached_result_window(self, module_key: str) -> None:
        """Discard a detached Results window during module page teardown.

        The page owning the section is being deleted, so the section is
        destroyed together with the window instead of being reattached.
        """
        window = self._detached_result_windows.pop(module_key, None)
        if window is None:
            return
        window.force_closing = True
        window.close()
        window.deleteLater()

    def _module_status_widget(self) -> QPlainTextEdit:
        if self.calculation_log is None:
            fallback = QPlainTextEdit()
            fallback.setReadOnly(True)
            return fallback
        return self.calculation_log

    def _open_solidification_type_dialog(self, module: CalculationModule) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle(f"{module.name} Calculation Type")
        dialog.setModal(True)
        dialog.resize(840, 240)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(12)

        title = QLabel(module.name)
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        imported_condition = {
            key: list(values) for key, values in module.batch_condition.items()
        }
        imported_path = module.batch_condition_path
        current_model = _solidification_model_value(module.solidification_model)
        current_mode = "batch" if module.calculation_type == "batch" else "single"
        option_group = QButtonGroup(dialog)
        option_group.setExclusive(True)
        option_buttons: dict[tuple[str, str], QRadioButton] = {}

        import_button = QPushButton("Import Conditions")
        import_button.setObjectName("SmallActionButton")
        import_button.setIcon(QIcon(str(calculation_icon_path("upload"))))
        remove_button = QToolButton()
        remove_button.setObjectName("IconToolButton")
        remove_button.setIcon(
            self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon)
        )
        remove_button.setToolTip("Remove imported conditions")
        cpu_label = QLabel("CPU cores")
        cpu_label.setObjectName("DialogRowLabel")
        cpu_combo = QComboBox()
        _populate_batch_cpu_combo(cpu_combo, module.batch_cpu_count)

        def checked_option() -> QRadioButton | None:
            return option_group.checkedButton()

        def checked_mode() -> str:
            option = checked_option()
            if option is None:
                return current_mode
            return str(option.property("mode") or current_mode)

        def checked_model() -> str:
            option = checked_option()
            if option is None:
                return current_model
            return _solidification_model_value(str(option.property("model")))

        def sync_import_buttons() -> None:
            is_batch = checked_mode() == "batch"
            import_button.setEnabled(is_batch)
            remove_button.setEnabled(is_batch and bool(imported_condition))
            cpu_label.setEnabled(is_batch)
            cpu_combo.setEnabled(is_batch)

        def select_option(mode: str, model: str) -> None:
            option = option_buttons.get((mode, _solidification_model_value(model)))
            if option is not None:
                option.setChecked(True)
            sync_import_buttons()

        def add_row(label_text: str, mode: str) -> QHBoxLayout:
            row = QHBoxLayout()
            row.setSpacing(10)
            label = QLabel(label_text)
            label.setObjectName("DialogRowLabel")
            label.setMinimumWidth(130)
            row.addWidget(label)
            for value, option_label in _solidification_model_options():
                option = QRadioButton(option_label)
                option.setObjectName("DialogRadioOption")
                option.setAutoExclusive(True)
                option.setProperty("mode", mode)
                option.setProperty("model", value)
                option.setChecked(mode == current_mode and value == current_model)
                option.clicked.connect(lambda _checked=False: sync_import_buttons())
                option_group.addButton(option)
                option_buttons[(mode, value)] = option
                row.addWidget(option)
            row.addStretch(1)
            return row

        layout.addLayout(add_row("Single condition:", "single"))

        batch_row = add_row("Batch condition:", "batch")
        batch_row.addWidget(import_button)
        batch_row.addWidget(remove_button)
        batch_row.addSpacing(8)
        batch_row.addWidget(cpu_label)
        batch_row.addWidget(cpu_combo)
        layout.addLayout(batch_row)

        import_status = QLabel(
            Path(imported_path).name if imported_path else "No conditions imported"
        )
        import_status.setObjectName("DatabaseListLabel")
        layout.addWidget(import_status)
        sync_import_buttons()

        def import_conditions() -> None:
            nonlocal imported_condition, imported_path
            path, _selected_filter = QFileDialog.getOpenFileName(
                dialog,
                "Import Conditions",
                "",
                "CSV files (*.csv);;All files (*)",
            )
            if not path:
                return
            try:
                imported_condition = _read_condition_csv(path)
            except ValueError as exc:
                QMessageBox.warning(dialog, "Import Conditions", str(exc))
                return
            imported_path = path
            select_option("batch", checked_model())
            import_status.setText(Path(path).name)

        def remove_conditions() -> None:
            nonlocal imported_condition, imported_path
            imported_condition = {}
            imported_path = ""
            sync_import_buttons()
            import_status.setText("No conditions imported")

        import_button.clicked.connect(import_conditions)
        remove_button.clicked.connect(remove_conditions)

        actions = QHBoxLayout()
        actions.addStretch(1)
        apply = QPushButton("Apply")
        apply.setObjectName("PrimaryButton")
        cancel = QPushButton("Cancel")
        cancel.setObjectName("SmallActionButton")
        cancel.setDefault(True)
        cancel.setAutoDefault(True)
        apply.setAutoDefault(False)

        def apply_type() -> None:
            previous_signature = (
                module.calculation_type,
                module.solidification_model,
                module.batch_condition_path,
            )
            previous_model = _solidification_model_value(module.solidification_model)
            option = checked_option()
            if option is not None:
                module.solidification_model = _solidification_model_value(
                    str(option.property("model"))
                )
            is_batch = checked_mode() == "batch"

            if is_batch:
                module.calculation_type = "batch"
                module.batch_condition = {
                    key: list(values) for key, values in imported_condition.items()
                }
                module.batch_condition_path = imported_path
                module.batch_cpu_count = _batch_cpu_combo_value(cpu_combo)
            else:
                module.calculation_type = "single"
                module.batch_condition = {}
                module.batch_condition_path = ""
                module.batch_cpu_count = _batch_cpu_combo_value(cpu_combo)

            current_signature = (
                module.calculation_type,
                module.solidification_model,
                module.batch_condition_path,
            )
            if current_signature != previous_signature:
                module.phase_names = []
                phase_tree = module.runtime.get("phase_tree")
                if isinstance(phase_tree, QTreeWidget):
                    phase_tree.clear()
            if (
                module.solidification_model == "nucleoscheil"
                and previous_model != "nucleoscheil"
            ):
                delta_t = module.runtime.get("delta_t")
                if isinstance(delta_t, QLineEdit):
                    delta_t.setText("0.1")
            self._apply_calculation_type_to_runtime(module)
            dialog.accept()

        apply.clicked.connect(apply_type)
        cancel.clicked.connect(dialog.reject)
        actions.addWidget(apply)
        actions.addWidget(cancel)
        layout.addLayout(actions)
        dialog.exec()

    def _open_calculation_type_dialog(self, module: CalculationModule) -> None:
        if module.kind == "solidification":
            self._open_solidification_type_dialog(module)
            return
        if module.kind != "equilibrium":
            return

        dialog = QDialog(self)
        dialog.setWindowTitle(f"{module.name} Calculation Type")
        dialog.setModal(True)
        dialog.resize(680 if module.kind == "solidification" else 640, 240)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(18, 16, 18, 16)
        layout.setSpacing(12)

        title = QLabel(module.name)
        title.setObjectName("SectionTitle")
        layout.addWidget(title)

        single = QRadioButton("Single condition")
        batch = QRadioButton("Batch condition")
        mode_group = QButtonGroup(dialog)
        mode_group.addButton(single)
        mode_group.addButton(batch)
        single.setChecked(module.calculation_type != "batch")
        batch.setChecked(module.calculation_type == "batch")

        imported_condition = {
            key: list(values) for key, values in module.batch_condition.items()
        }
        imported_path = module.batch_condition_path
        current_model = _solidification_model_value(module.solidification_model)

        single_row = QHBoxLayout()
        single_row.addWidget(single)
        single_model_group = QButtonGroup(dialog)
        batch_model_group = QButtonGroup(dialog)
        single_model_buttons: dict[str, QRadioButton] = {}
        batch_model_buttons: dict[str, QRadioButton] = {}

        def add_model_options(
            row: QHBoxLayout,
            owner: QRadioButton,
            group: QButtonGroup,
            buttons: dict[str, QRadioButton],
        ) -> None:
            if module.kind != "solidification":
                return
            for value, label in _solidification_model_options():
                option = QRadioButton(label)
                option.setAutoExclusive(False)
                option.setProperty("model", value)
                option.setChecked(value == current_model)
                option.clicked.connect(
                    lambda _checked=False, owner=owner: owner.setChecked(True)
                )
                group.addButton(option)
                buttons[value] = option
                row.addWidget(option)

        add_model_options(single_row, single, single_model_group, single_model_buttons)
        single_row.addStretch(1)
        layout.addLayout(single_row)

        batch_row = QHBoxLayout()
        batch_row.setSpacing(8)
        batch_row.addWidget(batch)
        add_model_options(batch_row, batch, batch_model_group, batch_model_buttons)
        batch_row.addStretch(1)

        import_button = QPushButton("Import Conditions")
        import_button.setObjectName("SmallActionButton")
        import_button.setIcon(QIcon(str(calculation_icon_path("upload"))))
        remove_button = QToolButton()
        remove_button.setObjectName("IconToolButton")
        remove_button.setIcon(
            self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon)
        )
        remove_button.setToolTip("Remove imported conditions")
        cpu_label = QLabel("CPU cores")
        cpu_label.setObjectName("DialogRowLabel")
        cpu_combo = QComboBox()
        _populate_batch_cpu_combo(cpu_combo, module.batch_cpu_count)
        import_button.setEnabled(batch.isChecked())
        remove_button.setEnabled(batch.isChecked() and bool(imported_condition))
        cpu_label.setEnabled(batch.isChecked())
        cpu_combo.setEnabled(batch.isChecked())
        batch_row.addWidget(import_button)
        batch_row.addWidget(remove_button)
        batch_row.addSpacing(8)
        batch_row.addWidget(cpu_label)
        batch_row.addWidget(cpu_combo)
        layout.addLayout(batch_row)

        import_status = QLabel(
            Path(imported_path).name if imported_path else "No conditions imported"
        )
        import_status.setObjectName("DatabaseListLabel")
        layout.addWidget(import_status)

        def import_conditions() -> None:
            nonlocal imported_condition, imported_path
            path, _selected_filter = QFileDialog.getOpenFileName(
                dialog,
                "Import Conditions",
                "",
                "CSV files (*.csv);;All files (*)",
            )
            if not path:
                return
            try:
                imported_condition = _read_condition_csv(path)
            except ValueError as exc:
                QMessageBox.warning(dialog, "Import Conditions", str(exc))
                return
            imported_path = path
            batch.setChecked(True)
            remove_button.setEnabled(True)
            import_status.setText(Path(path).name)

        def remove_conditions() -> None:
            nonlocal imported_condition, imported_path
            imported_condition = {}
            imported_path = ""
            remove_button.setEnabled(False)
            import_status.setText("No conditions imported")

        import_button.clicked.connect(import_conditions)
        remove_button.clicked.connect(remove_conditions)
        single.toggled.connect(lambda checked: import_button.setEnabled(not checked))
        single.toggled.connect(
            lambda checked: remove_button.setEnabled(
                (not checked) and bool(imported_condition)
            )
        )
        batch.toggled.connect(lambda checked: import_button.setEnabled(checked))
        batch.toggled.connect(
            lambda checked: remove_button.setEnabled(
                checked and bool(imported_condition)
            )
        )
        batch.toggled.connect(lambda checked: cpu_label.setEnabled(checked))
        batch.toggled.connect(lambda checked: cpu_combo.setEnabled(checked))

        actions = QHBoxLayout()
        actions.addStretch(1)
        apply = QPushButton("Apply")
        apply.setObjectName("PrimaryButton")
        cancel = QPushButton("Cancel")
        cancel.setObjectName("SmallActionButton")
        cancel.setDefault(True)
        cancel.setAutoDefault(True)
        apply.setAutoDefault(False)

        def apply_type() -> None:
            previous_signature = (
                module.calculation_type,
                module.solidification_model,
                module.batch_condition_path,
            )
            if batch.isChecked():
                module.calculation_type = "batch"
                module.batch_condition = {
                    key: list(values) for key, values in imported_condition.items()
                }
                module.batch_condition_path = imported_path
                module.batch_cpu_count = _batch_cpu_combo_value(cpu_combo)
            else:
                module.calculation_type = "single"
                module.batch_condition = {}
                module.batch_condition_path = ""
                module.batch_cpu_count = _batch_cpu_combo_value(cpu_combo)
            if module.kind == "solidification":
                group = batch_model_group if batch.isChecked() else single_model_group
                checked_model = group.checkedButton()
                if checked_model is not None:
                    module.solidification_model = _solidification_model_value(
                        str(checked_model.property("model"))
                    )

            current_signature = (
                module.calculation_type,
                module.solidification_model,
                module.batch_condition_path,
            )
            if current_signature != previous_signature:
                module.phase_names = []
                phase_tree = module.runtime.get("phase_tree")
                if isinstance(phase_tree, QTreeWidget):
                    phase_tree.clear()
            self._apply_calculation_type_to_runtime(module)
            dialog.accept()

        apply.clicked.connect(apply_type)
        cancel.clicked.connect(dialog.reject)
        actions.addWidget(apply)
        actions.addWidget(cancel)
        layout.addLayout(actions)
        dialog.exec()

    def _apply_calculation_type_to_runtime(self, module: CalculationModule) -> None:
        runtime = module.runtime
        is_batch = (
            module.kind in {"equilibrium", "solidification"}
            and module.calculation_type == "batch"
        )
        has_imported_batch = is_batch and bool(module.batch_condition)
        for widget in runtime.get("single_condition_widgets", []):
            if isinstance(widget, QWidget):
                widget.setVisible(not has_imported_batch)
        for widget in runtime.get("single_condition_buttons", []):
            if isinstance(widget, QWidget):
                widget.setVisible(not has_imported_batch)
        for widget in runtime.get("batch_visible_widgets", []):
            if isinstance(widget, QWidget):
                widget.setVisible(module.kind == "solidification")

        transition_search = runtime.get("transition_search")
        temperature = runtime.get("temperature")
        if isinstance(transition_search, QCheckBox) and isinstance(
            temperature,
            QLineEdit,
        ):
            transition_search.setVisible(
                (not is_batch)
                and module.kind == "equilibrium"
                and _temperature_text_is_range(temperature.text())
            )
        if isinstance(temperature, QLineEdit):
            temperature.setPlaceholderText(
                "[initial] [Final] [Step]"
                if is_batch and not has_imported_batch
                else "1600"
            )

        table = runtime.get("batch_conditions_table")
        if isinstance(table, QTableWidget):
            if has_imported_batch:
                _populate_batch_condition_table(table, module.batch_condition)
            else:
                table.clear()
                table.setColumnCount(0)
                table.setRowCount(0)
                table.setVisible(False)
        composition_table = runtime.get("composition_table")
        if _is_composition_table(composition_table):
            composition_table.setProperty("amount_unit", module.amount_unit)
            composition_table.setProperty(
                "allow_amount_ranges",
                is_batch and not has_imported_batch,
            )
            _refresh_composition_amount_editors(composition_table)
        self._update_manual_batch_point_count(module)
        self._refresh_nucleation_undercooling_runtime(module)

    def _update_manual_batch_point_count(self, module: CalculationModule) -> None:
        label = module.runtime.get("batch_point_count")
        if not isinstance(label, QLabel):
            return
        is_manual_batch = (
            module.kind in {"equilibrium", "solidification"}
            and module.calculation_type == "batch"
            and not module.batch_condition
        )
        if not is_manual_batch:
            label.setVisible(False)
            label.setText("")
            return
        temperature = module.runtime.get("temperature")
        pressure = module.runtime.get("pressure")
        table = module.runtime.get("composition_table")
        if not (
            isinstance(temperature, QLineEdit)
            and isinstance(pressure, QLineEdit)
            and _is_composition_table(table)
        ):
            label.setVisible(False)
            label.setText("")
            return
        try:
            points = _manual_batch_point_count(
                table,
                temperature.text(),
                pressure.text(),
            )
        except ValueError as exc:
            label.setVisible(True)
            label.setText(f"Batch points: {exc}")
            return
        label.setVisible(True)
        label.setText(f"Batch points: {points}")

    def _refresh_nucleation_undercooling_runtime(
        self,
        module: CalculationModule,
    ) -> None:
        section = module.runtime.get("nucleation_undercooling_section")
        tree = module.runtime.get("nucleation_undercooling_tree")
        visible = (
            module.kind == "solidification"
            and _solidification_model_value(module.solidification_model)
            == "nucleoscheil"
        )
        if isinstance(section, QWidget):
            section.setVisible(visible)
            section.setMaximumHeight(16777215 if visible else 0)
        if not visible or not isinstance(tree, QTreeWidget):
            return
        phase_tree = module.runtime.get("phase_tree")
        phase_names = list(module.phase_names)
        if not phase_names:
            session_id = str(module.runtime.get("session_id") or "")
            session = self.calculation_sessions.get(session_id)
            active_database = (
                self._active_session_database(session)
                if session is not None
                else None
            )
            composition_table = module.runtime.get("composition_table")
            if (
                session is not None
                and active_database is not None
                and _is_composition_table(composition_table)
            ):
                try:
                    import equilipy as eq

                    elements = _module_condition_elements(
                        module,
                        composition_table,
                        active_database.database,
                    )
                    phase_names = [
                        phase
                        for phase in eq.list_phases(
                            active_database.database,
                            elements,
                        )
                        if not _is_reference_ser_phase(phase)
                    ]
                    module.phase_names = list(phase_names)
                except Exception:
                    phase_names = []
        if isinstance(phase_tree, QTreeWidget) and _all_phase_names(phase_tree):
            phase_names = _selected_phase_names(phase_tree)
        liquid_widget = module.runtime.get("liquid_phase")
        liquid_name = "LIQUID"
        if isinstance(liquid_widget, QLineEdit):
            liquid_name = liquid_widget.text().strip() or liquid_name
        _populate_nucleation_undercooling_tree(
            tree,
            phase_names,
            module.nucleation_undercooling,
            liquid_name,
            module.temperature_unit,
        )

    def _nucleation_placeholder_page(
        self,
        module: CalculationModule,
    ) -> QWidget:
        page = QScrollArea()
        page.setObjectName("EditorScroll")
        page.setWidgetResizable(True)
        page.setFrameShape(QFrame.Shape.NoFrame)

        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(14)

        section = _form_section("Nucleation")
        message = QLabel("To be implemented")
        message.setObjectName("DatabaseListLabel")
        message.setAlignment(Qt.AlignmentFlag.AlignCenter)
        message.setMinimumHeight(120)
        section.layout().addWidget(message)
        layout.addWidget(section)
        layout.addStretch(1)

        module.runtime = {"status": message}
        page.setWidget(content)
        return page

    def _add_composition_table_row(
        self,
        session: CalculationSession,
        table: QTableWidget,
        status: QPlainTextEdit,
    ) -> None:
        _add_composition_row(table)
        _refresh_composition_amount_editors(table)
        self._validate_composition_table(session, table, status, quiet=True)

    def _remove_composition_table_rows(
        self,
        session: CalculationSession,
        table: QTableWidget,
        status: QPlainTextEdit,
    ) -> None:
        _remove_selected_table_rows(table)
        self._validate_composition_table(session, table, status, quiet=True)

    def _validate_composition_table(
        self,
        session: CalculationSession,
        table: QTableWidget,
        status: QPlainTextEdit | None = None,
        *,
        quiet: bool = False,
    ) -> bool:
        active_database = self._active_session_database(session)
        database = active_database.database if active_database is not None else None
        diagnostics = _validate_composition_items(table, database)
        if status is not None and diagnostics and not quiet:
            status.setPlainText("Composition issue: " + diagnostics[0])
        return not diagnostics

    def _revalidate_session_composition_tables(
        self,
        session: CalculationSession,
    ) -> None:
        for module in session.modules.values():
            table = module.runtime.get("composition_table")
            status = module.runtime.get("status")
            if _is_composition_table(table):
                self._validate_composition_table(
                    session,
                    table,
                    status if isinstance(status, QPlainTextEdit) else None,
                    quiet=True,
                )

    def _ensure_calculation_module_runtime(
        self,
        session: CalculationSession,
        module: CalculationModule,
    ) -> dict[str, Any]:
        key = self._module_key(module.id)
        if key not in self.calculation_pages or not module.runtime:
            page = self._calculation_module_page(session, module)
            page.setProperty("session_id", session.id)
            page.setProperty("module_id", module.id)
            page.setProperty("module_kind", module.kind)
            self._set_calculation_page(key, page, make_current=False)
        return module.runtime
