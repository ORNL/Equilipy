"""Stylesheet helpers for the Equilipy GUI."""

from __future__ import annotations

import re
from pathlib import Path


def app_style_sheet(font_delta: int = 0) -> str:
    """Return the application stylesheet, adjusted by a font zoom delta."""
    _checkmark = str(
        Path(__file__).resolve().parent / "icons" / "checkmark.svg"
    ).replace("\\", "/")
    sheet = """
QMainWindow#DatabaseEditorWindow {
    background: #1f2526;
    color: #eceff0;
}

QWidget#RootSurface {
    background: #1f2526;
    color: #eceff0;
    font-size: 14px;
}

QFrame#AppHeader {
    background: #222829;
    border-bottom: 1px solid #15191a;
}

QFrame#AppStatusBar {
    background: #1b2021;
    border-top: 1px solid #15191a;
}

QLabel#AppStatusText {
    color: #cdd4d6;
    font-size: 12px;
    font-weight: 600;
}

QLabel#AppStatusSymbol {
    background: #7ac143;
    border-radius: 4px;
}

QLabel#AppStatusSymbol[status="running"] {
    background: #1687ff;
}

QLabel#AppStatusSymbol[status="warning"] {
    background: transparent;
}

QLabel#AppStatusSymbol[status="error"] {
    background: transparent;
}

QSvgWidget#StatusOrnlLogo {
    background: transparent;
}

QLabel#StatusOrnlText {
    color: #aeb8ba;
    font-size: 12px;
    font-weight: 700;
}

QFrame#HeaderBrand {
    background: transparent;
}

QFrame#HeaderBrandText {
    background: transparent;
}

QSvgWidget#HeaderLogo {
    background: transparent;
}

QLabel#HeaderBrandName {
    color: #eef2f2;
    font-size: 34px;
    font-weight: 800;
}

QLabel#HeaderBrandVersion {
    color: #aeb8ba;
    font-size: 13px;
    font-weight: 700;
}

QFrame#LayoutToggleGroup {
    background: transparent;
}

QPushButton#LayoutToggleButton {
    background: transparent;
    border: 0;
    border-radius: 0;
    padding: 0;
}

QPushButton#LayoutToggleButton:hover {
    background: transparent;
}

QPushButton#LayoutToggleButton:checked {
    background: transparent;
}

QFrame#ModeSelector {
    background: transparent;
    border: 0;
    border-radius: 0;
}

QPushButton#ModeButton {
    background: #263033;
    color: #b7bec0;
    border: 1px solid #3b4648;
    border-radius: 9px;
    padding: 5px 18px;
    font-weight: 500;
    font-size: 18px;
}

QPushButton#ModeButton:hover {
    background: #30393b;
    color: #edf2f4;
}

QPushButton#ModeButton:checked {
    background: #1687ff;
    border-color: #1687ff;
    color: white;
}

QPushButton#ThermoModeButton {
    background: #263033;
    color: #b7bec0;
    border: 1px solid #3b4648;
    border-radius: 8px;
    padding: 4px 14px;
    font-weight: 600;
    font-size: 14px;
}

QPushButton#ThermoModeButton:hover {
    background: #30393b;
    color: #edf2f4;
}

QPushButton#ThermoModeButton:checked {
    background: #1687ff;
    border-color: #1687ff;
    color: white;
}

QToolButton#GibbsBalanceCoefficientButton {
    background: #263033;
    color: #dce2e4;
    border: 1px solid #4b575a;
    border-radius: 5px;
    padding: 3px 7px;
    font-size: 12px;
    font-weight: 700;
}

QToolButton#GibbsBalanceCoefficientButton:hover {
    background: #30393b;
    color: #ffffff;
    border-color: #1687ff;
}

QLineEdit#SidebarSearch {
    background: #202728;
    color: #e8ebec;
    border: 1px solid #3f4b4e;
    border-radius: 18px;
    padding: 8px 14px;
    selection-background-color: #1687ff;
}

QSplitter#WorkspaceSplitter::handle {
    background: #1a1f20;
}

QFrame#SidebarPanel {
    background: #202627;
    border-right: 1px solid #171b1c;
}

QFrame#ChatPanel {
    background: #202627;
    border-left: 1px solid #171b1c;
}

QLabel#PanelTitle {
    color: #f0f2f3;
    font-size: 17px;
    font-weight: 700;
}

QLabel#InlineLabel {
    color: #aeb6b8;
}

QComboBox#AssistantProviderSelect {
    background: #202728;
    color: #e8ebec;
    border: 1px solid #465153;
    border-radius: 8px;
    padding: 6px 10px;
}

QTreeView#DatabaseTree,
QTreeWidget#DatabaseTree {
    background: #202627;
    color: #d8dcdd;
    border: 1px solid #4a5457;
    border-radius: 4px;
    padding: 8px;
    outline: 0;
    show-decoration-selected: 1;
}

QTreeView#DatabaseTree::item,
QTreeWidget#DatabaseTree::item {
    min-height: 28px;
    padding: 4px 8px;
    border-radius: 10px;
}

QTreeView#DatabaseTree::item:hover,
QTreeWidget#DatabaseTree::item:hover {
    background: #31383a;
}

QTreeView#DatabaseTree::item:selected,
QTreeWidget#DatabaseTree::item:selected {
    background: #4a4d4f;
    color: #ffffff;
}

QPushButton#SmallActionButton,
QPushButton#PrimaryButton,
QToolButton#SmallMenuButton {
    background: #30383a;
    color: #e7ebec;
    border: 1px solid #566063;
    border-radius: 8px;
    padding: 7px 12px;
    font-weight: 700;
    font-size: 14px;
}

QPushButton#SmallActionButton:hover,
QPushButton#PrimaryButton:hover,
QToolButton#SmallMenuButton:hover {
    border-color: #1687ff;
    color: #ffffff;
}

QPushButton#ResultViewButton {
    background: #30383a;
    color: #d5dcde;
    border: 1px solid #566063;
    border-radius: 8px;
    padding: 7px 12px;
    font-weight: 800;
    font-size: 14px;
}

QPushButton#ResultViewButton:hover {
    border-color: #1687ff;
    color: #ffffff;
}

QPushButton#ResultViewButton:checked {
    background: #1687ff;
    border-color: #1687ff;
    color: #ffffff;
}

QPushButton#ResultViewButton:disabled {
    background: #2b3234;
    color: #788285;
    border-color: #40494b;
}

QToolButton#SmallMenuButton::menu-indicator {
    image: none;
    width: 0;
}

QToolButton#CalculationModuleButton {
    background: #263033;
    border: 1px solid #3b4648;
    border-radius: 12px;
    padding: 2px;
}

QToolButton#CalculationModuleButton:hover {
    background: #30393b;
    border-color: #1687ff;
    color: #ffffff;
}

QToolButton#IconActionButton {
    background: #30383a;
    border: 1px solid #566063;
    border-radius: 8px;
    min-width: 30px;
    min-height: 30px;
    max-width: 30px;
    max-height: 30px;
    padding: 0;
}

QToolButton#IconActionButton:hover {
    border-color: #1687ff;
    background: #353f41;
}

QToolButton#DatabaseRemoveButton {
    background: #30383a;
    border: 1px solid #566063;
    border-radius: 5px;
    padding: 0;
    margin: 0;
}

QToolButton#DatabaseRemoveButton:hover {
    border-color: #1687ff;
    background: #353f41;
}

QFrame#ResultColumnGroup {
    background: transparent;
    border: 0;
}

QToolButton#ResultColumnGroupToggle {
    background: transparent;
    border: 0;
    color: #eef2f3;
    font-weight: 800;
    font-size: 16px;
    padding: 4px 0;
}

QToolButton#ResultColumnGroupToggle:hover {
    color: #ffffff;
}

QToolButton#ResultOverviewCollapseButton,
QToolButton#SectionCollapseButton {
    background: transparent;
    border: 0;
    color: #e7ebec;
    min-width: 22px;
    min-height: 22px;
    padding: 0;
}

QToolButton#ResultOverviewCollapseButton:hover,
QToolButton#SectionCollapseButton:hover {
    background: #30383a;
    border-radius: 5px;
}

QPushButton#NewSessionButton {
    background: #30383a;
    color: #f0f5f6;
    border: 1px solid #566063;
    border-radius: 10px;
    min-height: 36px;
    padding: 7px 12px;
    font-size: 14px;
    font-weight: 800;
    text-align: left;
}

QPushButton#NewSessionButton:hover {
    border-color: #1687ff;
    color: #ffffff;
}

QLabel#NewSessionLabel {
    color: #e5e9ea;
    font-size: 14px;
    font-weight: 700;
}

QLabel#DatabaseListLabel {
    color: #e5e9ea;
    font-size: 14px;
}

QLabel#DialogRowLabel,
QRadioButton#DialogRadioOption {
    color: #e5e9ea;
    font-size: 14px;
    font-weight: 700;
}

QScrollArea#DatabaseListScroll {
    background: transparent;
    border: 0;
}

QWidget#DatabaseListViewport {
    background: transparent;
}

QFrame#DatabaseListItem {
    background: #30383a;
    border: 1px solid #465153;
    border-radius: 8px;
}

QPushButton#PrimaryButton {
    background: #1687ff;
    border-color: #1687ff;
    color: white;
    padding-left: 18px;
    padding-right: 18px;
}

QPushButton#PrimaryButton:hover {
    background: #2d9cff;
    border-color: #64bdff;
    color: white;
}

QPushButton#PrimaryButton:pressed {
    background: #0f72d7;
    border-color: #0f72d7;
}

QFrame#EditorPanel {
    background: #252b2c;
    border: 1px solid #3e484a;
    border-radius: 8px;
    margin: 18px;
}

QLabel#EditorBadge {
    color: #1687ff;
    font-size: 13px;
    font-weight: 700;
}

QLabel#EditorTitle {
    color: #f0f2f3;
    font-size: 24px;
    font-weight: 800;
}

QLabel#EditorSubtitle {
    color: #aeb6b8;
    font-size: 14px;
}

QScrollArea#EditorScroll {
    background: transparent;
    border: 0;
}

QFrame#FormSection {
    background: #282f30;
    border: 1px solid #475254;
    border-radius: 8px;
}

QLabel#SectionTitle {
    color: #e9eded;
    font-size: 16px;
    font-weight: 700;
}

QLabel {
    color: #d7dcdd;
}

QLineEdit,
QPlainTextEdit,
QComboBox {
    background: #3a4143;
    color: #eef1f2;
    border: 1px solid #4e595c;
    border-radius: 7px;
    padding: 7px 10px;
    selection-background-color: #1687ff;
}

QLineEdit#EditorBadge {
    background: transparent;
    color: #1687ff;
    border: 0;
    border-radius: 0;
    padding: 0;
    font-size: 13px;
    font-weight: 700;
    selection-background-color: #2a5d83;
}

QLineEdit#EditorTitle {
    background: transparent;
    color: #f0f2f3;
    border: 0;
    border-radius: 0;
    padding: 0;
    font-size: 24px;
    font-weight: 800;
    selection-background-color: #2a5d83;
}

QLineEdit#CompoundPhaseNameInlineInput {
    background: transparent;
    color: #aeb6b8;
    border: 0;
    border-radius: 0;
    padding: 0;
    font-size: 14px;
    selection-background-color: #2a5d83;
}

QPlainTextEdit {
    font-family: Menlo;
}

QPlainTextEdit#ChatTranscript {
    background: #202627;
    border: 1px solid #465153;
}

QPlainTextEdit#ChatInput {
    border-radius: 12px;
}

QComboBox::drop-down {
    border: 0;
    width: 26px;
}

QCheckBox {
    color: #d9ddde;
    spacing: 10px;
}

QCheckBox::indicator {
    width: 20px;
    height: 20px;
    border-radius: 5px;
    background: #3b4345;
    border: 1px solid #4e595c;
}

QCheckBox::indicator:checked {
    background: #1687ff;
    border-color: #1687ff;
    image: url(""" + _checkmark + """);
}

QCheckBox#ResultSectionCheck::indicator:indeterminate {
    background: rgba(22, 135, 255, 51);
    border-color: rgba(22, 135, 255, 140);
}

QCheckBox#ResultColumnCheck::indicator,
QCheckBox#ResultSectionCheck::indicator,
QCheckBox#ResultUnitCheck::indicator {
    width: 14px;
    height: 14px;
    border-radius: 4px;
}

QCheckBox#ResultColumnCheck,
QCheckBox#ResultSectionCheck,
QCheckBox#ResultUnitCheck {
    spacing: 7px;
    font-size: 15px;
    padding-top: 3px;
    padding-bottom: 3px;
}

QTabWidget#InspectorTabs::pane {
    background: #252b2c;
    border: 1px solid #475254;
    border-radius: 8px;
}

QTabWidget#InspectorTabs::tab-bar {
    alignment: left;
    left: 8px;
}

QTabBar::tab {
    background: #30383a;
    color: #c5cbcd;
    border: 1px solid #475254;
    border-bottom: 0;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    padding: 8px 16px;
    margin-right: 3px;
}

QTabBar::tab:selected {
    background: #1687ff;
    color: white;
    border-color: #1687ff;
}

QPlainTextEdit#SourcePreview {
    background: #202627;
    border: 0;
}

QTableView,
QTableWidget {
    background: #202627;
    alternate-background-color: #242b2c;
    color: #e5e9ea;
    gridline-color: #3c4648;
    border: 1px solid #465153;
    border-radius: 6px;
    selection-background-color: #4a4d4f;
    selection-color: white;
}

QHeaderView::section {
    background: #30383a;
    color: #f0f2f3;
    border: 0;
    border-right: 1px solid #465153;
    padding: 7px;
    font-weight: 700;
}

QLabel#DevelopersNameCell {
    background: transparent;
    color: #e5e9ea;
}

QFrame#CompositionEditor {
    background: #202627;
    border: 1px solid #465153;
    border-radius: 6px;
}

QFrame#CompositionEditorHeader {
    background: #30383a;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    border-bottom: 1px solid #465153;
}

QLabel#CompositionEditorHeaderLabel {
    color: #f0f2f3;
    font-weight: 700;
    font-size: 13px;
}

QFrame#CompositionEditorHeaderDivider {
    background: #465153;
    border: 0;
    min-height: 24px;
    max-height: 24px;
}

QFrame#CompositionRow {
    background: transparent;
    border: 1px solid transparent;
    border-radius: 7px;
}

QFrame#CompositionRow[selected="true"] {
    background: #30383a;
    border: 1px solid #566063;
}

QLineEdit#CompositionCellInput {
    background: #30383a;
    color: #f0f3f4;
    border: 1px solid #566063;
    border-radius: 6px;
    padding: 6px 10px;
    min-height: 30px;
    selection-background-color: #1687ff;
    placeholder-text-color: #aab3b6;
}

QLineEdit#CompositionCellInput:focus {
    border: 1px solid #1687ff;
}

QWidget#CompositionAmountInputContainer,
QWidget#CompositionRangeInputContainer {
    background: transparent;
    border: 0;
}

QTableWidget#SolutionSublatticeTable {
    background: #202627;
    alternate-background-color: #202627;
    border: 1px solid #465153;
    border-radius: 6px;
    gridline-color: transparent;
    selection-background-color: transparent;
}

QTableWidget#SolutionSublatticeTable::item {
    background: transparent;
    border: 0;
    padding: 8px;
}

QTableWidget#SolutionSublatticeTable QHeaderView::section {
    background: #30383a;
    color: #f0f2f3;
    border: 0;
    border-right: 1px solid #465153;
    border-bottom: 1px solid #465153;
    padding: 12px 10px;
    font-weight: 700;
    font-size: 13px;
}

QLineEdit#SolutionSublatticeNameInput,
QLineEdit#SolutionSublatticeCoefficientInput,
QLineEdit#SolutionSublatticeConstituentsInput {
    background: #30383a;
    color: #f0f3f4;
    border: 1px solid #566063;
    border-radius: 6px;
    padding: 6px 10px;
    selection-background-color: #1687ff;
    placeholder-text-color: #aab3b6;
}

QLineEdit#SolutionSublatticeNameInput:focus,
QLineEdit#SolutionSublatticeCoefficientInput:focus,
QLineEdit#SolutionSublatticeConstituentsInput:focus {
    border: 1px solid #1687ff;
}

QToolButton#BalanceAmountButton {
    background: #30383a;
    color: #d8dee0;
    border: 1px solid #566063;
    border-radius: 6px;
    min-width: 48px;
    min-height: 30px;
    max-width: 48px;
    padding: 0;
    font-weight: 800;
}

QToolButton#BalanceAmountButton:checked {
    background: #1687ff;
    color: #ffffff;
    border-color: #1687ff;
}

QMenuBar {
    background: #222829;
    color: #d6dbdc;
}

QMenuBar::item:selected {
    background: #30383a;
}

QMenu {
    background: #252b2c;
    color: #e8ebec;
    border: 1px solid #475254;
}

QMenu::item:selected {
    background: #1687ff;
    color: white;
}

QScrollBar:vertical {
    background: #252b2c;
    width: 14px;
    margin: 3px;
}

QScrollBar::handle:vertical {
    background: #8c9496;
    border-radius: 6px;
    min-height: 42px;
}

QScrollBar::add-line:vertical,
QScrollBar::sub-line:vertical {
    height: 0;
}
"""
    return re.sub(
        r"font-size: (\d+)px;",
        lambda match: (
            f"font-size: {_zoomed_font_size(int(match.group(1)), font_delta)}px;"
        ),
        sheet,
    )


def _zoomed_font_size(base_size: int, font_delta: int) -> int:
    """Return a bounded zoomed font size for stylesheet generation."""
    return max(9, min(52, base_size + font_delta))
