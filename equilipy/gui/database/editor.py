"""Database workspace construction for the main GUI window."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QFrame,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ..models import DiagnosticsModel

SIDEBAR_MAX_WIDTH = 900


class DatabaseWorkspaceMixin:
    """Build the database workspace shell for ``DatabaseEditorWindow``."""

    def _build_database_workspace(self) -> QWidget:
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setObjectName("WorkspaceSplitter")
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self._build_database_sidebar())
        splitter.addWidget(self._build_database_editor())
        splitter.addWidget(self._build_chat_sidebar("database"))
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 0)
        splitter.setSizes([330, 850, 320])
        return splitter

    def _build_database_sidebar(self) -> QWidget:
        self.database_sidebar = QFrame()
        self.database_sidebar.setObjectName("SidebarPanel")
        self.database_sidebar.setMinimumWidth(290)
        self.database_sidebar.setMaximumWidth(SIDEBAR_MAX_WIDTH)
        layout = QVBoxLayout(self.database_sidebar)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        actions = QHBoxLayout()
        actions.setSpacing(8)
        for label_text in ("Open", "Validate", "Export"):
            button = QPushButton(label_text)
            button.setObjectName("SmallActionButton")
            button.clicked.connect(
                lambda _checked=False, text=label_text: self._sidebar_action(text)
            )
            actions.addWidget(button)
        layout.addLayout(actions)

        search = QLineEdit()
        search.setObjectName("SidebarSearch")
        search.setPlaceholderText("Filter database objects")
        layout.addWidget(search)

        self.tree.setObjectName("DatabaseTree")
        self.tree.setModel(self.tree_model)
        self.tree.setHeaderHidden(True)
        self.tree.setAlternatingRowColors(False)
        self._expand_database_roots()
        self.tree.setDragEnabled(True)
        self.tree.setAcceptDrops(True)
        self.tree.setDropIndicatorShown(True)
        self.tree.setDefaultDropAction(Qt.DropAction.CopyAction)
        self.tree.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)
        self.tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._show_database_tree_menu)
        self.tree.selectionModel().currentChanged.connect(self._show_selected)
        self.editor_title.editingFinished.connect(self._commit_function_name)
        layout.addWidget(self.tree, 1)
        return self.database_sidebar

    def _build_database_editor(self) -> QWidget:
        self.database_editor = QFrame()
        self.database_editor.setObjectName("EditorPanel")
        layout = QVBoxLayout(self.database_editor)
        layout.setContentsMargins(22, 20, 22, 20)
        layout.setSpacing(14)

        header = QHBoxLayout()
        title_block = QVBoxLayout()
        title_block.setSpacing(4)

        self.editor_badge.setObjectName("EditorBadge")
        self.editor_title.setObjectName("EditorTitle")
        self.editor_title.setReadOnly(True)
        self.editor_title.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.editor_subtitle.setObjectName("EditorSubtitle")
        self.editor_subtitle.setWordWrap(True)
        self.editor_subtitle.setMinimumWidth(0)
        self.editor_subtitle.setSizePolicy(
            QSizePolicy.Policy.Ignored,
            QSizePolicy.Policy.Preferred,
        )
        self.editor_subtitle_control.setObjectName("EditorSubtitleControl")
        self.editor_subtitle_control_layout = QHBoxLayout(self.editor_subtitle_control)
        self.editor_subtitle_control_layout.setContentsMargins(0, 0, 0, 0)
        self.editor_subtitle_control_layout.setSpacing(8)
        self.editor_subtitle_control.hide()
        title_block.addWidget(self.editor_badge)
        title_block.addWidget(self.editor_title)
        title_block.addWidget(self.editor_subtitle)
        title_block.addWidget(self.editor_subtitle_control)
        header.addLayout(title_block, 1)

        layout.addLayout(header)

        vertical_splitter = QSplitter(Qt.Orientation.Vertical)
        vertical_splitter.setChildrenCollapsible(False)

        scroll = QScrollArea()
        scroll.setObjectName("EditorScroll")
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.editor_form_layout.setContentsMargins(0, 0, 0, 0)
        self.editor_form_layout.setSpacing(14)
        self.editor_form_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        scroll.setWidget(self.editor_form)
        vertical_splitter.addWidget(scroll)

        self.database_inspector_tabs = QTabWidget()
        self.database_inspector_tabs.setObjectName("InspectorTabs")
        self.source_preview.setReadOnly(True)
        self.source_preview.setObjectName("SourcePreview")
        self.source_preview.setPlaceholderText("Source preview will appear here.")
        self.database_inspector_tabs.addTab(self.source_preview, "Source")

        diagnostics_table = QTableView()
        self.database_diagnostics_table = diagnostics_table
        diagnostics_table.setObjectName("DiagnosticsTable")
        diagnostics_table.setModel(DiagnosticsModel(self._database_diagnostics()))
        diagnostics_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows
        )
        diagnostics_table.resizeColumnsToContents()
        self.database_inspector_tabs.addTab(diagnostics_table, "Diagnostics")
        vertical_splitter.addWidget(self.database_inspector_tabs)
        self.database_inspector_tabs.setVisible(False)
        vertical_splitter.setStretchFactor(0, 4)
        vertical_splitter.setStretchFactor(1, 1)
        vertical_splitter.setSizes([520, 180])
        layout.addWidget(vertical_splitter, 1)
        return self.database_editor
