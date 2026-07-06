"""Context snapshots sent to right-sidebar assistant providers."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .conversation_log import recent_assistant_log


@dataclass(frozen=True)
class AssistantContext:
    """Small, serializable snapshot of the current GUI state."""

    workspace: str
    database_name: str = ""
    loaded_databases: tuple[str, ...] = ()
    selected_object: str = ""
    diagnostics: tuple[str, ...] = ()
    calculation_sessions: tuple[str, ...] = ()
    active_calculation_item: str = ""
    status_message: str = ""
    extra: dict[str, str] = field(default_factory=dict)

    def to_prompt_text(self) -> str:
        """Return concise text context for a provider prompt."""
        lines = [
            "Equilipy GUI context:",
            f"- Workspace: {self.workspace or 'Unknown'}",
        ]
        if self.database_name:
            lines.append(f"- Active database: {self.database_name}")
        if self.loaded_databases:
            lines.append(f"- Loaded databases: {', '.join(self.loaded_databases)}")
        if self.selected_object:
            lines.append(f"- Selected object: {self.selected_object}")
        if self.calculation_sessions:
            lines.append(
                f"- Calculation sessions: {', '.join(self.calculation_sessions)}"
            )
        if self.active_calculation_item:
            lines.append(f"- Active calculation item: {self.active_calculation_item}")
        if self.status_message:
            lines.append(f"- Status: {self.status_message}")
        if self.diagnostics:
            lines.append("- Validation diagnostics:")
            lines.extend(f"  - {diagnostic}" for diagnostic in self.diagnostics[:20])
            if len(self.diagnostics) > 20:
                lines.append(
                    f"  - ... {len(self.diagnostics) - 20} additional diagnostic(s)"
                )
        for key, value in sorted(self.extra.items()):
            if value:
                lines.append(f"- {key}: {value}")
        return "\n".join(lines)


def build_assistant_context(window: Any, workspace: str) -> AssistantContext:
    """Build a defensive context snapshot from the main GUI window."""
    database = getattr(window, "current_database", None)
    database_name = str(getattr(database, "name", "") or "")
    loaded_databases = tuple(
        str(getattr(item, "name", "") or "Untitled database")
        for item in getattr(window, "databases", []) or []
    )
    diagnostics = _diagnostic_texts(window, database)
    selected_object = (
        _selected_database_object_text(window) if workspace == "database" else ""
    )
    active_calculation_item = (
        _selected_calculation_item_text(window) if workspace == "calculation" else ""
    )
    calculation_sessions = tuple(
        str(getattr(session, "name", "") or session_id)
        for session_id, session in (
            getattr(window, "calculation_sessions", {}) or {}
        ).items()
    )
    status_message = _status_text(window)
    extra = _extra_context(window)
    return AssistantContext(
        workspace=workspace,
        database_name=database_name,
        loaded_databases=loaded_databases,
        selected_object=selected_object,
        diagnostics=tuple(diagnostics),
        calculation_sessions=calculation_sessions,
        active_calculation_item=active_calculation_item,
        status_message=status_message,
        extra=extra,
    )


def _extra_context(window: Any) -> dict[str, str]:
    extra: dict[str, str] = {}
    recent_log = recent_assistant_log(window)
    if recent_log:
        extra["Recent assistant conversation"] = recent_log
    active_directory = getattr(window, "active_calculation_directory", None)
    if active_directory is None:
        extra["Active calculation directory"] = "Not set"
        return extra
    root = Path(active_directory).expanduser().resolve()
    database_dir = root / "database"
    extra["Active calculation directory"] = str(root)
    extra["Project database folder"] = str(database_dir)
    files = _database_files_in_project(root)
    if files:
        extra["Project database files"] = ", ".join(files)
    else:
        extra["Project database files"] = "None found"
    return extra


def _database_files_in_project(root: Path) -> list[str]:
    names: list[str] = []
    for folder_name in ("database", "databases"):
        directory = root / folder_name
        if not directory.exists():
            continue
        for path in sorted(directory.iterdir(), key=lambda item: item.name.lower()):
            if path.is_file() and path.suffix.lower() in {".dat", ".tdb"}:
                names.append(path.name)
    return names


def _diagnostic_texts(window: Any, database: Any) -> list[str]:
    diagnostics: list[Any] = []
    if hasattr(window, "_database_diagnostics"):
        try:
            diagnostics = list(window._database_diagnostics(database))
        except Exception:
            diagnostics = []
    texts: list[str] = []
    for diagnostic in diagnostics:
        severity = str(getattr(diagnostic, "severity", "") or "").strip()
        message = str(getattr(diagnostic, "message", diagnostic) or "").strip()
        if severity and message:
            texts.append(f"{severity}: {message}")
        elif message:
            texts.append(message)
    return texts


def _selected_database_object_text(window: Any) -> str:
    tree = getattr(window, "tree", None)
    if tree is None:
        return ""
    try:
        index = tree.currentIndex()
        if not index.isValid():
            return ""
        return str(index.data() or "")
    except Exception:
        return ""


def _selected_calculation_item_text(window: Any) -> str:
    item = getattr(window, "last_calculation_item", None)
    if item is None:
        tree = getattr(window, "calculation_tree", None)
        try:
            item = tree.currentItem() if tree is not None else None
        except Exception:
            item = None
    if item is None:
        return ""
    try:
        return str(item.text(0) or "")
    except Exception:
        return ""


def _status_text(window: Any) -> str:
    status = getattr(window, "status_text", None)
    if status is None:
        return ""
    try:
        return str(status.text() or "")
    except Exception:
        return ""
