"""Persistent conversation log for GUI assistant providers."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

LOG_DIR_NAME = ".equilipy"
LOG_FILE_NAME = "assistant_conversation.md"
MAX_PROMPT_HISTORY_CHARS = 6000


def assistant_log_path(owner: Any) -> Path:
    """Return the shared assistant log path for the current GUI context."""
    active_directory = getattr(owner, "active_calculation_directory", None)
    if active_directory:
        root = Path(active_directory).expanduser().resolve()
    else:
        root = Path.cwd().resolve()
    return root / LOG_DIR_NAME / LOG_FILE_NAME


def append_assistant_log(
    owner: Any,
    *,
    workspace: str,
    provider: str,
    role: str,
    text: str,
) -> None:
    """Append one assistant conversation entry to the shared markdown log."""
    content = str(text or "").strip()
    if not content:
        return
    path = assistant_log_path(owner)
    path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    entry = (
        f"\n## {timestamp} | {workspace or 'unknown'} | {provider or 'unknown'} | "
        f"{role}\n\n"
        f"{content}\n"
    )
    with path.open("a", encoding="utf-8") as handle:
        handle.write(entry)


def recent_assistant_log(
    owner: Any, *, max_chars: int = MAX_PROMPT_HISTORY_CHARS
) -> str:
    """Return a bounded recent conversation tail for provider prompts."""
    try:
        text = assistant_log_path(owner).read_text(encoding="utf-8")
    except OSError:
        return ""
    text = text.strip()
    if len(text) <= max_chars:
        return text
    return "... previous assistant conversation truncated ...\n" + text[-max_chars:]
