"""Qt right-sidebar assistant panel."""

from __future__ import annotations

import os
from typing import Any

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtGui import QKeyEvent, QTextCursor
from PySide6.QtWidgets import (
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)

from .actions import (
    AssistantActionError,
    execute_assistant_actions,
    parse_action_requests,
    strip_action_markup,
)
from .context import build_assistant_context
from .conversation_log import append_assistant_log
from .providers import (
    AssistantProvider,
    AssistantProviderError,
    LocalAssistantProvider,
    available_providers,
    create_cli_provider,
)

ASSISTANT_PANEL_MAX_WIDTH = 900


class AssistantRequestWorker(QObject):
    """Run one assistant provider request off the GUI thread."""

    chunk = Signal(str)
    error = Signal(str)
    finished = Signal()

    def __init__(
        self,
        provider: AssistantProvider,
        context,
        message: str,
        cwd: str | os.PathLike[str] | None,
    ) -> None:
        super().__init__()
        self._provider = provider
        self._context = context
        self._message = message
        self._cwd = cwd
        self._canceled = False

    @Slot()
    def run(self) -> None:
        """Stream provider output into main-thread signals."""
        try:
            for text in self._provider.stream_reply(
                self._context,
                self._message,
                cwd=self._cwd,
            ):
                if self._canceled:
                    break
                self.chunk.emit(str(text))
        except AssistantProviderError as exc:
            if not self._canceled:
                self.error.emit(str(exc))
        except Exception as exc:  # pragma: no cover - defensive GUI path.
            if not self._canceled:
                self.error.emit(f"Assistant request failed: {exc}")
        finally:
            self.finished.emit()

    @Slot()
    def cancel(self) -> None:
        """Cancel the request and terminate any active subprocess."""
        self._canceled = True
        self._provider.cancel()


def _composer_placeholder_for_status(status: str) -> str:
    """Map request status to compact composer placeholder text."""
    normalized = status.strip().lower()
    if normalized.startswith("request failed"):
        return "Request failed - Ask Equilipy"
    if "canceled" in normalized:
        return "Canceled - Ask Equilipy"
    return "Ready - Ask Equilipy"


class AssistantPromptEdit(QPlainTextEdit):
    """Small multi-line composer where Enter submits and Shift+Enter inserts."""

    submitRequested = Signal()

    def __init__(self) -> None:
        super().__init__()
        self.setMinimumHeight(128)
        self.setMaximumHeight(128)
        self.setTabChangesFocus(True)

    def text(self) -> str:
        """Return the composer text, matching the old QLineEdit call site."""
        return self.toPlainText()

    def setText(self, text: str) -> None:
        """Set composer text, matching the old QLineEdit call site."""
        self.setPlainText(text)

    def keyPressEvent(self, event: QKeyEvent) -> None:
        """Submit on Enter, insert a newline on Shift+Enter."""
        if event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                super().keyPressEvent(event)
                return
            self.submitRequested.emit()
            event.accept()
            return
        super().keyPressEvent(event)


class AssistantPanel(QFrame):
    """Right-side assistant panel shared by Database and Calculation workspaces."""

    def __init__(self, owner: Any, workspace: str) -> None:
        super().__init__()
        self.owner = owner
        self.workspace = workspace
        self._thread: QThread | None = None
        self._worker: AssistantRequestWorker | None = None
        self._provider_objects: dict[str, AssistantProvider] = {
            "local": LocalAssistantProvider(),
        }
        self._request_outcome_status = "Ready."
        self._composer_placeholder = "Ready - Ask Equilipy"
        self._active_status_text = ""
        self._active_user_message = ""
        self._active_provider_label = ""
        self._assistant_reply_buffer = ""
        self._assistant_reply_start: int | None = None

        self.setObjectName("ChatPanel")
        self.setMinimumWidth(280)
        self.setMaximumWidth(ASSISTANT_PANEL_MAX_WIDTH)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(12)

        title = QLabel("Assistant")
        title.setObjectName("PanelTitle")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        provider_row = QHBoxLayout()
        provider_row.setSpacing(8)
        provider_label = QLabel("Provider")
        provider_label.setObjectName("InlineLabel")
        self.provider_select = QComboBox()
        self.provider_select.setObjectName("AssistantProviderSelect")
        provider_row.addWidget(provider_label)
        provider_row.addWidget(self.provider_select, 1)
        layout.addLayout(provider_row)

        self.transcript = QPlainTextEdit()
        self.transcript.setObjectName("ChatTranscript")
        self.transcript.setReadOnly(True)
        self.transcript.setPlaceholderText("Chat transcript")
        layout.addWidget(self.transcript, 1)

        input_row = QHBoxLayout()
        input_row.setSpacing(8)
        self.prompt = AssistantPromptEdit()
        self.prompt.setObjectName("ChatInput")
        self.prompt.setPlaceholderText(self._composer_placeholder)
        self.send = QPushButton("Send")
        self.send.setObjectName("SmallActionButton")
        self.refresh = QPushButton("Refresh")
        self.refresh.setObjectName("SmallActionButton")
        self.stop = QPushButton("Cancel")
        self.stop.setObjectName("SmallActionButton")
        self.stop.setEnabled(False)
        button_column = QVBoxLayout()
        button_column.setSpacing(8)
        button_column.addWidget(self.send)
        button_column.addWidget(self.refresh)
        button_column.addWidget(self.stop)
        input_row.addWidget(self.prompt, 1)
        input_row.addLayout(button_column)
        layout.addLayout(input_row)

        self.send.clicked.connect(self.send_message)
        self.prompt.submitRequested.connect(self.send_message)
        self.refresh.clicked.connect(self._refresh_providers)
        self.stop.clicked.connect(self.cancel_request)
        self._refresh_providers()
        self.setVisible(False)

    def _refresh_providers(self) -> None:
        selected_name = str(self.provider_select.currentData() or "local")
        self.provider_select.clear()
        self.provider_select.addItem("Local", "local")
        self._provider_objects = {
            "local": self._provider_objects.get("local", LocalAssistantProvider())
        }
        for info in available_providers():
            if info.available:
                self._provider_objects[info.name] = create_cli_provider(info)
                self.provider_select.addItem(info.label, info.name)
        for index in range(self.provider_select.count()):
            if str(self.provider_select.itemData(index) or "") == selected_name:
                self.provider_select.setCurrentIndex(index)
                break
        self._update_status("Ready.")

    def _selected_provider(self) -> AssistantProvider:
        name = str(self.provider_select.currentData() or "local")
        return self._provider_objects.get(name, self._provider_objects["local"])

    @Slot()
    def send_message(self) -> None:
        """Send the prompt to the selected provider."""
        message = self.prompt.text().strip()
        if not message or self._thread is not None:
            return
        self._active_user_message = message
        self.prompt.clear()
        prefix = "\n" if self.transcript.toPlainText().strip() else ""
        self._append_transcript(f"{prefix}You: {message}\n\n")
        provider = self._selected_provider()
        self._active_provider_label = provider.label
        context = build_assistant_context(self.owner, self.workspace)
        self._log_conversation("user", message, provider.label)
        if isinstance(provider, LocalAssistantProvider):
            self._append_transcript(f"{provider.label}: ")
            reply_chunks: list[str] = []
            for chunk in provider.stream_reply(
                context, message, cwd=self._working_dir()
            ):
                reply_chunks.append(str(chunk))
                self._append_transcript(chunk)
            self._append_transcript("\n")
            self._log_conversation(
                "assistant",
                "".join(reply_chunks).strip(),
                provider.label,
            )
            self._update_status("Ready.")
            return

        self._append_transcript(f"{provider.label}: ")
        self._assistant_reply_buffer = ""
        self._assistant_reply_start = len(self.transcript.toPlainText())
        self._request_outcome_status = "Ready."
        self._active_status_text = f"{provider.label} running..."
        self._set_running(True)
        thread = QThread(self)
        worker = AssistantRequestWorker(provider, context, message, self._working_dir())
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.chunk.connect(self._append_provider_chunk)
        worker.error.connect(self._append_error)
        worker.finished.connect(thread.quit)
        worker.finished.connect(self._finish_request)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._request_thread_finished)
        self._thread = thread
        self._worker = worker
        thread.start()
        self._update_status(f"{provider.label} running...")

    @Slot()
    def cancel_request(self) -> None:
        """Cancel the active assistant request."""
        if self._worker is not None:
            self._request_outcome_status = "Assistant request canceled."
            self._worker.cancel()
        self._update_status("Assistant request canceled.")

    @Slot(str)
    def _append_transcript(self, text: str) -> None:
        cursor = self.transcript.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        cursor.insertText(str(text))
        self.transcript.setTextCursor(cursor)
        self.transcript.ensureCursorVisible()

    @Slot(str)
    def _append_provider_chunk(self, text: str) -> None:
        self._assistant_reply_buffer += str(text)
        self._append_transcript(text)

    @Slot(str)
    def _append_error(self, text: str) -> None:
        self._append_transcript(f"\n[{text}]\n")
        self._log_conversation(
            "assistant-error", str(text), self._active_provider_label
        )
        self._request_outcome_status = "Request failed. See transcript."
        self._update_status(self._request_outcome_status)

    @Slot()
    def _finish_request(self) -> None:
        if self._request_outcome_status == "Ready.":
            self._handle_completed_assistant_reply()
        else:
            self._append_transcript("\n")
        self._set_running(False)

    @Slot()
    def _request_thread_finished(self) -> None:
        """Release request objects after Qt has actually stopped the thread."""
        self._thread = None
        self._worker = None
        self._update_status(self._request_outcome_status)

    def _handle_completed_assistant_reply(self) -> None:
        """Execute provider-requested action blocks after a clean response."""
        reply = self._strip_echoed_prompt_from_reply(self._assistant_reply_buffer)
        visible_reply = strip_action_markup(reply).rstrip()
        try:
            actions = parse_action_requests(reply)
        except AssistantActionError as exc:
            self._append_transcript(f"\n[Assistant action request failed: {exc}]\n")
            self._log_conversation(
                "assistant-error",
                f"{visible_reply}\nAssistant action request failed: {exc}".strip(),
                self._active_provider_label,
            )
            self._request_outcome_status = "Action failed. See transcript."
            return
        if not actions:
            if visible_reply != self._assistant_reply_buffer.rstrip():
                self._replace_assistant_reply(visible_reply)
            self._append_transcript("\n")
            self._log_conversation(
                "assistant",
                visible_reply,
                self._active_provider_label,
            )
            return

        self._replace_assistant_reply(visible_reply)
        if visible_reply:
            self._append_transcript("\n\n")
        try:
            summaries = execute_assistant_actions(self.owner, actions)
        except AssistantActionError as exc:
            self._append_transcript(f"Equilipy: Action failed: {exc}\n")
            self._log_conversation(
                "assistant",
                f"{visible_reply}\n\nEquilipy: Action failed: {exc}".strip(),
                self._active_provider_label,
            )
            self._request_outcome_status = "Action failed. See transcript."
            return
        if summaries:
            summary_text = "Equilipy: " + " ".join(summaries)
            self._append_transcript(summary_text + "\n")
        else:
            summary_text = "Equilipy: No GUI action was requested."
            self._append_transcript(summary_text + "\n")
        self._log_conversation(
            "assistant",
            f"{visible_reply}\n\n{summary_text}".strip(),
            self._active_provider_label,
        )

    def _strip_echoed_prompt_from_reply(self, text: str) -> str:
        """Remove a provider prompt echo before display and action parsing."""
        reply = str(text or "")
        marker = "User request:\n"
        marker_index = reply.rfind(marker)
        if marker_index < 0:
            return reply
        prefix = reply[:marker_index]
        if (
            "You are the right-side assistant inside the Equilipy GUI" not in prefix
            and "Equilipy GUI context:" not in prefix
        ):
            return reply
        after_marker = reply[marker_index + len(marker) :].lstrip()
        user_message = self._active_user_message.strip()
        if user_message and after_marker.startswith(user_message):
            after_marker = after_marker[len(user_message) :].lstrip()
        return after_marker

    def _replace_assistant_reply(self, text: str) -> None:
        start = self._assistant_reply_start
        if start is None:
            return
        current_text = self.transcript.toPlainText()
        self.transcript.setPlainText(current_text[:start] + str(text))
        cursor = self.transcript.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self.transcript.setTextCursor(cursor)
        self.transcript.ensureCursorVisible()

    def _set_running(self, running: bool) -> None:
        self.send.setEnabled(not running)
        if running:
            self.prompt.setText(self._active_status_text)
        else:
            self.prompt.clear()
            self.prompt.setPlaceholderText(self._composer_placeholder)
        self.prompt.setEnabled(not running)
        self.stop.setEnabled(running)
        self.refresh.setEnabled(not running)
        self.provider_select.setEnabled(not running)

    def _update_status(self, status: str) -> None:
        if self._thread is not None:
            return
        self._composer_placeholder = _composer_placeholder_for_status(status)
        if not self.prompt.text():
            self.prompt.setPlaceholderText(self._composer_placeholder)

    def _working_dir(self) -> str | None:
        active_directory = getattr(self.owner, "active_calculation_directory", None)
        if active_directory:
            return str(active_directory)
        database = getattr(self.owner, "current_database", None)
        metadata = getattr(database, "metadata", {}) or {}
        for key in (
            "_loaded_path",
            "source_file",
            "loaded_path",
            "source_path",
            "path",
        ):
            value = metadata.get(key)
            if value:
                path = os.path.abspath(os.path.expanduser(str(value)))
                return os.path.dirname(path) if os.path.isfile(path) else path
        return os.getcwd()

    def _log_conversation(
        self,
        role: str,
        text: str,
        provider_label: str | None = None,
    ) -> None:
        try:
            append_assistant_log(
                self.owner,
                workspace=self.workspace,
                provider=provider_label or self._active_provider_label,
                role=role,
                text=text,
            )
        except OSError:
            return
