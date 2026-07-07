"""GUI stdout/stderr capture with targeted platform-noise filtering."""

from __future__ import annotations

import atexit
import os
import sys
import threading
from pathlib import Path

from PySide6.QtCore import QObject, Signal

_FILTERED_LOG_PATTERNS = (
    "TSMSendMessageToUIServer: CFMessagePortSendRequest FAILED",
)


def filter_gui_log_text(text: str) -> str:
    """Remove known GUI-platform noise from complete log text."""
    return "".join(
        line for line in text.splitlines(True) if not _is_filtered_log_line(line)
    )


def _is_filtered_log_line(line: str) -> bool:
    return any(pattern in line for pattern in _FILTERED_LOG_PATTERNS)


def _flush_stream(stream) -> None:
    """Flush a stdio stream, tolerating windowed builds where it is None."""
    if stream is None:
        return
    try:
        stream.flush()
    except (OSError, ValueError):
        pass


class GuiLogCapture(QObject):
    """Tee process stdout/stderr to the GUI log and an optional file."""

    text_received = Signal(str)

    def __init__(self):
        super().__init__()
        self._started = False
        self._reader_thread: threading.Thread | None = None
        self._read_fd: int | None = None
        self._original_stdout_fd: int | None = None
        self._original_stderr_fd: int | None = None
        self._pending_text = ""
        self._lock = threading.Lock()
        self._log_file = None
        self._log_file_path: Path | None = None

    def start(self) -> None:
        """Start capturing stdout and stderr at the file-descriptor level."""
        if self._started:
            return

        _flush_stream(sys.stdout)
        _flush_stream(sys.stderr)
        try:
            self._original_stdout_fd = os.dup(1)
            self._original_stderr_fd = os.dup(2)
            self._read_fd, write_fd = os.pipe()
            os.dup2(write_fd, 1)
            os.dup2(write_fd, 2)
            os.close(write_fd)
        except OSError:
            # Windowed builds on Windows run without a console:
            # sys.stdout/stderr are None and fds 1/2 are not usable.
            # Run without console capture instead of failing startup.
            for fd_attr in ("_original_stdout_fd", "_original_stderr_fd", "_read_fd"):
                fd = getattr(self, fd_attr)
                if fd is not None:
                    try:
                        os.close(fd)
                    except OSError:
                        pass
                    setattr(self, fd_attr, None)
            return
        self._started = True
        self._reader_thread = threading.Thread(
            target=self._read_loop,
            name="EquilipyGuiLogCapture",
            daemon=True,
        )
        self._reader_thread.start()
        atexit.register(self.stop)

    def stop(self) -> None:
        """Restore stdout/stderr and close the log file."""
        if not self._started:
            self._close_log_file()
            return

        _flush_stream(sys.stdout)
        _flush_stream(sys.stderr)
        if self._original_stdout_fd is not None:
            os.dup2(self._original_stdout_fd, 1)
            os.close(self._original_stdout_fd)
            self._original_stdout_fd = None
        if self._original_stderr_fd is not None:
            os.dup2(self._original_stderr_fd, 2)
            os.close(self._original_stderr_fd)
            self._original_stderr_fd = None
        self._started = False
        self._close_log_file()

    def set_log_file(self, path: Path | None) -> None:
        """Write future captured output to path, or disable file logging."""
        with self._lock:
            self._close_log_file_locked()
            self._log_file_path = path
            if path is None:
                return
            path.parent.mkdir(parents=True, exist_ok=True)
            self._log_file = path.open("ab", buffering=0)

    def write_message(self, message: str) -> None:
        """Write an internal GUI status message to terminal and log file."""
        if not message:
            return
        if not message.endswith("\n"):
            message += "\n"
        filtered = filter_gui_log_text(message)
        if not filtered:
            return
        data = filtered.encode(errors="replace")
        self._write_to_original_stdout(data)
        self._write_to_log_file(data)

    @property
    def log_file_path(self) -> Path | None:
        """Return the active log file path."""
        return self._log_file_path

    def _read_loop(self) -> None:
        if self._read_fd is None:
            return
        with os.fdopen(self._read_fd, "rb", buffering=0) as stream:
            while True:
                try:
                    data = stream.read(4096)
                except OSError:
                    break
                if not data:
                    break
                self._handle_captured_bytes(data)

    def _handle_captured_bytes(self, data: bytes) -> None:
        text = data.decode(errors="replace")
        filtered = self._filter_captured_text(text)
        if not filtered:
            return
        encoded = filtered.encode(errors="replace")
        self._write_to_original_stdout(encoded)
        self._write_to_log_file(encoded)
        self.text_received.emit(filtered)

    def _filter_captured_text(self, text: str) -> str:
        self._pending_text += text
        lines = self._pending_text.splitlines(True)
        if lines and not lines[-1].endswith(("\n", "\r")):
            self._pending_text = lines.pop()
        else:
            self._pending_text = ""
        return filter_gui_log_text("".join(lines))

    def _write_to_original_stdout(self, data: bytes) -> None:
        if self._original_stdout_fd is None:
            return
        try:
            os.write(self._original_stdout_fd, data)
        except OSError:
            pass

    def _write_to_log_file(self, data: bytes) -> None:
        with self._lock:
            if self._log_file is None:
                return
            try:
                self._log_file.write(data)
            except OSError:
                pass

    def _close_log_file(self) -> None:
        with self._lock:
            self._close_log_file_locked()

    def _close_log_file_locked(self) -> None:
        if self._log_file is None:
            return
        try:
            self._log_file.close()
        except OSError:
            pass
        finally:
            self._log_file = None
