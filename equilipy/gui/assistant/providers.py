"""Assistant provider detection, prompt building, and CLI streaming."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Protocol

from .action_contract import (
    action_contract_prompt,
    action_response_schema_json,
    structured_response_to_text,
)
from .context import AssistantContext
from .skill_loader import agent_skills_prompt


class AssistantProviderError(RuntimeError):
    """Raised when an assistant provider cannot produce a response."""


class AssistantProvider(Protocol):
    """Streaming assistant provider protocol."""

    name: str
    label: str

    def stream_reply(
        self,
        context: AssistantContext,
        message: str,
        *,
        cwd: str | os.PathLike[str] | None = None,
    ) -> Iterator[str]:
        """Yield response chunks for one user message."""

    def cancel(self) -> None:
        """Cancel the active provider request if one is running."""


@dataclass(frozen=True)
class ProviderInfo:
    """Detected provider metadata for UI status and tests."""

    name: str
    label: str
    executable: str | None
    status: str
    command_preview: tuple[str, ...] = ()

    @property
    def available(self) -> bool:
        """Return whether the provider executable is available."""
        return self.status == "available"


class LocalAssistantProvider:
    """Deterministic fallback provider that never calls a network service."""

    name = "local"
    label = "Local"

    def stream_reply(
        self,
        context: AssistantContext,
        message: str,
        *,
        cwd: str | os.PathLike[str] | None = None,
    ) -> Iterator[str]:
        """Yield a canned reply built from the GUI context, offline."""
        del cwd
        text = message.strip().lower()
        if "validation" in text or "diagnostic" in text or "error" in text:
            if context.diagnostics:
                yield "Current validation diagnostics:\n"
                for diagnostic in context.diagnostics[:10]:
                    yield f"- {diagnostic}\n"
                if len(context.diagnostics) > 10:
                    yield f"- ... {len(context.diagnostics) - 10} more\n"
            else:
                yield "No validation diagnostics are cached for the current database."
            return
        if "selected" in text or "phase" in text or "function" in text:
            if context.selected_object:
                yield f"Selected object: {context.selected_object}"
            elif context.active_calculation_item:
                yield f"Selected calculation item: {context.active_calculation_item}"
            else:
                yield "No selected object is available in the current workspace."
            return
        yield (
            "I can summarize the current GUI state, explain validation errors, "
            "and suggest safe next actions. Select Codex, Claude, or Gemini for "
            "online CLI-backed assistance."
        )

    def cancel(self) -> None:
        """Local responses are immediate, so there is nothing to cancel."""


class CliProvider:
    """Provider that shells out to a logged-in local assistant CLI."""

    def __init__(
        self,
        name: str,
        label: str,
        executable: str,
        command_kind: str,
    ) -> None:
        self.name = name
        self.label = label
        self.executable = executable
        self.command_kind = command_kind
        self._process: subprocess.Popen[str] | None = None

    def stream_reply(
        self,
        context: AssistantContext,
        message: str,
        *,
        cwd: str | os.PathLike[str] | None = None,
    ) -> Iterator[str]:
        """Stream reply chunks from the provider CLI subprocess."""
        prompt = (
            build_claude_user_prompt(context, message, cwd=cwd)
            if self.command_kind == "claude"
            else build_prompt(context, message, cwd=cwd)
        )
        command, stdin_text = cli_command(self.command_kind, self.executable, prompt)
        environment = cli_environment(self.command_kind)
        output_tail: list[str] = []
        try:
            self._process = subprocess.Popen(
                command,
                cwd=str(cwd) if cwd else None,
                env=environment,
                stdin=subprocess.PIPE if stdin_text is not None else None,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )
        except OSError as exc:
            raise AssistantProviderError(
                f"Could not start {self.label}: {exc}"
            ) from exc

        assert self._process.stdout is not None
        if stdin_text is not None and self._process.stdin is not None:
            try:
                self._process.stdin.write(stdin_text)
                self._process.stdin.close()
            except BrokenPipeError:
                pass

        suppress_provider_dump = False
        try:
            for line in self._process.stdout:
                if self.command_kind == "gemini" and suppress_provider_dump:
                    if _is_assistant_stream_json_line(line):
                        suppress_provider_dump = False
                    else:
                        continue
                chunk = parse_stream_chunk(line)
                if chunk:
                    if _is_provider_setup_error(self.command_kind, chunk):
                        output_tail.append(chunk.strip())
                        output_tail = output_tail[-5:]
                        suppress_provider_dump = self.command_kind == "gemini"
                        continue
                    if _should_suppress_provider_text(self.command_kind, chunk):
                        continue
                    output_tail.append(chunk.strip())
                    output_tail = output_tail[-5:]
                    yield chunk
            return_code = self._process.wait()
            if return_code != 0:
                tail = " ".join(part for part in output_tail if part)
                raise AssistantProviderError(
                    _provider_exit_error_message(
                        self.command_kind,
                        self.label,
                        return_code,
                        tail,
                    )
                )
        finally:
            self._process = None

    def cancel(self) -> None:
        """Terminate the running CLI request, if any."""
        process = self._process
        if process is None or process.poll() is not None:
            return
        process.terminate()


def _is_provider_setup_error(kind: str, text: str) -> bool:
    """Return whether a provider line is setup/auth noise, not model output."""
    normalized = str(text or "").lower()
    if kind == "gemini":
        return (
            "please set an auth method" in normalized
            or "gemini_api_key" in normalized
            or "google_genai_use_vertexai" in normalized
            or "google_genai_use_gca" in normalized
            or _is_gemini_capacity_error_text(normalized)
        )
    return False


def _should_suppress_provider_text(kind: str, text: str) -> bool:
    """Return whether provider text is local tool/status chatter."""
    normalized = " ".join(str(text or "").lower().split())
    if kind == "gemini":
        short_line = len(normalized) < 160
        return (
            _is_gemini_raw_error_dump_text(normalized)
            or (
                short_line
                and (
                    "ripgrep is not available" in normalized
                    or "falling back to greptool" in normalized
                    or normalized.startswith("equilipy gui context:")
                    or normalized.startswith("allowed action names are")
                )
            )
        )
    return False


def _is_gemini_capacity_error_text(normalized: str) -> bool:
    """Return whether Gemini reported temporary API capacity/rate limiting."""
    return (
        "no capacity available for model" in normalized
        or "too many requests" in normalized
        or "ratelimitexceeded" in normalized
        or "status 429" in normalized
        or "status: 429" in normalized
        or "code: 429" in normalized
    )


def _is_gemini_raw_error_dump_text(normalized: str) -> bool:
    """Return whether a Gemini line is part of a verbose JS/Gaxios dump."""
    return (
        "gaxioserror" in normalized
        or "gaxios-gaxios-error" in normalized
        or "cloudcode-pa.googleapis.com" in normalized
        or "cloudcodevscode" in normalized
        or "google-api-nodejs-client" in normalized
        or normalized.startswith("at gaxios.")
        or normalized.startswith("at async ")
        or normalized.startswith("config: {")
        or normalized.startswith("response: {")
        or normalized.startswith("headers: {")
        or normalized.startswith("request: {")
        or normalized.startswith("retryconfig: {")
        or normalized.startswith("body: '{\"model\"")
        or "\"systeminstruction\"" in normalized
        or "\"<session_context>" in normalized
    )


def _is_assistant_stream_json_line(line: str) -> bool:
    """Return whether a raw provider line is assistant JSON text."""
    stripped = str(line or "").strip()
    if not stripped.startswith("{"):
        return False
    try:
        payload = json.loads(stripped)
    except json.JSONDecodeError:
        return False
    return bool(_text_from_json_payload(payload))


def _provider_exit_error_message(
    kind: str,
    label: str,
    return_code: int,
    tail: str,
) -> str:
    """Return a concise provider failure message suitable for the GUI."""
    normalized_tail = tail.lower()
    if kind == "gemini" and (
        "please set an auth method" in normalized_tail
        or "gemini_api_key" in normalized_tail
    ):
        return (
            "Gemini is installed but not authenticated. Open Terminal and run "
            "`gemini` once to configure its auth method, or set `GEMINI_API_KEY` "
            "or configure `~/.gemini/settings.json`. Then retry Gemini in "
            "Equilipy."
        )
    if kind == "gemini" and _is_gemini_capacity_error_text(normalized_tail):
        return (
            "Gemini is installed and authenticated, but Google returned a "
            "temporary capacity/rate-limit error for the selected model. Retry "
            "Gemini in a moment, or switch the sidebar provider to Codex or "
            "Claude for this request."
        )
    detail = f". Last output: {tail}" if tail else ""
    return f"{label} exited with status {return_code}{detail}"


def build_prompt(
    context: AssistantContext,
    message: str,
    *,
    cwd: str | os.PathLike[str] | None = None,
) -> str:
    """Build the prompt sent to online/CLI providers."""
    return (
        f"{assistant_system_prompt()}\n\n"
        f"{build_claude_user_prompt(context, message, cwd=cwd)}"
    )


def assistant_system_prompt() -> str:
    """Return provider instructions for the Equilipy GUI action bridge."""
    return action_contract_prompt()


def build_claude_user_prompt(
    context: AssistantContext,
    message: str,
    *,
    cwd: str | os.PathLike[str] | None = None,
) -> str:
    """Build Claude's user prompt with the GUI state and request."""
    skills = agent_skills_prompt(context, message, cwd=cwd)
    parts = []
    if skills:
        parts.append(skills)
    parts.append(context.to_prompt_text())
    parts.append(f"User request:\n{message.strip()}")
    return "\n\n".join(parts)


def parse_stream_chunk(line: str) -> str:
    """Parse plain text, JSONL, or stream-json provider output into display text."""
    text = str(line or "")
    if not text.strip():
        return ""
    stripped = text.strip()
    if not stripped.startswith("{"):
        return text
    try:
        payload = json.loads(stripped)
    except json.JSONDecodeError:
        return text
    return _text_from_json_payload(payload)


def _text_from_json_payload(payload: object) -> str:
    if isinstance(payload, str):
        return payload
    if not isinstance(payload, dict):
        return ""
    structured = payload.get("structured_output")
    text = structured_response_to_text(structured)
    if text:
        return text
    if _is_non_assistant_json_payload(payload):
        return ""
    for key in (
        "text",
        "content",
        "delta",
        "message",
        "output",
        "value",
        "response",
        "item",
        "result",
    ):
        value = payload.get(key)
        text = _coerce_text(value)
        if text:
            return text
    if payload.get("type") in {"assistant", "assistant_message", "message"}:
        return _coerce_text(payload.get("data"))
    return ""


def _is_non_assistant_json_payload(payload: dict) -> bool:
    """Return whether a JSON stream item is clearly not assistant text."""
    role = str(
        payload.get("role")
        or payload.get("author")
        or payload.get("speaker")
        or ""
    ).strip().lower()
    if role in {"user", "system", "tool"}:
        return True
    event = str(payload.get("event") or "").strip().lower()
    if event in {"input", "prompt", "user", "system", "tool"}:
        return True
    payload_type = str(payload.get("type") or "").strip().lower()
    return any(
        token in payload_type
        for token in (
            "user",
            "system",
            "tool_call",
            "tool_result",
            "debug",
            "thought",
        )
    )


def _coerce_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        for key in ("text", "content", "delta", "message"):
            text = _coerce_text(value.get(key))
            if text:
                return text
        return ""
    if isinstance(value, list):
        parts = [_coerce_text(item) for item in value]
        return "".join(part for part in parts if part)
    return ""


def available_providers(
    *,
    path_lookup=shutil.which,
    extra_paths: dict[str, str | None] | None = None,
) -> list[ProviderInfo]:
    """Return provider availability, including known macOS app/Homebrew paths."""
    extra_paths = extra_paths or {}
    providers = [
        _provider_info(
            "codex",
            "Codex",
            "codex",
            (
                "exec",
                "--json",
                "--skip-git-repo-check",
                "--ephemeral",
                "--color",
                "never",
                "-",
            ),
            path_lookup=path_lookup,
            extra_path=extra_paths.get("codex"),
            known_paths=(
                "/Applications/Codex.app/Contents/Resources/codex",
                "/opt/homebrew/bin/codex",
                "/usr/local/bin/codex",
            ),
        ),
        _provider_info(
            "claude",
            "Claude",
            "claude",
            (
                "-p",
                "<prompt>",
                "--system-prompt",
                "<instructions>",
                "--output-format",
                "json",
                "--json-schema",
                "<schema>",
                "--tools",
                '""',
                "--no-session-persistence",
            ),
            path_lookup=path_lookup,
            extra_path=extra_paths.get("claude"),
            known_paths=(
                "/opt/homebrew/bin/claude",
                "/usr/local/bin/claude",
                "~/.local/bin/claude",
            ),
        ),
        _provider_info(
            "gemini",
            "Gemini",
            "gemini",
            ("-p", "<prompt>", "--output-format", "stream-json"),
            path_lookup=path_lookup,
            extra_path=extra_paths.get("gemini"),
            known_paths=(
                "/opt/homebrew/bin/gemini",
                "/usr/local/bin/gemini",
                "~/.local/bin/gemini",
            ),
        ),
    ]
    return providers


def create_cli_provider(info: ProviderInfo) -> CliProvider:
    """Create a CLI provider from available provider metadata."""
    if not info.executable:
        raise AssistantProviderError(f"{info.label} is not available.")
    return CliProvider(info.name, info.label, info.executable, info.name)


def cli_command(
    kind: str, executable: str, prompt: str
) -> tuple[list[str], str | None]:
    """Return a provider command and optional stdin payload."""
    if kind == "codex":
        return [
            executable,
            "exec",
            "--json",
            "--skip-git-repo-check",
            "--ephemeral",
            "--color",
            "never",
            "-",
        ], prompt
    if kind == "gemini":
        return [executable, "-p", prompt, "--output-format", "stream-json"], None
    if kind == "claude":
        return [
            executable,
            "-p",
            prompt,
            "--system-prompt",
            assistant_system_prompt(),
            "--output-format",
            "json",
            "--json-schema",
            action_response_schema_json(),
            "--tools",
            "",
            "--no-session-persistence",
        ], None
    raise AssistantProviderError(f"Unknown CLI provider: {kind}")


def cli_environment(kind: str) -> dict[str, str] | None:
    """Return provider-specific subprocess environment overrides."""
    if kind != "codex":
        return None
    env = os.environ.copy()
    codex_home = _prepare_equilipy_codex_home()
    env["CODEX_HOME"] = str(codex_home)
    tmp_dir = codex_home / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    env["TMPDIR"] = f"{tmp_dir}{os.sep}"
    return env


def _prepare_equilipy_codex_home() -> Path:
    """Prepare a writable Codex home for embedded GUI sidebar launches."""
    configured = os.environ.get("EQUILIPY_CODEX_HOME")
    if configured:
        codex_home = Path(configured).expanduser()
    else:
        codex_home = (
            Path.home()
            / "Library"
            / "Application Support"
            / "Equilipy"
            / "assistant"
            / "codex"
        )
    codex_home.mkdir(parents=True, exist_ok=True)

    user_home = Path(
        os.environ.get("CODEX_USER_HOME", Path.home() / ".codex")
    ).expanduser()
    if user_home.resolve() != codex_home.resolve():
        for name in (
            "auth.json",
            "config.toml",
            "AGENTS.md",
            "cloud-requirements-cache.json",
            "models_cache.json",
            "rules",
            "plugins",
        ):
            _bridge_codex_home_item(user_home / name, codex_home / name)
    return codex_home


def _bridge_codex_home_item(source: Path, destination: Path) -> None:
    """Expose read-mostly Codex config/auth files inside the writable home."""
    if not source.exists() or destination.exists() or destination.is_symlink():
        return
    try:
        destination.symlink_to(source, target_is_directory=source.is_dir())
    except OSError:
        try:
            if source.is_dir():
                shutil.copytree(source, destination)
            else:
                shutil.copy2(source, destination)
        except OSError:
            return


def _provider_info(
    name: str,
    label: str,
    binary: str,
    preview_args: Iterable[str],
    *,
    path_lookup,
    extra_path: str | None,
    known_paths: Iterable[str],
) -> ProviderInfo:
    executable = _find_executable(
        binary,
        path_lookup=path_lookup,
        extra_path=extra_path,
        known_paths=known_paths,
    )
    status = "available" if executable else "not_found"
    preview = (Path(executable).name if executable else binary, *tuple(preview_args))
    return ProviderInfo(name, label, executable, status, preview)


def _find_executable(
    binary: str,
    *,
    path_lookup,
    extra_path: str | None,
    known_paths: Iterable[str],
) -> str | None:
    candidates: list[str] = []
    if extra_path:
        candidates.append(extra_path)
    found = path_lookup(binary)
    if found:
        candidates.append(found)
    candidates.extend(known_paths)
    for candidate in candidates:
        expanded = os.path.expanduser(str(candidate))
        if os.path.isfile(expanded) and os.access(expanded, os.X_OK):
            return expanded
    return None
