"""Repo-local assistant skill discovery for provider prompts."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .context import AssistantContext

MAX_SKILL_PROMPT_CHARS = 14000


@dataclass(frozen=True)
class AgentSkill:
    """One repo-local assistant skill."""

    name: str
    description: str
    path: Path
    body: str


def agent_skills_prompt(
    context: AssistantContext,
    message: str,
    *,
    cwd: str | os.PathLike[str] | None = None,
    max_chars: int = MAX_SKILL_PROMPT_CHARS,
) -> str:
    """Return applicable repo-local skill instructions for a provider prompt."""
    skills = applicable_agent_skills(context, message, cwd=cwd)
    if not skills:
        return ""
    sections = [
        "Equilipy agent skills loaded from agent_skills/.",
        "Follow these instructions before answering the user request.",
    ]
    remaining = max_chars - len("\n".join(sections))
    for skill in skills:
        if remaining <= 0:
            break
        text = _skill_prompt_section(skill, remaining)
        if not text:
            continue
        sections.append(text)
        remaining -= len(text)
    return "\n\n".join(sections).strip()


def applicable_agent_skills(
    context: AssistantContext,
    message: str,
    *,
    cwd: str | os.PathLike[str] | None = None,
) -> list[AgentSkill]:
    """Return repo-local skills that should guide this request."""
    skills = discover_agent_skills(cwd=cwd)
    if not skills:
        return []
    text = _selection_text(context, message)
    selected = [skill for skill in skills if _skill_matches(skill, text)]
    if selected:
        return selected
    return []


def discover_agent_skills(
    *,
    cwd: str | os.PathLike[str] | None = None,
) -> list[AgentSkill]:
    """Discover skill folders from known agent_skills roots."""
    skills: list[AgentSkill] = []
    seen: set[Path] = set()
    seen_names: set[str] = set()
    for root in _candidate_skill_roots(cwd=cwd):
        if not root.is_dir():
            continue
        for skill_md in sorted(
            root.glob("*/SKILL.md"), key=lambda item: item.as_posix()
        ):
            path = skill_md.parent.resolve()
            if path in seen:
                continue
            seen.add(path)
            skill = _load_agent_skill(skill_md)
            if skill is not None:
                if skill.name in seen_names:
                    continue
                seen_names.add(skill.name)
                skills.append(skill)
    return skills


def _candidate_skill_roots(
    *,
    cwd: str | os.PathLike[str] | None = None,
) -> Iterable[Path]:
    configured = os.environ.get("EQUILIPY_AGENT_SKILLS_DIR", "")
    for value in configured.split(os.pathsep):
        if value.strip():
            yield Path(value).expanduser()
    if cwd is not None:
        current = Path(cwd).expanduser().resolve()
        yield current / "agent_skills"
        for parent in current.parents:
            yield parent / "agent_skills"
    yield Path(__file__).resolve().parents[3] / "agent_skills"
    yield Path.cwd() / "agent_skills"


def _load_agent_skill(path: Path) -> AgentSkill | None:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    metadata, body = _split_frontmatter(text)
    name = metadata.get("name") or path.parent.name
    description = metadata.get("description") or ""
    return AgentSkill(
        name=name,
        description=description,
        path=path.parent,
        body=body.strip(),
    )


def _split_frontmatter(text: str) -> tuple[dict[str, str], str]:
    if not text.startswith("---\n"):
        return {}, text
    end = text.find("\n---\n", 4)
    if end < 0:
        return {}, text
    metadata: dict[str, str] = {}
    for line in text[4:end].splitlines():
        key, separator, value = line.partition(":")
        if separator:
            metadata[key.strip()] = value.strip().strip("\"'")
    return metadata, text[end + len("\n---\n") :]


def _selection_text(context: AssistantContext, message: str) -> str:
    parts = [
        context.workspace,
        context.database_name,
        context.selected_object,
        context.active_calculation_item,
        context.status_message,
        message,
    ]
    parts.extend(context.loaded_databases)
    parts.extend(context.calculation_sessions)
    parts.extend(context.diagnostics)
    parts.extend(context.extra.values())
    return " ".join(part for part in parts if part).lower()


def _skill_matches(skill: AgentSkill, text: str) -> bool:
    if not text:
        return False
    if skill.name == "equilipy-calculation-assist":
        return _calculation_skill_matches(text)
    haystack = f"{skill.name} {skill.description}".lower()
    tokens = set(re.findall(r"[a-z0-9][a-z0-9_.%-]*", haystack))
    informative = {
        token
        for token in tokens
        if (len(token) >= 4 or token == "tdb")
        and token
        not in {
            "assistant",
            "assist",
            "database",
            "databases",
            "equilipy",
            "guide",
            "actions",
            "action",
            "through",
            "workflow",
            "workflows",
            "when",
            "with",
            "that",
            "user",
            "asks",
        }
    }
    matches = {token for token in informative if token in text}
    return len(matches) >= 2


def _calculation_skill_matches(text: str) -> bool:
    calculation_terms = {
        "batch",
        "calculate",
        "calculation",
        "composition",
        "condition",
        "equilibrium",
        "grid",
        "module",
        "pressure",
        "scheil",
        "session",
        "solidification",
        "temperature",
    }
    action_terms = {
        "add",
        "configure",
        "create",
        "load",
        "populate",
        "run",
        "select",
        "set",
        "setup",
    }
    has_calculation_term = any(term in text for term in calculation_terms)
    has_action_term = any(term in text for term in action_terms)
    return has_calculation_term and has_action_term


def _skill_prompt_section(skill: AgentSkill, max_chars: int) -> str:
    body = skill.body
    references = _load_references(skill.path)
    section = (
        f"Skill: {skill.name}\n"
        f"Description: {skill.description}\n\n"
        f"{body}"
    )
    if references:
        section = f"{section}\n\nReferences:\n{references}"
    if len(section) <= max_chars:
        return section
    if max_chars < 200:
        return ""
    return section[: max_chars - 35].rstrip() + "\n[skill content truncated]"


def _load_references(skill_path: Path) -> str:
    reference_dir = skill_path / "references"
    if not reference_dir.is_dir():
        return ""
    sections: list[str] = []
    for path in sorted(reference_dir.glob("*.md"), key=lambda item: item.name.lower()):
        try:
            text = path.read_text(encoding="utf-8").strip()
        except OSError:
            continue
        if text:
            sections.append(f"Reference: {path.name}\n{text}")
    return "\n\n".join(sections)
