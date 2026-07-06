"""Right-sidebar assistant support for the Equilipy GUI."""

from __future__ import annotations

from .context import AssistantContext, build_assistant_context
from .panel import AssistantPanel
from .providers import (
    AssistantProvider,
    AssistantProviderError,
    CliProvider,
    ProviderInfo,
    available_providers,
    build_prompt,
    parse_stream_chunk,
)

__all__ = [
    "AssistantContext",
    "AssistantPanel",
    "AssistantProvider",
    "AssistantProviderError",
    "CliProvider",
    "ProviderInfo",
    "available_providers",
    "build_assistant_context",
    "build_prompt",
    "parse_stream_chunk",
]
