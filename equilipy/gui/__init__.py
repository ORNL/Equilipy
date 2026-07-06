"""PySide6 GUI package for Equilipy database editing."""

from __future__ import annotations

__all__ = ["main"]


def main(argv: list[str] | None = None) -> int:
    """Launch the Equilipy GUI."""
    from .app import main as run_app

    return run_app(argv)
