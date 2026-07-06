"""Asset helpers for the PySide6 GUI."""

from __future__ import annotations

import platform
from pathlib import Path


def equilipy_logo_path() -> Path:
    """Return the packaged Equilipy GUI logo path."""
    packaged_logo = Path(__file__).resolve().parent / "icons" / "equilipy.svg"
    if packaged_logo.exists():
        return packaged_logo
    return Path(__file__).resolve().parents[2] / "docs" / "logo" / "Equilipy.svg"


def equilipy_app_icon_path() -> Path:
    """Return the packaged native app/window icon path."""
    icon_dir = Path(__file__).resolve().parent / "icons"
    system_name = platform.system()
    if system_name == "Darwin":
        preferred = icon_dir / "equilipy.icns"
    elif system_name == "Windows":
        preferred = icon_dir / "equilipy.ico"
    else:
        preferred = icon_dir / "equilipy.png"
    if preferred.exists():
        return preferred
    fallback = icon_dir / "equilipy.png"
    if fallback.exists():
        return fallback
    return equilipy_logo_path()


def ornl_logo_path() -> Path:
    """Return the packaged ORNL wordmark path."""
    return Path(__file__).resolve().parent / "icons" / "ornl.svg"


def math_font_path() -> Path:
    """Return the packaged math font path."""
    return Path(__file__).resolve().parent / "fonts" / "STIXTwoMath.otf"


def calculation_icon_path(name: str) -> Path:
    """Return the vector icon path for a calculation module."""
    return Path(__file__).resolve().parent / "icons" / f"{name}.svg"


def status_icon_path(name: str) -> Path:
    """Return the vector status icon path."""
    return Path(__file__).resolve().parent / "icons" / f"status_{name}.svg"
