"""Primary command-line entry point for Equilipy."""

from __future__ import annotations

from multiprocessing import freeze_support


def main() -> int:
    """Launch the GUI; shared by ``python -m equilipy`` and the console script."""
    freeze_support()
    from .gui.app import main as gui_main

    return gui_main()


if __name__ == "__main__":
    raise SystemExit(main())
