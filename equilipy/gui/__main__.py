"""Command-line entry point for ``python -m equilipy.gui``."""

from __future__ import annotations

from multiprocessing import freeze_support

if __name__ == "__main__":
    freeze_support()
    from .app import main

    raise SystemExit(main())
