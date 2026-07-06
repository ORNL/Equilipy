"""Launch or build the Equilipy documentation."""

from __future__ import annotations

import argparse
import subprocess
import sys
import webbrowser
from http.server import SimpleHTTPRequestHandler
from pathlib import Path
from socketserver import TCPServer

DOCS_URL = "https://ornl.github.io/Equilipy/"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Open the Equilipy documentation.")
    parser.add_argument(
        "--online",
        action="store_true",
        help="Open the published online documentation instead of local docs.",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Print the documentation location without opening a browser.",
    )
    parser.add_argument(
        "--serve",
        action="store_true",
        help="Serve the local HTML documentation with a small HTTP server.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8000,
        help="Port to use with --serve. Default: 8000.",
    )
    return parser


def _find_docs_source() -> Path | None:
    for parent in Path(__file__).resolve().parents:
        docs_source = parent / "docs"
        if (docs_source / "conf.py").exists():
            return docs_source
    return None


def _build_local_docs(docs_source: Path) -> Path | None:
    build_dir = docs_source / "_build" / "html"
    index = build_dir / "index.html"
    try:
        subprocess.run(
            [
                sys.executable,
                "-m",
                "sphinx",
                "-b",
                "html",
                str(docs_source),
                str(build_dir),
            ],
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        if index.exists():
            return index
        return None
    return index


def _serve_docs(index: Path, port: int) -> int:
    class DocsHandler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs) -> None:
            super().__init__(*args, directory=str(index.parent), **kwargs)

    with TCPServer(("", port), DocsHandler) as server:
        url = f"http://localhost:{port}"
        print(f"Serving Equilipy docs at {url}")
        webbrowser.open(url)
        server.serve_forever()
    return 0


def main(argv: list[str] | None = None) -> int:
    """Open, build, or serve the Equilipy documentation."""
    args = _build_parser().parse_args(argv)

    if args.online:
        print(DOCS_URL)
        if not args.no_browser:
            webbrowser.open(DOCS_URL)
        return 0

    docs_source = _find_docs_source()
    if docs_source is None:
        print(DOCS_URL)
        if not args.no_browser:
            webbrowser.open(DOCS_URL)
        return 0

    index = _build_local_docs(docs_source)
    if index is None:
        print(
            "Could not build local docs. Install the docs dependencies with "
            "`pip install -r docs/requirements.txt`, or use "
            "`python -m equilipy.docs --online`.",
            file=sys.stderr,
        )
        return 1

    if args.serve:
        return _serve_docs(index, args.port)

    location = index.resolve().as_uri()
    print(location)
    if not args.no_browser:
        webbrowser.open(location)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
