"""Sphinx configuration for Equilipy documentation."""

from __future__ import annotations

from datetime import datetime
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path
import re

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    tomllib = None

# ROOT is used to read the version from pyproject.toml. It is deliberately
# NOT inserted into sys.path: autodoc must import the installed equilipy
# (which has a loadable equilifort extension), not the source checkout —
# the in-repo package would shadow it and fail to import on machines whose
# platform/Python does not match the locally built .so.
ROOT = Path(__file__).resolve().parents[1]

project = "Equilipy"
author = "Sunyong Kwon"
copyright = f"2024-{datetime.now().year}, U.S. Department of Energy"


def _read_project_version() -> str:
    """Read the package version from pyproject.toml or installed metadata."""
    pyproject_path = ROOT / "pyproject.toml"
    if pyproject_path.exists():
        if tomllib is not None:
            with pyproject_path.open("rb") as pyproject_file:
                pyproject = tomllib.load(pyproject_file)
            project_table = pyproject.get("project", {})
            if project_table.get("version"):
                return project_table["version"]

        pyproject_text = pyproject_path.read_text(encoding="utf-8")
        match = re.search(r'(?m)^version\s*=\s*["\']([^"\']+)["\']', pyproject_text)
        if match:
            return match.group(1)

    try:
        return package_version("equilipy")
    except PackageNotFoundError:
        return "0+unknown"


release = _read_project_version()
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinxcontrib.mermaid",
    "sphinx_copybutton",
]

# Copy button: strip shell/Python prompts so pasted code runs as-is.
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "_site",
    "Thumbs.db",
    ".DS_Store",
    "Gemfile",
    "Gemfile.lock",
    # Developer documentation stays in the repo but is not published yet.
    # Restore the hidden toctree in index.md when re-enabling.
    "architecture.md",
    "developer/**",
]

myst_heading_anchors = 3
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

html_theme = "sphinx_rtd_theme"
html_logo = "logo/Equilipy.svg"
html_favicon = "logo/Equilipy.svg"
html_title = "Equilipy"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
# Wraps content images in links that open the full-size file in a new tab.
html_js_files = ["custom.js"]
html_extra_path = [".nojekyll"]
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "sticky_navigation": True,
}
html_context = {
    "display_github": True,
    "github_user": "ORNL",
    "github_repo": "Equilipy",
    "github_version": "main",
    "conf_py_path": "/docs/",
}


# ---------------------------------------------------------------------------
# Custom Pygments lexer for Thermo-Calc TDB snippets (```tdb code fences).
# Pygments has no built-in "tdb" lexer; without this every TDB block warns
# with misc.highlighting_failure and renders unhighlighted.
# ---------------------------------------------------------------------------
from pygments.lexer import RegexLexer  # noqa: E402
from pygments.token import (  # noqa: E402
    Comment,
    Keyword,
    Name,
    Number,
    Punctuation,
    Text,
    Whitespace,
)


class TdbLexer(RegexLexer):
    """Minimal highlighter for Thermo-Calc database (TDB) files."""

    name = "TDB"
    aliases = ["tdb"]
    filenames = ["*.tdb", "*.TDB"]

    tokens = {
        "root": [
            (r"\$[^\n]*", Comment.Single),
            (
                r"(?i)^\s*(ELEMENT|SPECIES|FUNCTION|PHASE|CONSTITUENT|CONST|"
                r"PARAMETER|PAR|TYPE_DEFINITION|TYPE_DEF|"
                r"DEFINE_SYSTEM_DEFAULT|DEFAULT_COMMAND|DEF_SYS_ELEMENT|"
                r"ASSESSED_SYSTEMS|LIST_OF_REFERENCES|REFERENCE_FILE|"
                r"DATABASE_INFO|VERSION_DATE|TEMPERATURE_LIMITS)\b",
                Keyword,
            ),
            (
                r"(?i)\b(GES|AMEND_PHASE_DESCRIPTION|MAGNETIC|DIS_PART|"
                r"DISORDER_PART|SEQ)\b",
                Keyword.Type,
            ),
            (r"[+-]?(\d+\.\d*|\.\d+|\d+)([EeDd][+-]?\d+)?", Number),
            (r"[!;]", Punctuation),
            (r"[:,#%@*<>()=+\-/]", Punctuation),
            (r"[A-Za-z_][A-Za-z0-9_.'\-]*", Name),
            (r"\s+", Whitespace),
            (r".", Text),
        ],
    }


def setup(app):
    """Register project-specific Sphinx customizations."""
    app.add_lexer("tdb", TdbLexer)
