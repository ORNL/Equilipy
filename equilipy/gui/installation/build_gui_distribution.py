"""Build one local Equilipy GUI distribution artifact for the current platform."""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
APP_NAME = "Equilipy"
DIST_GUI = PROJECT_ROOT / "dist-gui"
PYINSTALLER_DIST = PROJECT_ROOT / "pyinstaller-dist"
PYINSTALLER_BUILD = PROJECT_ROOT / "pyinstaller-build"
DMG_ROOT = PROJECT_ROOT / "dmg-root"
ENTRY_POINT = PROJECT_ROOT / "equilipy" / "gui" / "_pyinstaller_entry.py"
ICON_DIR = PROJECT_ROOT / "equilipy" / "gui" / "icons"
ICON_ICO = ICON_DIR / "equilipy.ico"
ICON_ICNS = ICON_DIR / "equilipy.icns"
MACOS_ENTITLEMENTS = (
    PROJECT_ROOT / "equilipy" / "gui" / "installation" / "macos_entitlements.plist"
)
SMOKE_DATABASE = PROJECT_ROOT / "database" / "AlFeSi_99Liu.dat"
MACOS_CODESIGN_IDENTITY_ENV = "EQUILIPY_MACOS_CODESIGN_IDENTITY"
MACOS_NOTARY_PROFILE_ENV = "EQUILIPY_MACOS_NOTARY_PROFILE"


def _default_label(system_name: str) -> str:
    """Return a platform label for the generated artifact name."""
    machine = platform.machine().lower() or "unknown"
    if system_name == "Darwin":
        return f"macos-{machine}"
    if system_name == "Windows":
        return f"windows-{machine}"
    if system_name == "Linux":
        return f"linux-{machine}"
    return f"{system_name.lower()}-{machine}"


def _run(command: list[str]) -> None:
    """Run a subprocess command from the project root."""
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, check=True)


def _require_tool(name: str) -> None:
    """Fail early when a required release-build command is unavailable."""
    if shutil.which(name) is None:
        raise SystemExit(f"{name!r} is required for macOS signed release builds.")


def _remove_path(path: Path) -> None:
    """Remove one file, symlink, or directory if it exists."""
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def _clean_build_dirs() -> None:
    """Remove generated GUI distribution build directories."""
    for path in (DIST_GUI, PYINSTALLER_DIST, PYINSTALLER_BUILD, DMG_ROOT):
        _remove_path(path)


def _clean_temporary_build_dirs() -> None:
    """Remove temporary GUI distribution build directories."""
    for path in (PYINSTALLER_DIST, PYINSTALLER_BUILD, DMG_ROOT):
        _remove_path(path)


def _pyinstaller_hidden_imports(system_name: str) -> list[str]:
    """Return hidden imports needed by the GUI bundle on this platform."""
    loky_backend = (
        "joblib.externals.loky.backend.popen_loky_win32"
        if system_name == "Windows"
        else "joblib.externals.loky.backend.popen_loky_posix"
    )
    return [
        "equilipy.equilifort",
        "equilipy._parallel",
        loky_backend,
        "joblib.externals.loky.backend.resource_tracker",
    ]


def _pyinstaller_data_specs() -> list[str]:
    """Return data files that must be bundled for frozen-app smoke tests."""
    candidates = [
        SMOKE_DATABASE,
        Path.cwd() / "database" / "AlFeSi_99Liu.dat",
    ]
    for candidate in candidates:
        if candidate.exists():
            return [f"{candidate}{os.pathsep}database"]
    return []


def _pyinstaller_command(*, onefile: bool, icon_path: Path | None = None) -> list[str]:
    """Return the PyInstaller command for the current platform."""
    command = [
        sys.executable,
        "-m",
        "PyInstaller",
        "--noconfirm",
        "--clean",
        "--windowed",
        "--name",
        APP_NAME,
        "--distpath",
        str(PYINSTALLER_DIST),
        "--workpath",
        str(PYINSTALLER_BUILD),
        "--specpath",
        str(PYINSTALLER_BUILD),
        "--collect-data",
        "equilipy.gui",
        str(ENTRY_POINT),
    ]
    for module_name in _pyinstaller_hidden_imports(platform.system()):
        command[-1:-1] = ["--hidden-import", module_name]
    for data_spec in _pyinstaller_data_specs():
        command[-1:-1] = ["--add-data", data_spec]
    if onefile:
        command.insert(5, "--onefile")
    if icon_path is not None:
        command[-1:-1] = ["--icon", str(icon_path)]
    return command


def _pyinstaller_artifact_path(system_name: str) -> Path:
    """Return the platform-specific PyInstaller build product path."""
    if system_name == "Darwin":
        return PYINSTALLER_DIST / f"{APP_NAME}.app"
    if system_name == "Windows":
        return PYINSTALLER_DIST / f"{APP_NAME}.exe"
    return PYINSTALLER_DIST / APP_NAME


def _package_linux(label: str) -> Path:
    """Copy the Linux one-file executable into dist-gui."""
    source = PYINSTALLER_DIST / APP_NAME
    target = DIST_GUI / f"Equilipy-{label}"
    shutil.copy2(source, target)
    target.chmod(0o755)
    return target


def _package_windows(label: str) -> Path:
    """Copy the Windows one-file executable into dist-gui."""
    source = PYINSTALLER_DIST / f"{APP_NAME}.exe"
    target = DIST_GUI / f"Equilipy-{label}.exe"
    shutil.copy2(source, target)
    return target


def _package_macos_dmg(label: str) -> Path:
    """Wrap the macOS app bundle in a single distributable dmg file."""
    settings = (
        PROJECT_ROOT
        / "equilipy"
        / "gui"
        / "installation"
        / "dmgbuild_settings.py"
    )
    try:
        import dmgbuild  # noqa: F401
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "dmgbuild is required for the styled macOS installer. "
            "Install it with `python -m pip install dmgbuild`."
        ) from exc

    target = DIST_GUI / f"Equilipy-{label}.dmg"

    _run(
        [
            sys.executable,
            "-m",
            "dmgbuild",
            "-s",
            str(settings),
            APP_NAME,
            str(target),
        ]
    )
    return target


def _macos_codesign_identity(identity: str | None = None) -> str:
    """Return the Developer ID identity used to sign macOS release artifacts."""
    value = identity or os.environ.get(MACOS_CODESIGN_IDENTITY_ENV, "")
    value = value.strip()
    if not value:
        raise SystemExit(
            "macOS signing requires a Developer ID Application identity. "
            f"Pass --codesign-identity or set {MACOS_CODESIGN_IDENTITY_ENV}."
        )
    return value


def _macos_notary_profile(profile: str | None = None) -> str:
    """Return the notarytool keychain profile used for notarization."""
    value = profile or os.environ.get(MACOS_NOTARY_PROFILE_ENV, "")
    value = value.strip()
    if not value:
        raise SystemExit(
            "macOS notarization requires an xcrun notarytool keychain profile. "
            f"Pass --notary-profile or set {MACOS_NOTARY_PROFILE_ENV}."
        )
    return value


def _codesign_macos_app(app_path: Path, identity: str) -> None:
    """Sign the PyInstaller app bundle for Gatekeeper notarization."""
    _require_tool("codesign")
    command = [
        "codesign",
        "--force",
        "--deep",
        "--options",
        "runtime",
        "--timestamp",
        "--entitlements",
        str(MACOS_ENTITLEMENTS),
        "--sign",
        identity,
        str(app_path),
    ]
    _run(command)
    _run(
        [
            "codesign",
            "--verify",
            "--deep",
            "--strict",
            "--verbose=2",
            str(app_path),
        ]
    )


def _codesign_macos_dmg(dmg_path: Path, identity: str) -> None:
    """Sign the distributable DMG before notarization."""
    _require_tool("codesign")
    _run(
        [
            "codesign",
            "--force",
            "--timestamp",
            "--sign",
            identity,
            str(dmg_path),
        ]
    )
    _run(["codesign", "--verify", "--verbose=2", str(dmg_path)])


def _submit_macos_notarization(path: Path, notary_profile: str) -> None:
    """Submit one macOS artifact with Apple's notary service."""
    _require_tool("xcrun")
    _run(
        [
            "xcrun",
            "notarytool",
            "submit",
            str(path),
            "--keychain-profile",
            notary_profile,
            "--wait",
        ]
    )


def _staple_macos_artifact(path: Path) -> None:
    """Staple and validate the notarization ticket on one macOS artifact."""
    _require_tool("xcrun")
    _run(["xcrun", "stapler", "staple", str(path)])
    _run(["xcrun", "stapler", "validate", str(path)])


def _notarize_macos_artifact(path: Path, notary_profile: str) -> None:
    """Submit, staple, and validate one macOS artifact with Apple's notary service."""
    _submit_macos_notarization(path, notary_profile)
    _staple_macos_artifact(path)


def _notarize_macos_app(app_path: Path, notary_profile: str) -> None:
    """Notarize and staple the app bundle before it is wrapped in the DMG."""
    _require_tool("ditto")
    _require_tool("spctl")
    archive = PYINSTALLER_BUILD / f"{APP_NAME}-notary.zip"
    _run(["ditto", "-c", "-k", "--keepParent", str(app_path), str(archive)])
    _submit_macos_notarization(archive, notary_profile)
    _staple_macos_artifact(app_path)
    _run(["spctl", "--assess", "--type", "execute", "--verbose=2", str(app_path)])


def _assess_macos_dmg(dmg_path: Path) -> None:
    """Verify the final DMG passes Gatekeeper assessment."""
    _require_tool("spctl")
    _run(["spctl", "--assess", "--type", "open", "--verbose=2", str(dmg_path)])


def build_distribution(
    label: str | None = None,
    *,
    clean: bool = True,
    app_only: bool = False,
    sign: bool = False,
    notarize: bool = False,
    codesign_identity: str | None = None,
    notary_profile: str | None = None,
) -> Path:
    """Build one GUI distribution artifact for the current platform."""
    system_name = platform.system()
    artifact_label = label or _default_label(system_name)

    if clean:
        _clean_build_dirs()
    DIST_GUI.mkdir(parents=True, exist_ok=True)

    is_macos = system_name == "Darwin"
    sign_macos = is_macos and (sign or notarize or codesign_identity is not None)
    notarize_macos = is_macos and notarize
    macos_identity = (
        _macos_codesign_identity(codesign_identity) if sign_macos else None
    )
    macos_notary_profile = (
        _macos_notary_profile(notary_profile) if notarize_macos else None
    )
    icon_path = None
    if system_name == "Darwin":
        icon_path = ICON_ICNS
    elif system_name == "Windows":
        icon_path = ICON_ICO
    _run(_pyinstaller_command(onefile=not is_macos, icon_path=icon_path))
    pyinstaller_artifact = _pyinstaller_artifact_path(system_name)
    if is_macos and macos_identity is not None:
        _codesign_macos_app(pyinstaller_artifact, macos_identity)
    if is_macos and macos_notary_profile is not None:
        _notarize_macos_app(pyinstaller_artifact, macos_notary_profile)
    if app_only:
        artifact = pyinstaller_artifact
        print(f"Created {artifact}", flush=True)
        return artifact

    try:
        if system_name == "Darwin":
            artifact = _package_macos_dmg(artifact_label)
            if macos_identity is not None:
                _codesign_macos_dmg(artifact, macos_identity)
            if macos_notary_profile is not None:
                _notarize_macos_artifact(artifact, macos_notary_profile)
                _assess_macos_dmg(artifact)
        elif system_name == "Windows":
            artifact = _package_windows(artifact_label)
        elif system_name == "Linux":
            artifact = _package_linux(artifact_label)
        else:
            raise SystemExit(f"Unsupported GUI distribution platform: {system_name}")
    except Exception:
        raise
    else:
        _clean_temporary_build_dirs()

    print(f"Created {artifact}", flush=True)
    return artifact


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Build one Equilipy GUI distribution artifact."
    )
    parser.add_argument(
        "--label",
        help="Artifact label used in the output filename.",
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Reuse existing PyInstaller build directories.",
    )
    parser.add_argument(
        "--app-only",
        action="store_true",
        help="Build only the PyInstaller app/executable, without packaging it.",
    )
    parser.add_argument(
        "--sign",
        action="store_true",
        help=(
            "Sign the macOS app/DMG with Developer ID. Also enabled by "
            "--notarize or --codesign-identity."
        ),
    )
    parser.add_argument(
        "--notarize",
        action="store_true",
        help=(
            "Submit and staple the macOS app and DMG using xcrun notarytool. "
            "Requires a Developer ID signature."
        ),
    )
    parser.add_argument(
        "--codesign-identity",
        help=(
            "Developer ID Application identity for macOS signing. Defaults to "
            f"${MACOS_CODESIGN_IDENTITY_ENV}."
        ),
    )
    parser.add_argument(
        "--notary-profile",
        help=(
            "xcrun notarytool keychain profile for macOS notarization. Defaults "
            f"to ${MACOS_NOTARY_PROFILE_ENV}."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the GUI distribution build command."""
    args = _build_parser().parse_args(argv)
    build_distribution(
        args.label,
        clean=not args.no_clean,
        app_only=args.app_only,
        sign=args.sign,
        notarize=args.notarize,
        codesign_identity=args.codesign_identity,
        notary_profile=args.notary_profile,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
