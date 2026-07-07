"""Tests for the --setup-linux-libs helper logic."""

from equilipy.gui.linux_libs import _XCB_SONAMES, ensure_xcb_libs


def test_already_bundled_lib_is_left_alone(tmp_path):
    """A library already present in Qt/lib is not overwritten."""
    qt_lib = tmp_path / "qtlib"
    qt_lib.mkdir()
    bundled = qt_lib / _XCB_SONAMES[0]
    bundled.write_bytes(b"bundled")

    actions, missing = ensure_xcb_libs(qt_lib, set(), [])

    assert bundled.read_bytes() == b"bundled"
    assert any("already" in action for action in actions)
    assert _XCB_SONAMES[0] not in missing


def test_system_provided_lib_is_not_copied(tmp_path):
    """Libraries the system loader resolves are left to the system."""
    qt_lib = tmp_path / "qtlib"
    qt_lib.mkdir()

    actions, missing = ensure_xcb_libs(qt_lib, set(_XCB_SONAMES), [])

    assert missing == []
    assert not any(qt_lib.iterdir())
    assert all("system" in action for action in actions)


def test_missing_lib_is_copied_from_candidate_dir(tmp_path):
    """A library absent from Qt/lib and the system is copied in."""
    qt_lib = tmp_path / "qtlib"
    qt_lib.mkdir()
    conda_lib = tmp_path / "conda" / "lib"
    conda_lib.mkdir(parents=True)
    source = conda_lib / _XCB_SONAMES[0]
    source.write_bytes(b"conda copy")

    system = set(_XCB_SONAMES) - {_XCB_SONAMES[0]}
    actions, missing = ensure_xcb_libs(qt_lib, system, [conda_lib])

    assert missing == []
    assert (qt_lib / _XCB_SONAMES[0]).read_bytes() == b"conda copy"
    assert any("copied" in action for action in actions)


def test_symlinked_source_is_copied_as_regular_file(tmp_path):
    """Symlinked sources are materialized as regular files."""
    qt_lib = tmp_path / "qtlib"
    qt_lib.mkdir()
    conda_lib = tmp_path / "conda" / "lib"
    conda_lib.mkdir(parents=True)
    real = conda_lib / (_XCB_SONAMES[0] + ".0.0")
    real.write_bytes(b"real")
    (conda_lib / _XCB_SONAMES[0]).symlink_to(real)

    system = set(_XCB_SONAMES) - {_XCB_SONAMES[0]}
    _actions, missing = ensure_xcb_libs(qt_lib, system, [conda_lib])

    copied = qt_lib / _XCB_SONAMES[0]
    assert missing == []
    assert copied.is_file() and not copied.is_symlink()
    assert copied.read_bytes() == b"real"


def test_unlocatable_lib_is_reported_missing(tmp_path):
    """A library found nowhere is reported as missing."""
    qt_lib = tmp_path / "qtlib"
    qt_lib.mkdir()

    system = set(_XCB_SONAMES) - {_XCB_SONAMES[0]}
    _actions, missing = ensure_xcb_libs(qt_lib, system, [])

    assert missing == [_XCB_SONAMES[0]]
