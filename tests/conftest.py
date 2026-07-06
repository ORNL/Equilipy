"""Shared fixtures for the lightweight public test suite.

Author:
    S.Y. Kwon
Date:
    06/30/2026
Tests:
    1. database_dir and path fixtures: Resolve public databases from the repo
       ``database`` folder and fail clearly when a required fixture is missing.
    2. parsed database fixtures: Share lightweight parsed FS73, Liu Al-Fe-Si,
       CAS FS83, and CAS MQM FS73 databases across the public tests.
"""

from pathlib import Path

import pytest

import equilipy as eq

REPO_ROOT = Path(__file__).resolve().parents[1]
DATABASE_DIR = REPO_ROOT / "database"


def _database_path(database_dir: Path, file_name: str) -> Path:
    path = database_dir / file_name
    if not path.is_file():
        raise FileNotFoundError(f"Expected test database in database/: {path}")
    return path


@pytest.fixture(scope="session")
def database_dir() -> Path:
    """Return the repository database fixture directory."""
    return DATABASE_DIR


@pytest.fixture(scope="session")
def alcumgsi_fs73_database_path(database_dir: Path) -> Path:
    """Return the small Al-Cu-Mg-Si FactSage 7.3 DAT database path."""
    return _database_path(database_dir, "AlCuMgSi_ORNL_FS73.dat")


@pytest.fixture(scope="session")
def alcumgsi_fs83_database_path(database_dir: Path) -> Path:
    """Return the small Al-Cu-Mg-Si FactSage 8.3 DAT database path."""
    return _database_path(database_dir, "AlCuMgSi_ORNL_FS83.dat")


@pytest.fixture(scope="session")
def liu_alfesi_dat_database_path(database_dir: Path) -> Path:
    """Return the Liu Al-Fe-Si ChemSage DAT database path."""
    return _database_path(database_dir, "AlFeSi_99Liu.dat")


@pytest.fixture(scope="session")
def liu_alfesi_tdb_database_path(database_dir: Path) -> Path:
    """Return the Liu Al-Fe-Si TDB database path."""
    return _database_path(database_dir, "AlFeSi_99Liu.tdb")


@pytest.fixture(scope="session")
def cas_fs83_database_path(database_dir: Path) -> Path:
    """Return the CAS FactSage 8.3 DAT database path."""
    return _database_path(database_dir, "CAS_FS83.dat")


@pytest.fixture(scope="session")
def cas_mqm_fs73_database_path(database_dir: Path) -> Path:
    """Return the CAS MQM FactSage 7.3 DAT database path."""
    return _database_path(database_dir, "CAS_MQM_FS73.dat")


@pytest.fixture(scope="module")
def alcumgsi_fs73_database(alcumgsi_fs73_database_path: Path) -> dict:
    """Return parsed Al-Cu-Mg-Si FS73 database data."""
    return eq.read_dat(str(alcumgsi_fs73_database_path), factsage_8_plus=False)


@pytest.fixture(scope="module")
def liu_alfesi_dat_database(liu_alfesi_dat_database_path: Path) -> dict:
    """Return parsed Liu Al-Fe-Si DAT database data."""
    return eq.read_dat(str(liu_alfesi_dat_database_path))


@pytest.fixture(scope="module")
def liu_alfesi_tdb_database(liu_alfesi_tdb_database_path: Path) -> dict:
    """Return parsed Liu Al-Fe-Si TDB database data."""
    return eq.read_tdb(liu_alfesi_tdb_database_path)


@pytest.fixture(scope="module")
def cas_fs83_database(cas_fs83_database_path: Path) -> dict:
    """Return parsed CAS FS83 database data."""
    return eq.read_dat(str(cas_fs83_database_path), factsage_8_plus=True)


@pytest.fixture(scope="module")
def cas_mqm_fs73_database(cas_mqm_fs73_database_path: Path) -> dict:
    """Return parsed CAS MQM FS73 database data."""
    return eq.read_dat(str(cas_mqm_fs73_database_path), factsage_8_plus=False)
