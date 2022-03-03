"""Test versiondb.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import os

from argparse import Namespace
from pathlib import Path
from typing import List, NamedTuple, Tuple

import pandas as pd
import pytest
import unittest
import shutil
import filecmp

from pandas.util.testing import assert_frame_equal

from pyani import versiondb, pyani_files, pyani_tools


# Test get_version()
# Test case 0: no executable location is specified
def test_get_version_nonetype():
    """Test behaviour when no location for the executable is given."""
    test_file_0 = None

    assert versiondb.get_version(test_file_0) == f"{test_file_0} is not found in $PATH"


# Test case 1: there is no executable
def test_get_version_no_exe(executable_missing, monkeypatch):
    """Test behaviour when there is no file at the specified executable location."""
    test_file_1 = Path("/non/existent/alembic")
    assert versiondb.get_version(test_file_1) == f"No alembic at {test_file_1}"


# Test case 2: there is a file, but it is not executable
def test_get_version_exe_not_executable(executable_not_executable, monkeypatch):
    """Test behaviour when the file at the executable location is not executable."""
    test_file_2 = Path("/non/executable/alembic")
    assert (
        versiondb.get_version(test_file_2)
        == f"alembic exists at {test_file_2} but not executable"
    )


# Test case 3: there is an executable file, but the version can't be retrieved
def test_get_version_exe_no_version(executable_without_version, monkeypatch):
    """Test behaviour when the version for the executable can not be retrieved."""
    test_file_3 = Path("/missing/version/alembic")
    assert (
        versiondb.get_version(test_file_3)
        == f"alembic exists at {test_file_3} but could not retrieve version"
    )


@pytest.fixture
def versiondb_namespaces(dir_versiondb_in):
    {
        "upgrade": Namespace(
            dbpath=dir_versiondb_in / "pyanidb_upgrade",
            upgrade="head",
            downgrade=None,
            dry_run=None,
            direction=None,
        ),
        "downgrade": Namespace(
            dbpath=dir_versiondb_in / "pyanidb_downgrade",
            upgrade=None,
            downgrade="base",
            dry_run=None,
            direction=None,
        ),
        "dry_down": Namespace(
            dbpath=dir_versiondb_in / "pyanidb_dry_down",
            upgrade=None,
            downgrade=None,
            dry_run="head:base",
            direction="downgrade",
        ),
        "dry_up": Namespace(
            dbpath=dir_versiondb_in / "pyanidb_dry_up",
            upgrade=None,
            downgrade=None,
            dry_run="base:head",
            direction="upgrade",
        ),
        "altdb": Namespace(
            dbpath=dir_versiondb_in / "pyanidb_altdb",
            upgrade="head",
            downgrade=None,
            dry_run=None,
            direction=None,
            name="altdb",
        ),
        "alt_config": Namespace(
            dbpath=dir_versiondb_in / "pyanidb_alt_config",
            upgrade="head",
            downgrade=None,
            dry_run=None,
            direction=None,
            config="alt_config",
        ),
    }


# Test alembic command generation
def test_alembic_cmdline_generation():
    """Generate single alembic command line."""
    pass
    # alembic_cmd = versiondb.construct_alembic_cmdline()
    # dir_alembic = tmp_path / "versiondb_output"
    # expected = "alembic upgrade"


# Test upgrade
def test_versiondb_upgrade(dir_versiondb_in):
    """ """
    args = versiondb_namespaces["upgrade"]
    timestamp = "testing"
    shutil.copy(dir_versiondb_in / "base_pyanidb", dir_versiondb_in / "pyanidb_upgrade")
    versiondb.migrate_database(args.direction, args, timestamp)

    assert filecmp.cmp(
        dir_versiondb_in / "pyanidb_upgrade", dir_versiondb_in / "head_pyanidb"
    )


# Test downgrade
def test_versiondb_downgrade(dir_versiondb_in):
    """ """
    args = versiondb_namespaces["downgrade"]
    timestamp = "testing"
    shutil.copy(
        dir_versiondb_in / "head_pyanidb", dir_versiondb_in / "pyanidb_downgrade"
    )
    versiondb.migrate_database(args.direction, args, timestamp)

    assert filecmp.cmp(
        dir_versiondb_in / "pyanidb_downgrade", dir_versiondb_in / "base_pyanidb"
    )


# Test dry-run upgrade result
def test_versiondb_dry_upgrade(dir_versiondb_in):
    args = versiondb_namespaces["dry_up"]
    timestamp = "testing"
    shutil.copy(dir_versiondb_in / "base_pyanidb", dir_versiondb_in / "pyanidb_dry_up")
    versiondb.migrate_database(args.direction, args, timestamp)

    assert filecmp.cmp(
        dir_versiondb_in / "pyanidb_dry_up", dir_versiondb_in / "head_pyanidb"
    )


# Test dry-run upgrade result
def test_versiondb_dry_downgrade(dir_versiondb_in):
    args = versiondb_namespaces["dry_down"]
    timestamp = "testing"
    shutil.copy(
        dir_versiondb_in / "head_pyanidb", dir_versiondb_in / "pyanidb_dry_down"
    )
    versiondb.migrate_database(args.direction, args, timestamp)

    assert filecmp.cmp(
        dir_versiondb_in / "pyanidb_dry_down", dir_versiondb_in / "base_pyanidb"
    )


# Test dry-run upgrade result
def test_versiondb_altname(dir_versiondb_in):
    args = versiondb_namespaces["altname"]
    timestamp = "testing"
    shutil.copy(dir_versiondb_in / "head_pyanidb", dir_versiondb_in / "pyanidb_altdb")
    versiondb.migrate_database(args.direction, args, timestamp)

    assert filecmp.cmp(
        dir_versiondb_in / "pyanidb_altdb", dir_versiondb_in / "head_pyanidb"
    )


# Test dry-run upgrade result
def test_versiondb_alt_config(dir_versiondb_in):
    args = versiondb_namespaces["alt_config"]
    timestamp = "testing"
    shutil.copy(
        dir_versiondb_in / "head_pyanidb", dir_versiondb_in / "pyanidb_alt_config"
    )
    versiondb.migrate_database(args.direction, args, timestamp)

    assert filecmp.cmp(
        dir_versiondb_in / "pyanidb_alt_config", dir_versiondb_in / "head_pyanidb"
    )
