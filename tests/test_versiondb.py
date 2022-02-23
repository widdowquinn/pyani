"""Test versiondb.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import os

from pathlib import Path
from typing import List, NamedTuple, Tuple

import pandas as pd
import pytest
import unittest

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


# Test alembic command generation
def test_alembic_cmdline_generation():
    """Generate single alembic command line."""
    pass
    # alembic_cmd = versiondb.construct_alembic_cmdline()
    # dir_alembic = tmp_path / "versiondb_output"
    # expected = "alembic upgrade"


# Test upgrade
def test_versiondb_upgrade():
    """ """
    pass


# Test downgrade


# Test dry-run result
