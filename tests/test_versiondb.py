"""Test versiondb.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import os
import sys
import subprocess
import platform

from argparse import Namespace
from pathlib import Path
from typing import List, NamedTuple, Tuple

import pandas as pd
import pytest
import unittest
import shutil

from pandas.util.testing import assert_frame_equal

from pyani import versiondb, pyani_files, pyani_tools

from pyani.pyani_orm import PyaniORMException, get_session, add_alembic
from tools import modify_namespace

from pyani import pyani_config


# Create environment variables for alembic to access
def setenv(dir_versiondb_in, dbfile: Path):
    abs_path = Path(dir_versiondb_in / dbfile).resolve()
    os.environ["PYANI_DATABASE"] = str(abs_path)


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
def generic_versiondb_namespace(dir_versiondb_in):
    """Generic namespace for the pyani versiondb subcommand."""
    return Namespace(
        dbpath="pyanidb_upgrade",
        upgrade="head",
        downgrade=None,
        dry_run=None,
        direction="upgrade",
        dbname=None,
        alembic_config=None,
        start=dir_versiondb_in / "base_pyanidb",
        target=dir_versiondb_in / "head_pyanidb",
    )


@pytest.fixture
def downgrade_namespace(generic_versiondb_namespace, dir_versiondb_in):
    """Namespace for pyani versiondb downgrade."""
    return modify_namespace(
        generic_versiondb_namespace,
        dbpath="pyanidb_downgrade",
        upgrade=None,
        downgrade="base",
        direction="downgrade",
        start=dir_versiondb_in / "head_pyanidb",
        target=dir_versiondb_in / "base_pyanidb",
    )


@pytest.fixture
def altdb_namespace(generic_versiondb_namespace):
    """Namespace for pyani versiondb -n altdb."""
    return modify_namespace(
        generic_versiondb_namespace,
        dbpath="pyanidb_altdb",
        dbname="pyanidb_altdb",
        alembic_config="alt_alembic_config.ini",
    )


@pytest.fixture
def dry_up_namespace(generic_versiondb_namespace):
    """Namespace for pyani versiondb dry-run upgrade."""
    return modify_namespace(
        generic_versiondb_namespace,
        dbpath="pyanidb_dry_up",
        upgrade=None,
        dry_run="base:head",
    )


@pytest.fixture
def dry_down_namespace(generic_versiondb_namespace):
    """Namespace for pyani versiondb dry-run downgrade."""
    return modify_namespace(
        generic_versiondb_namespace,
        dbpath="pyanidb_dry_down",
        direction="downgrade",
        downgrade=None,
        dry_run="head:base",
    )


def expected_diffs(namespace):
    """Expected (acceptable) differences between output and target databases."""
    return {
        "upgrade": "",
        "downgrade": "DROP TABLE alembic_version;\n",
        "altdb": "",
    }.get(namespace, None)


# Create database dump—a version that can be edited using sed
def dumpdb(abs_dbpath):
    """Dump contents of database to a plain-text file."""

    cmdline = [pyani_config.SQLITE_DEFAULT, f"{abs_dbpath}", ".dump"]
    with open(f"{abs_dbpath}.sql", "w") as outfile:
        subprocess.run(
            cmdline,
            shell=False,
            stdout=outfile,
            stderr=subprocess.PIPE,
        )
    return f"{abs_dbpath}.sql"


def name_base_reqs(startdb_dump):
    """Name unique constraint in comparisons table of old database schema."""
    old_constraint = (
        "UNIQUE (query_id, subject_id, program, version, fragsize, maxmatch),"
    )
    new_constraint = "CONSTRAINT base_reqs UNIQUE (query_id, subject_id, program, version, fragsize, maxmatch),"

    # Edit .dump file so that the unique constraint is named
    # This is required in order to subsequently modify it
    # In-place usage differs on macOs vs Linux
    if platform.system() == "Darwin":
        sed_cmd = [
            "sed",
            "-i",
            ".bak",
            f"s/{old_constraint}/{new_constraint}/",
            startdb_dump,
        ]
    else:
        sed_cmd = [
            "sed",
            "-i",
            f"s/{old_constraint}/{new_constraint}/",
            startdb_dump,
        ]
    subprocess.run(
        sed_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )


def cleanup(abs_dbpath, test, dir_versiondb_out, args):
    """Remove files created for test."""

    dir_versiondb_out.mkdir(exist_ok=True)
    Path(f"{dir_versiondb_out}/{test}").mkdir(exist_ok=True)

    # The files must be moved, but shutil.move() does not work
    # if the file already exists, so this is a two-step process
    # Copy files to new location
    shutil.copy(abs_dbpath, dir_versiondb_out / test)
    shutil.copy(f"{args.start}.sql", dir_versiondb_out / test)

    # Remove old files
    os.remove(abs_dbpath)
    os.remove(f"{args.start}.sql")

    # This file is not generated in the downgrade test
    try:
        shutil.copy(f"{args.start}.sql.bak", dir_versiondb_out / test)
        os.remove(f"{args.start}.sql.bak")
    except FileNotFoundError:
        pass


# Test alembic command generation
def test_alembic_cmdline_generation(
    generic_versiondb_namespace,
    downgrade_namespace,
    altdb_namespace,
    dir_versiondb_in,
):
    """Generate alembic command lines."""

    alembic_cmds = []
    upgrade_args = generic_versiondb_namespace
    alembic_cmds.append(
        " ".join(
            versiondb.construct_alembic_cmdline(upgrade_args.direction, upgrade_args)
        )
    )

    downgrade_args = downgrade_namespace
    alembic_cmds.append(
        " ".join(
            versiondb.construct_alembic_cmdline(
                downgrade_args.direction, downgrade_args
            )
        )
    )

    altdb_args = altdb_namespace
    alembic_cmds.append(
        " ".join(versiondb.construct_alembic_cmdline(altdb_args.direction, altdb_args))
    )

    assert alembic_cmds == [
        "alembic upgrade head",
        "alembic downgrade base",
        f"alembic -n {altdb_args.dbname} -c {altdb_args.alembic_config} upgrade head",
    ]


# Test upgrade
@pytest.mark.skip_if_exe_missing("sqldiff")
def test_versiondb_upgrade(
    generic_versiondb_namespace,
    dir_versiondb_in,
    dir_versiondb_out,
):
    """Test upgrade of database."""
    # Test setup
    # Retrieve test namespace and
    # Set environment variables and resolve absolute path of database
    args = generic_versiondb_namespace
    setenv(dir_versiondb_in, args.dbpath)
    abs_dbpath = os.environ.get("PYANI_DATABASE")

    # Create and edit dump file to fix constraint problem
    startdb_dump = dumpdb(args.start)
    name_base_reqs(startdb_dump)

    # Run `sqlite3 -init <file>
    init_cmd = [pyani_config.SQLITE_DEFAULT, abs_dbpath]
    subprocess.run(
        init_cmd,
        stdin=open(startdb_dump),
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    # Run test migration
    versiondb.migrate_database(args.direction, args, timestamp="testing")

    # Run diff
    diff_cmd = [pyani_config.SQLDIFF_DEFAULT, "--schema", abs_dbpath, args.target]
    result = subprocess.run(
        diff_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    expected_diff = ""

    sys.stdout.write(f"Expected_diff: {expected_diff}\n\n")
    sys.stdout.write(f"Actual diff: {result.stdout.decode()}\n\n")

    # Move files
    cleanup(abs_dbpath, "upgrade", dir_versiondb_out, args)

    assert result.stdout.decode() == expected_diff


@pytest.mark.skip_if_exe_missing("sqldiff")
def test_versiondb_downgrade(downgrade_namespace, dir_versiondb_in, dir_versiondb_out):
    """Test downgrade of database."""
    # Test setup
    # Retrieve test namespace and
    # Set environment variables and resolve absolute path of database
    args = downgrade_namespace
    setenv(dir_versiondb_in, args.dbpath)
    abs_dbpath = os.environ.get("PYANI_DATABASE")

    # Create dump file
    startdb_dump = dumpdb(args.start)

    # Run `sqlite3 -init <file>
    init_cmd = [pyani_config.SQLITE_DEFAULT, abs_dbpath]
    subprocess.run(
        init_cmd,
        stdin=open(startdb_dump),
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    # Run test migration
    versiondb.migrate_database(args.direction, args, timestamp="testing")

    # Run diff
    diff_cmd = [pyani_config.SQLDIFF_DEFAULT, "--schema", abs_dbpath, args.target]
    result = subprocess.run(
        diff_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    expected_diff = "DROP TABLE alembic_version;\n"

    sys.stdout.write(f"Expected_diff: {expected_diff}\n\n")
    sys.stdout.write(f"Actual diff: {result.stdout.decode()}\n\n")

    # Move output files
    cleanup(abs_dbpath, "downgrade", dir_versiondb_out, args)

    assert result.stdout.decode() == expected_diff


# Test alternate dbname
@pytest.mark.skip_if_exe_missing("sqldiff")
def test_versiondb_altdb(altdb_namespace, dir_versiondb_in, dir_versiondb_out):
    """Test upgrade of database using an alternate database name and config file, such as in a multidb situation."""
    # Test setup
    # Retrieve test namespace and
    # Set environment variables and resolve absolute path of database
    args = altdb_namespace
    setenv(dir_versiondb_in, args.dbpath)
    abs_dbpath = os.environ.get("PYANI_DATABASE")

    # Create dump file
    startdb_dump = dumpdb(args.start)
    name_base_reqs(startdb_dump)

    # assert False
    # Run `sqlite3 -init <file>
    init_cmd = [pyani_config.SQLITE_DEFAULT, abs_dbpath]
    subprocess.run(
        init_cmd,
        stdin=open(startdb_dump),
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    # Run test migration
    versiondb.migrate_database(args.direction, args, timestamp="testing")

    # Run diff
    diff_cmd = [pyani_config.SQLDIFF_DEFAULT, "--schema", abs_dbpath, args.target]
    result = subprocess.run(
        diff_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    expected_diff = ""

    sys.stdout.write(f"Expected_diff: {expected_diff}\n\n")
    sys.stdout.write(f"Actual diff: {result.stdout.decode()}\n\n")

    # Move files
    cleanup(abs_dbpath, "altdb", dir_versiondb_out, args)

    assert result.stdout.decode() == expected_diff


# Dry-run tests still to be done