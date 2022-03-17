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
def alt_config_namespace(generic_versiondb_namespace):
    """Namespace for pyani versiondb -c alt_config."""
    return modify_namespace(
        generic_versiondb_namespace,
        dbpath="pyanidb_alt_config",
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
        "upgrade": b"2a3,7\n> CREATE TABLE alembic_version (\n> \tversion_num VARCHAR(32) NOT NULL, \n> \tCONSTRAINT alembic_version_pkc PRIMARY KEY (version_num)\n> );\n> INSERT INTO alembic_version VALUES('92f7f6b1626e');\n54,65c59\n< CREATE TABLE runs_comparisons (\n< \tcomparison_id INTEGER, \n< \trun_id INTEGER, \n< \tFOREIGN KEY(comparison_id) REFERENCES comparisons (comparison_id), \n< \tFOREIGN KEY(run_id) REFERENCES runs (run_id)\n< );\n< CREATE TABLE alembic_version (\n< \tversion_num VARCHAR(32) NOT NULL, \n< \tCONSTRAINT alembic_version_pkc PRIMARY KEY (version_num)\n< );\n< INSERT INTO alembic_version VALUES('92f7f6b1626e');\n< CREATE TABLE IF NOT EXISTS \"comparisons\" (\n---\n> CREATE TABLE comparisons (\n81d74\n< \tCHECK (maxmatch IN (0, 1)), \n85a79,84\n> CREATE TABLE runs_comparisons (\n> \tcomparison_id INTEGER, \n> \trun_id INTEGER, \n> \tFOREIGN KEY(comparison_id) REFERENCES comparisons (comparison_id), \n> \tFOREIGN KEY(run_id) REFERENCES runs (run_id)\n> );\n",
        "downgrade": b'3,6d2\n< CREATE TABLE alembic_version (\n< \tversion_num VARCHAR(32) NOT NULL, \n< \tCONSTRAINT alembic_version_pkc PRIMARY KEY (version_num)\n< );\n58,64c54\n< CREATE TABLE runs_comparisons (\n< \tcomparison_id INTEGER, \n< \trun_id INTEGER, \n< \tFOREIGN KEY(comparison_id) REFERENCES comparisons (comparison_id), \n< \tFOREIGN KEY(run_id) REFERENCES runs (run_id)\n< );\n< CREATE TABLE IF NOT EXISTS "comparisons" (\n---\n> CREATE TABLE comparisons (\n78,79c68\n< \tCHECK (maxmatch IN (0, 1)), \n< \tCONSTRAINT base_reqs UNIQUE (query_id, subject_id, program, version, fragsize, maxmatch), \n---\n> \tUNIQUE (query_id, subject_id, program, version, fragsize, maxmatch), \n82a72,77\n> CREATE TABLE runs_comparisons (\n> \tcomparison_id INTEGER, \n> \trun_id INTEGER, \n> \tFOREIGN KEY(comparison_id) REFERENCES comparisons (comparison_id), \n> \tFOREIGN KEY(run_id) REFERENCES runs (run_id)\n> );\n',
        "altdb": b"2a3,7\n> CREATE TABLE alembic_version (\n> \tversion_num VARCHAR(32) NOT NULL, \n> \tCONSTRAINT alembic_version_pkc PRIMARY KEY (version_num)\n> );\n> INSERT INTO alembic_version VALUES('92f7f6b1626e');\n66a72,73\n> \tkmersize INTEGER, \n> \tminmatch FLOAT, \n68c75\n< \tCONSTRAINT base_reqs UNIQUE (query_id, subject_id, program, version, fragsize, maxmatch), \n---\n> \tCONSTRAINT fastani_reqs UNIQUE (query_id, subject_id, program, version, fragsize, maxmatch, kmersize, minmatch), \n",
        "alt_config": b"2a3,7\n> CREATE TABLE alembic_version (\n> \tversion_num VARCHAR(32) NOT NULL, \n> \tCONSTRAINT alembic_version_pkc PRIMARY KEY (version_num)\n> );\n> INSERT INTO alembic_version VALUES('92f7f6b1626e');\n66a72,73\n> \tkmersize INTEGER, \n> \tminmatch FLOAT, \n68c75\n< \tCONSTRAINT base_reqs UNIQUE (query_id, subject_id, program, version, fragsize, maxmatch), \n---\n> \tCONSTRAINT fastani_reqs UNIQUE (query_id, subject_id, program, version, fragsize, maxmatch, kmersize, minmatch), \n",
    }.get(namespace, None)


# Create database dump
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
    # db = abs_dbpath.name
    # target = args.target.name
    # start = args.start.name

    dir_versiondb_out.mkdir(exist_ok=True)
    Path(f"{dir_versiondb_out}/{test}").mkdir(exist_ok=True)

    # The files must be moved, but shutil.move() does not work
    # if the file already exists, so this is a two-step process
    # Copy files to new location
    shutil.copy(abs_dbpath, dir_versiondb_out / test)
    shutil.copy(f"{abs_dbpath}.sql", dir_versiondb_out / test)
    shutil.copy(f"{args.target}.sql", dir_versiondb_out / test)
    shutil.copy(f"{args.start}.sql", dir_versiondb_out / test)

    # Remove old files
    os.remove(abs_dbpath)
    os.remove(f"{abs_dbpath}.sql")
    os.remove(f"{args.target}.sql")
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
def test_versiondb_upgrade(
    generic_versiondb_namespace,
    dir_versiondb_in,
    dir_versiondb_out,
    dir_versiondb_targets,
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

    # Dump altered and target databases
    enddb_dump = dumpdb(abs_dbpath)

    targetdb_dump = dumpdb(args.target)

    # Run diff
    diff_cmd = ["diff", "-y", "--suppress-common-lines", enddb_dump, targetdb_dump]
    result = subprocess.run(
        diff_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    actual_diff = "".join(
        open(dir_versiondb_targets / "upgrade_minus_head.diff", "r").readlines()
    )

    # Move files
    cleanup(abs_dbpath, "upgrade", dir_versiondb_out, args)

    assert result.stdout.decode() == actual_diff


def test_versiondb_downgrade(
    downgrade_namespace, dir_versiondb_in, dir_versiondb_out, dir_versiondb_targets
):
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

    # Dump altered and target databases
    enddb_dump = dumpdb(abs_dbpath)
    targetdb_dump = dumpdb(args.target)

    # Run diff
    diff_cmd = ["diff", "-y", "--suppress-common-lines", enddb_dump, targetdb_dump]
    result = subprocess.run(
        diff_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    actual_diff = "".join(
        open(dir_versiondb_targets / "downgrade_minus_base.diff", "r").readlines()
    )

    # Move output files
    cleanup(abs_dbpath, "downgrade", dir_versiondb_out, args)

    assert result.stdout.decode() == actual_diff


# Test alternate dbname
def test_versiondb_altdb(
    altdb_namespace, dir_versiondb_in, dir_versiondb_out, dir_versiondb_targets
):
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

    # Dump altered and target databases
    enddb_dump = dumpdb(abs_dbpath)
    targetdb_dump = dumpdb(args.target)

    # Run diff
    diff_cmd = ["diff", "-y", "--suppress-common-lines", enddb_dump, targetdb_dump]
    result = subprocess.run(
        diff_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    actual_diff = "".join(
        open(dir_versiondb_targets / "altdb_minus_head.diff", "r").readlines()
    )

    # Move files
    cleanup(abs_dbpath, "altdb", dir_versiondb_out, args)

    assert result.stdout.decode() == actual_diff


# Test alt_config result
@pytest.mark.skip(reason="May no be needed.")
def test_versiondb_alt_config(
    alt_config_namespace, dir_versiondb_in, dir_versiondb_out
):
    """Test upgrade of database using an alternate config file."""
    # Test setup
    # Retrieve test namespace and
    # Set environment variables and resolve absolute path of database
    args = alt_config_namespace
    setenv(dir_versiondb_in, args.dbpath)
    abs_dbpath = os.environ.get("PYANI_DATABASE")

    # Create dump file
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

    # Dump altered and target databases
    enddb_dump = dumpdb(abs_dbpath)
    targetdb_dump = dumpdb(args.target)

    # Run diff
    diff_cmd = ["diff", "--suppress-common-lines", enddb_dump, targetdb_dump]
    result = subprocess.run(
        diff_cmd,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
    )

    actual_diff = "".join(
        open(
            "/Users/baileythegreen/Software/pyani/tests/fixtures/versiondb/upgrade_minus_head.diff",
            "r",
        ).readlines()
    )

    # Move files
    cleanup(abs_dbpath, "alt_config", dir_versiondb_out, args)

    assert result.stdout.decode() == actual_diff
