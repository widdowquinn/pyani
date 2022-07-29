"""Test fastani.py module.

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

from pyani import fastani, pyani_files, pyani_tools


### Some classes... to be decided


class ComparisonResult(NamedTuple):
    reference: Path
    query: Path
    ani: float
    matches: int
    fragments: int


class fastANIParsed(NamedTuple):

    """Convenience struct for fastANI .fastani file and associated parsed output."""

    filename: Path
    data: ComparisonResult


class FastANIExample(NamedTuple):

    """Convenience struct for fastANI command-line examples."""

    infiles: List[Path]
    fastcmds: List[str]


@pytest.fixture
def fastanifile_parsed(dir_fastani_in):  # works
    """Example parsed fastANI file."""
    return fastANIParsed(
        dir_fastani_in / "ecoli_vs_shiga.fastani",
        ComparisonResult("ecoli.fna", "shiga.fna", 0.9766400000000001, 1322, 1547),
    )


@pytest.fixture
def fastani_cmds_four(path_file_four):  # works
    """Example fastANI commands (four files)."""
    return FastANIExample(
        path_file_four,
        [
            "fastANI -q file1.fna -r file1.fna -o fastani_output/file1_vs_file1.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file1.fna -r file2.fna -o fastani_output/file1_vs_file2.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file1.fna -r file3.fna -o fastani_output/file1_vs_file3.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file1.fna -r file4.fna -o fastani_output/file1_vs_file4.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file2.fna -r file1.fna -o fastani_output/file2_vs_file1.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file2.fna -r file2.fna -o fastani_output/file2_vs_file2.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file2.fna -r file3.fna -o fastani_output/file2_vs_file3.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file2.fna -r file4.fna -o fastani_output/file2_vs_file4.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file3.fna -r file1.fna -o fastani_output/file3_vs_file1.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file3.fna -r file2.fna -o fastani_output/file3_vs_file2.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file3.fna -r file3.fna -o fastani_output/file3_vs_file3.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file3.fna -r file4.fna -o fastani_output/file3_vs_file4.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file4.fna -r file1.fna -o fastani_output/file4_vs_file1.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file4.fna -r file2.fna -o fastani_output/file4_vs_file2.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file4.fna -r file3.fna -o fastani_output/file4_vs_file3.fastani --fragLen 3000 -k 16 --minFraction 0.2",
            "fastANI -q file4.fna -r file4.fna -o fastani_output/file4_vs_file4.fastani --fragLen 3000 -k 16 --minFraction 0.2",
        ],
    )


# Create object for accessing unittest assertions
assertions = unittest.TestCase("__init__")


# Test get_version()
# Test case 0: no executable location is specified
def test_get_version_nonetype():
    """Test behaviour when no location for the executable is given."""
    test_file_0 = None

    assert (
        fastani.get_version(test_file_0)
        == f"expected path to fastANI executable; received {test_file_0}"
    )


# Test case 1: no such file exists
def test_get_version_random_string():
    """Test behaviour when the given 'file' is not one."""
    test_file_1 = "string"

    assert fastani.get_version(test_file_1) == f"{test_file_1} is not found in $PATH"


# Test case 2: there is no executable
def test_get_version_no_exe(executable_missing, monkeypatch):
    """Test behaviour when there is no file at the specified executable location."""
    test_file_2 = Path("/non/existent/fastani")
    assert fastani.get_version(test_file_2) == f"No fastANI executable at {test_file_2}"


# Test case 3: there is a file, but it is not executable
def test_get_version_exe_not_executable(executable_not_executable, monkeypatch):
    """Test behaviour when the file at the executable location is not executable."""
    test_file_3 = Path("/non/executable/fastani")
    assert (
        fastani.get_version(test_file_3)
        == f"fastANI exists at {test_file_3} but not executable"
    )


# Test case 4: there is an executable file, but the version can't be retrieved
def test_get_version_exe_no_version(executable_without_version, monkeypatch):
    """Test behaviour when the version for the executable can not be retrieved."""
    test_file_4 = Path("/missing/version/fastani")
    assert (
        fastani.get_version(test_file_4)
        == f"fastANI exists at {test_file_4} but could not retrieve version"
    )


def test_fastanifile_parsing(fastanifile_parsed):  # works
    """Check parsing of test fastANI .fastani file."""
    result = fastani.parse_fastani_file(fastanifile_parsed.filename)
    assertions.assertEqual(result, fastanifile_parsed.data)


# Test fastANI command generation
def test_fastani_single(tmp_path, path_file_two):  # works
    """Generate single fastANI command line."""
    fastcmd = fastani.construct_fastani_cmdline(
        path_file_two[0], path_file_two[1], outdir=tmp_path
    )
    dir_fastani = tmp_path / "fastani_output"
    expected = (
        f"fastANI -q {path_file_two[0]} -r {path_file_two[1]} "
        f"-o {dir_fastani / str(path_file_two[0].stem + '_vs_' + path_file_two[1].stem + '.fastani')} "
        f"--fragLen 3000 -k 16 --minFraction 0.2"
    )
    assert fastcmd == expected


def test_fastani_multiple(fastani_cmds_four):
    """Generate multiple abstract fastANI command-liens.

    Tests that all the input files are correctly-paired.
    """
    cmds = fastani.generate_fastani_commands(fastani_cmds_four.infiles)
    print(f"\n{cmds}")
    print((fastani_cmds_four.fastcmds))
    assert cmds == (fastani_cmds_four.fastcmds)


def test_fastani_job_generation(fastani_cmds_four):  # works
    """Generate job names"""

    joblist = fastani.generate_fastani_jobs(fastani_cmds_four.infiles, jobprefix="test")
    assert len(joblist) == 16

    for idx, job in enumerate(joblist):
        assert job.name == f"test_{idx:06d}"
