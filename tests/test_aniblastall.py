#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) University of Strathclyde 2019-2021
# Author: Bailey Harrington
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# 161 Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2019-2021 University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Test aniblastall.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path
from typing import List, NamedTuple

import pandas as pd
import pytest  # noqa: F401  # pylint: disable=unused-import

from pandas.util.testing import assert_frame_equal

from pyani import anib, pyani_files  # probably don't need anib
from pyani import aniblastall


class ANIblastallOutput(NamedTuple):

    """Convenience struct for ANIblastall output."""

    fragfile: Path
    legacytabfile: Path


class ANIblastallOutputDir(NamedTuple):

    """Convenience struct for ANIblastall output."""

    infiles: List[Path]
    fragfiles: List[Path]
    legacyblastdir: Path
    legacyblastresult: pd.DataFrame


@pytest.fixture
def aniblastall_output(dir_aniblastall_in):
    """Namedtuple of example ANIblastall output.

    fragfile - fragmented FASTA query file
    legacytabfile - blastall tabular output
    """
    return ANIblastallOutput(
        dir_aniblastall_in / "NC_002696-fragments.fna",
        dir_aniblastall_in / "NC_002696_vs_NC_010338.blast_tab",
    )


@pytest.fixture
def aniblastall_output_dir(dir_aniblastall_in):
    """Namedtuple of example ANIblastall output - full directory.

    infiles - list of FASTA query files
    fragfiles - list of fragmented FASTA query files
    legacyblastdir - path to blastall output data
    legacyblastresult - pd.DataFrame result for blastall
    """
    return ANIblastallOutputDir(
        [
            _
            for _ in (dir_aniblastall_in / "sequences").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        [
            _
            for _ in (dir_aniblastall_in / "fragfiles").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        dir_aniblastall_in / "blastall",
        pd.read_csv(
            dir_aniblastall_in / "dataframes" / "blastall_result.csv", index_col=0
        ),
    )


# Test get_version()
# Test case 1: there is no executable
def test_get_version_missing_exe(executable_missing):
    """Test behaviour when there is no file at the specified executable location."""
    test_file_1 = Path("/non/existent/blastall")
    assert aniblastall.get_version(test_file_1) == f"No blastall at {test_file_1}"


# Test case 2: there is a file, but it is not executable
def test_get_version_not_executable(executable_not_executable):
    """Test behaviour when the file at the executable location is not executable."""
    test_file_2 = Path("/non/executable/blastall")
    assert (
        aniblastall.get_version(test_file_2)
        == f"blastall exists at {test_file_2} but not executable"
    )


# Test case 3: there is an executable file, but the version can't be retrieved
def test_get_version_no_version(executable_without_version):
    """Test behaviour when the version for the executable can not be retrieved."""
    test_file_3 = Path("/missing/version/blastall")
    assert (
        aniblastall.get_version(test_file_3)
        == f"blastall exists at {test_file_3} but could not retrieve version"
    )


# Test case 4: there is an executable file, but it will not run on the OS
def test_get_version_os_incompatible(executable_incompatible_with_os):
    """Test behaviour when the program can't run on the operating system.
    This will happen with newer versions of MacOS."""
    test_file_4 = Path("/os/incompatible/blastall")
    assert (
        aniblastall.get_version(test_file_4)
        == f"blastall exists at {test_file_4} but could not be executed"
    )


# Test legacy BLAST (blastall) command generation
def test_blastall_dbjobdict(path_fna_all, tmp_path):
    """Generate dictionary of legacy BLASTN database jobs."""
    blastcmds = aniblastall.make_blastcmd_builder(tmp_path)
    jobdict = aniblastall.build_db_jobs(path_fna_all, blastcmds)
    expected = [
        (tmp_path / _.name, f"formatdb -p F -i {tmp_path / _.name} -t {_.stem}")
        for _ in path_fna_all
    ]
    assert sorted([(k, v.script) for (k, v) in jobdict.items()]) == sorted(expected)


def test_blastall_graph(path_fna_all, tmp_path, fragment_length):
    """Create jobgraph for legacy BLASTN jobs."""
    fragresult = aniblastall.fragment_fasta_files(
        path_fna_all, tmp_path, fragment_length
    )
    blastcmds = aniblastall.make_blastcmd_builder(tmp_path)
    jobgraph = aniblastall.make_job_graph(path_fna_all, fragresult[0], blastcmds)
    # We check that the main script job is a blastn job, and that there
    # is a single dependency, which is a makeblastdb job
    for job in jobgraph:
        assert job.script.startswith("blastall -p blastn")
        assert len(job.dependencies) == 1
        assert job.dependencies[0].script.startswith("formatdb")


@pytest.mark.skip(reason="unsure this is needed")
def test_blastall_multiple(path_fna_two, tmp_path):
    """Generate legacy BLASTALL commands."""
    cmds = aniblastall.generate_blastall_commands(path_fna_two, tmp_path)
    expected = [
        (
            "blastall -p blastn -o "
            f"{tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
            f"-i {path_fna_two[0]} "
            f"-d {path_fna_two[1]} "
            "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
        ),
        (
            "blastall -p blastn -o "
            f"{tmp_path / str(path_fna_two[1].stem + '_vs_' + path_fna_two[0].stem + '.blast_tab')} "
            f"-i {path_fna_two[1]} "
            f"-d {path_fna_two[0]} "
            "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
        ),
    ]
    assert cmds == expected


def test_blastall_single(path_fna_two, tmp_path):
    """Generate legacy BLASTALL command-line."""
    cmd = aniblastall.construct_blastall_cmdline(
        path_fna_two[0], path_fna_two[1], tmp_path
    )
    expected = (
        f"blastall -p blastn -o {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
        f"-i {path_fna_two[0]} "
        f"-d {path_fna_two[1]} "
        "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
    )
    assert cmd == expected


# Test legacy BLAST database formatting (formatdb) command generation
@pytest.mark.skip(reason="unsure this is needed")
def test_formatdb_multiple(path_fna_two, tmp_path):
    """Generate legacy BLAST db creation commands."""
    cmds = aniblastall.generate_blastdb_commands(path_fna_two, tmp_path)
    expected = [
        (
            f"formatdb -p F -i {tmp_path / path_fna_two[0].name} -t {path_fna_two[0].stem}",
            tmp_path / path_fna_two[0].name,
        ),
        (
            f"formatdb -p F -i {tmp_path / path_fna_two[1].name} -t {path_fna_two[1].stem}",
            tmp_path / path_fna_two[1].name,
        ),
    ]
    assert cmds == expected


def test_formatdb_single(path_fna, tmp_path):
    """Generate legacy BLAST formatdb command-line."""
    cmd = aniblastall.construct_formatdb_cmd(path_fna, tmp_path)
    expected = f"formatdb -p F -i {tmp_path / path_fna.name} -t {path_fna.stem}"
    assert cmd[0] == expected


# Test output file parsing for ANIb methods
def test_parse_legacy_blastdir(aniblastall_output_dir):
    """Parses directory of legacy BLAST output."""
    orglengths = pyani_files.get_sequence_lengths(aniblastall_output_dir.infiles)
    fraglengths = aniblastall.get_fraglength_dict(aniblastall_output_dir.fragfiles)
    result = aniblastall.process_blast(
        aniblastall_output_dir.legacyblastdir, orglengths, fraglengths
    )
    assert_frame_equal(
        result.percentage_identity.sort_index(1).sort_index(),
        aniblastall_output_dir.legacyblastresult.sort_index(1).sort_index(),
    )


def test_parse_legacy_blasttab(aniblastall_output):
    """Parses ANIB legacy .blast_tab output."""
    fragdata = aniblastall.get_fraglength_dict([aniblastall_output.fragfile])
    result = aniblastall.parse_blast_tab(aniblastall_output.legacytabfile, fragdata)
    assert (
        a == b for a, b in zip(result, [1_966_922, 406_104, 78.578_978_313_253_018])
    )
