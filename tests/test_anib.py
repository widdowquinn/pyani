#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2022
# Author: Leighton Pritchard
#
# Contact: leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2019-2022 University of Strathclyde
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
"""Test anib.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path
from typing import List, NamedTuple

import pandas as pd
import pytest  # noqa: F401  # pylint: disable=unused-import
import unittest

from pandas.testing import assert_frame_equal

from pyani import anib, pyani_files


class ANIbOutput(NamedTuple):

    """Convenience struct for ANIb output."""

    fragfile: Path
    tabfile: Path
    legacytabfile: Path


class ANIbOutputDir(NamedTuple):

    """Convenience struct for ANIb output."""

    infiles: List[Path]
    fragfiles: List[Path]
    blastdir: Path
    legacyblastdir: Path
    blastresult: pd.DataFrame
    legacyblastresult: pd.DataFrame


# Create object for accessing unittest assertions
assertions = unittest.TestCase("__init__")


@pytest.fixture
def anib_output(dir_anib_in):
    """Namedtuple of example ANIb output.

    fragfile - fragmented FASTA query file
    tabfile  - BLAST+ tabular output
    legacytabfile - blastall tabular output
    """
    return ANIbOutput(
        dir_anib_in / "NC_002696-fragments.fna",
        dir_anib_in / "NC_002696_vs_NC_011916.blast_tab",
        dir_anib_in / "NC_002696_vs_NC_010338.blast_tab",
    )


@pytest.fixture
def anib_output_dir(dir_anib_in):
    """Namedtuple of example ANIb output - full directory.

    infiles - list of FASTA query files
    fragfiles - list of fragmented FASTA query files
    blastdir - path to BLAST+ output data
    legacyblastdir - path to blastall output data
    blastresult - pd.DataFrame result for BLAST+
    legacyblastresult - pd.DataFrame result for blastall
    """
    return ANIbOutputDir(
        [
            _
            for _ in (dir_anib_in / "sequences").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        [
            _
            for _ in (dir_anib_in / "fragfiles").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        dir_anib_in / "blastn",
        dir_anib_in / "blastall",
        pd.read_csv(dir_anib_in / "dataframes" / "blastn_result.csv", index_col=0),
        pd.read_csv(dir_anib_in / "dataframes" / "blastall_result.csv", index_col=0),
    )


# Test get_version()
# Test case 0: no executable location is specified
def test_get_version_nonetype():
    """Test behaviour when no location for the executable is given."""
    test_file_0 = None

    assert anib.get_version(test_file_0) == f"{test_file_0} is not found in $PATH"


# Test case 1: there is no executable
def test_get_version_no_exe(executable_missing, monkeypatch):
    """Test behaviour when there is no file at the specified executable location."""
    test_file_1 = Path("/non/existent/blastn")
    assert anib.get_version(test_file_1) == f"No blastn executable at {test_file_1}"


# Test case 2: there is a file, but it is not executable
def test_get_version_exe_not_executable(executable_not_executable, monkeypatch):
    """Test behaviour when the file at the executable location is not executable."""
    test_file_2 = Path("/non/executable/blastn")
    assert (
        anib.get_version(test_file_2)
        == f"blastn exists at {test_file_2} but not executable"
    )


# Test case 3: there is an executable file, but the version can't be retrieved
def test_get_version_exe_no_version(executable_without_version, monkeypatch):
    """Test behaviour when the version for the executable can not be retrieved."""
    test_file_3 = Path("/missing/version/blastn")
    assert (
        anib.get_version(test_file_3)
        == f"blastn exists at {test_file_3} but could not retrieve version"
    )


# Test legacy BLAST (blastall) command generation
def test_blastall_dbjobdict(path_fna_all, tmp_path):
    """Generate dictionary of legacy BLASTN database jobs."""
    blastcmds = anib.make_blastcmd_builder("ANIblastall", tmp_path)
    jobdict = anib.build_db_jobs(path_fna_all, blastcmds)
    expected = [
        (tmp_path / _.name, f"formatdb -p F -i {tmp_path / _.name} -t {_.stem}")
        for _ in path_fna_all
    ]
    assert sorted([(k, v.script) for (k, v) in jobdict.items()]) == sorted(expected)


def test_blastall_graph(path_fna_all, tmp_path, fragment_length):
    """Create jobgraph for legacy BLASTN jobs."""
    fragresult = anib.fragment_fasta_files(path_fna_all, tmp_path, fragment_length)
    blastcmds = anib.make_blastcmd_builder("ANIblastall", tmp_path)
    jobgraph = anib.make_job_graph(path_fna_all, fragresult[0], blastcmds)
    # We check that the main script job is a blastn job, and that there
    # is a single dependency, which is a makeblastdb job
    for job in jobgraph:
        assert job.script.startswith("blastall -p blastn")
        assert len(job.dependencies) == 1
        assert job.dependencies[0].script.startswith("formatdb")


def test_blastall_multiple(path_fna_two, tmp_path):
    """Generate legacy BLASTN commands."""
    cmds = anib.generate_blastn_commands(path_fna_two, tmp_path, mode="ANIblastall")
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
    """Generate legacy BLASTN command-line."""
    cmd = anib.construct_blastall_cmdline(path_fna_two[0], path_fna_two[1], tmp_path)
    expected = (
        f"blastall -p blastn -o {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
        f"-i {path_fna_two[0]} "
        f"-d {path_fna_two[1]} "
        "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
    )
    assert cmd == expected


# Test BLAST+ (blastn) command generation
def test_blastn_dbjobdict(path_fna_all, tmp_path):
    """Generate dictionary of BLASTN+ database jobs."""
    blastcmds = anib.make_blastcmd_builder("ANIb", tmp_path)
    jobdict = anib.build_db_jobs(path_fna_all, blastcmds)
    expected = [
        (
            tmp_path / _.name,
            f"makeblastdb -dbtype nucl -in {_} -title {_.stem} -out {tmp_path / _.name}",
        )
        for _ in path_fna_all
    ]
    assert sorted([(k, v.script) for (k, v) in jobdict.items()]) == sorted(expected)


def test_blastn_graph(path_fna_all, tmp_path, fragment_length):
    """Create jobgraph for BLASTN+ jobs."""
    fragresult = anib.fragment_fasta_files(path_fna_all, tmp_path, fragment_length)
    blastcmds = anib.make_blastcmd_builder("ANIb", tmp_path)
    jobgraph = anib.make_job_graph(path_fna_all, fragresult[0], blastcmds)
    # We check that the main script job is a blastn job, and that there
    # is a single dependency, which is a makeblastdb job
    for job in jobgraph:
        assert job.script.startswith("blastn")
        assert len(job.dependencies) == 1
        assert job.dependencies[0].script.startswith("makeblastdb")


def test_blastn_multiple(path_fna_two, tmp_path):
    """Generate BLASTN+ commands."""
    # BLAST+
    cmds = anib.generate_blastn_commands(path_fna_two, tmp_path, mode="ANIb")
    expected = [
        (
            f"blastn -out {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
            f"-query {path_fna_two[0]} "
            f"-db {path_fna_two[1]} "
            "-xdrop_gap_final 150 -dust no -evalue 1e-15 "
            "-max_target_seqs 1 -outfmt '6 qseqid sseqid "
            "length mismatch pident nident qlen slen qstart "
            "qend sstart send positive ppos gaps' -task blastn"
        ),
        (
            f"blastn -out {tmp_path / str(path_fna_two[1].stem + '_vs_' + path_fna_two[0].stem + '.blast_tab')} "
            f"-query {path_fna_two[1]} "
            f"-db {path_fna_two[0]} "
            "-xdrop_gap_final 150 -dust no -evalue 1e-15 "
            "-max_target_seqs 1 -outfmt '6 qseqid sseqid length "
            "mismatch pident nident qlen slen qstart qend "
            "sstart send positive ppos gaps' -task blastn"
        ),
    ]
    assert cmds == expected


def test_blastn_single(path_fna_two, tmp_path):
    """Generate BLASTN+ command-line."""
    cmd = anib.construct_blastn_cmdline(path_fna_two[0], path_fna_two[1], tmp_path)
    expected = (
        f"blastn -out {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
        f"-query {path_fna_two[0]} "
        f"-db {path_fna_two[1]} "
        "-xdrop_gap_final 150 -dust no -evalue 1e-15 -max_target_seqs 1 "
        "-outfmt '6 qseqid sseqid length mismatch pident nident qlen slen "
        "qstart qend sstart send positive ppos gaps' "
        "-task blastn"
    )
    assert cmd == expected


# Test legacy BLAST database formatting (formatdb) command generation
def test_formatdb_multiple(path_fna_two, tmp_path):
    """Generate legacy BLAST db creation commands."""
    cmds = anib.generate_blastdb_commands(path_fna_two, tmp_path, mode="ANIblastall")
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
    cmd = anib.construct_formatdb_cmd(path_fna, tmp_path)
    expected = f"formatdb -p F -i {tmp_path / path_fna.name} -t {path_fna.stem}"
    assert cmd[0] == expected


# Test FASTA file fragmentation for ANIb methods
def test_fragment_files(path_fna_all, tmp_path, dir_tgt_fragments, fragment_length):
    """Fragment files for ANIb/ANIblastall."""
    result = anib.fragment_fasta_files(path_fna_all, tmp_path, fragment_length)

    # # Test fragment lengths are in bounds
    for _, fragdict in result[-1].items():
        for _, fraglen in fragdict.items():
            assert fraglen <= fragment_length


# Test BLAST+ database formatting (makeblastdb) command generation
def test_makeblastdb_multiple(path_fna_two, tmp_path):
    """Generate multiple BLAST+ makeblastdb command-lines."""
    cmds = anib.generate_blastdb_commands(path_fna_two, tmp_path, mode="ANIb")
    expected = [
        (
            (
                f"makeblastdb -dbtype nucl -in {path_fna_two[0]} "
                f"-title {path_fna_two[0].stem} -out {tmp_path / path_fna_two[0].name}"
            ),
            tmp_path / path_fna_two[0].name,
        ),
        (
            (
                f"makeblastdb -dbtype nucl -in {path_fna_two[1]} "
                f"-title {path_fna_two[1].stem} -out {tmp_path / path_fna_two[1].name}"
            ),
            tmp_path / path_fna_two[1].name,
        ),
    ]
    assert cmds == expected


def test_makeblastdb_single(path_fna, tmp_path):
    """Generate single BLAST+ makeblastdb command-line."""
    cmd = anib.construct_makeblastdb_cmd(path_fna, tmp_path)
    expected = (
        f"makeblastdb -dbtype nucl -in {path_fna} "
        f"-title {path_fna.stem} -out {tmp_path / path_fna.name}"
    )
    assert cmd[0] == expected


# Test output file parsing for ANIb methods
def test_parse_legacy_blastdir(anib_output_dir):
    """Parses directory of legacy BLAST output."""
    orglengths = pyani_files.get_sequence_lengths(anib_output_dir.infiles)
    fraglengths = anib.get_fraglength_dict(anib_output_dir.fragfiles)
    result = anib.process_blast(
        anib_output_dir.legacyblastdir, orglengths, fraglengths, mode="ANIblastall"
    )
    assert_frame_equal(
        result.percentage_identity.sort_index(axis=1).sort_index(),
        anib_output_dir.legacyblastresult.sort_index(axis=1).sort_index(),
    )


def test_parse_blastdir(anib_output_dir):
    """Parse directory of BLAST+ output."""
    orglengths = pyani_files.get_sequence_lengths(anib_output_dir.infiles)
    fraglengths = anib.get_fraglength_dict(anib_output_dir.fragfiles)
    result = anib.process_blast(
        anib_output_dir.blastdir, orglengths, fraglengths, mode="ANIb"
    )
    assert_frame_equal(
        result.percentage_identity.sort_index(axis=1).sort_index(),
        anib_output_dir.blastresult.sort_index(axis=1).sort_index(),
    )


def test_parse_blasttab(anib_output):
    """Parse ANIb BLAST+ .blast_tab output."""
    fragdata = anib.get_fraglength_dict([anib_output.fragfile])
    result = anib.parse_blast_tab(anib_output.tabfile, fragdata, mode="ANIb")
    assert (a == b for a, b in zip(result, [4_016_551, 93, 99.997_693_577_050_029]))


def test_parse_legacy_blasttab(anib_output):
    """Parses ANIB legacy .blast_tab output."""
    fragdata = anib.get_fraglength_dict([anib_output.fragfile])
    result = anib.parse_blast_tab(
        anib_output.legacytabfile, fragdata, mode="ANIblastall"
    )
    assert (
        a == b for a, b in zip(result, [1_966_922, 406_104, 78.578_978_313_253_018])
    )
