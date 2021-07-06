#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
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
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2019-2020 University of Strathclyde
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
"""Test anim.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path
from typing import List, NamedTuple, Tuple

import pandas as pd
import pytest

from pandas.testing import assert_frame_equal

from pyani import anim, pyani_files


class DeltaDir(NamedTuple):

    """Convenience struct for MUMmer .delta file and associated parsed output."""

    seqdir: Path
    deltadir: Path
    deltaresult: pd.DataFrame


class DeltaParsed(NamedTuple):

    """Convenience struct for MUMmer .delta file and associated parsed output."""

    filename: Path
    data: Tuple[int]


class MUMmerExample(NamedTuple):

    """Convenience struct for MUMmer command-line examples."""

    infiles: List[Path]
    ncmds: List[str]
    fcmds: List[str]


@pytest.fixture
def delta_output_dir(dir_anim_in):
    """Namedtuple of example MUMmer .delta file output."""
    return DeltaDir(
        dir_anim_in / "sequences",
        dir_anim_in / "deltadir",
        pd.read_csv(dir_anim_in / "dataframes" / "deltadir_result.csv", index_col=0),
    )


@pytest.fixture
def deltafile_parsed(dir_anim_in):
    """Example parsed deltafile data."""
    return DeltaParsed(dir_anim_in / "test.delta", (4074148, 2191))


@pytest.fixture
def mummer_cmds_four(path_file_four):
    """Example MUMmer commands (four files)."""
    return MUMmerExample(
        path_file_four,
        [
            "nucmer --mum -p nucmer_output/file1/file1_vs_file2 file1.fna file2.fna",
            "nucmer --mum -p nucmer_output/file1/file1_vs_file3 file1.fna file3.fna",
            "nucmer --mum -p nucmer_output/file1/file1_vs_file4 file1.fna file4.fna",
            "nucmer --mum -p nucmer_output/file2/file2_vs_file3 file2.fna file3.fna",
            "nucmer --mum -p nucmer_output/file2/file2_vs_file4 file2.fna file4.fna",
            "nucmer --mum -p nucmer_output/file3/file3_vs_file4 file3.fna file4.fna",
        ],
        [
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file1/file1_vs_file2.delta "
                "nucmer_output/file1/file1_vs_file2.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file1/file1_vs_file3.delta "
                "nucmer_output/file1/file1_vs_file3.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file1/file1_vs_file4.delta "
                "nucmer_output/file1/file1_vs_file4.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file2/file2_vs_file3.delta "
                "nucmer_output/file2/file2_vs_file3.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file2/file2_vs_file4.delta "
                "nucmer_output/file2/file2_vs_file4.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file3/file3_vs_file4.delta "
                "nucmer_output/file3/file3_vs_file4.filter"
            ),
        ],
    )


# Test .delta output file processing
def test_deltadir_parsing(delta_output_dir):
    """Process test directory of .delta files into ANIResults."""
    seqfiles = pyani_files.get_fasta_files(delta_output_dir.seqdir)
    orglengths = pyani_files.get_sequence_lengths(seqfiles)
    result = anim.process_deltadir(delta_output_dir.deltadir, orglengths)
    assert_frame_equal(
        result.percentage_identity.sort_index(1).sort_index(),
        delta_output_dir.deltaresult.sort_index(1).sort_index(),
    )


def test_deltafile_parsing(deltafile_parsed):
    """Check parsing of test NUCmer .delta/.filter file."""
    result = anim.parse_delta(deltafile_parsed.filename)
    assert result == deltafile_parsed.data


# Test MUMmer command generation
def test_maxmatch_single(tmp_path, path_file_two):
    """Generate NUCmer command line with maxmatch."""
    ncmd, _ = anim.construct_nucmer_cmdline(
        path_file_two[0], path_file_two[1], outdir=tmp_path, maxmatch=True
    )
    dir_nucmer = tmp_path / "nucmer_output" / path_file_two[0].stem
    expected = (
        "nucmer --maxmatch -p "
        f"{dir_nucmer / str(path_file_two[0].stem + '_vs_' + path_file_two[1].stem)} "
        f"{path_file_two[0]} {path_file_two[1]}"
    )
    assert ncmd == expected


def test_mummer_multiple(mummer_cmds_four):
    """Generate multiple abstract NUCmer/delta-filter command-lines.

    Tests that all the input files are correctly-paired
    """
    cmds = anim.generate_nucmer_commands(mummer_cmds_four.infiles)
    # print(f"\n{cmds}")
    # print((mummer_cmds_four.ncmds, mummer_cmds_four.fcmds))
    assert cmds == (mummer_cmds_four.ncmds, mummer_cmds_four.fcmds)


def test_mummer_single(tmp_path, path_file_two):
    """Generate single NUCmer/delta-filter command-line."""
    cmds = anim.construct_nucmer_cmdline(
        path_file_two[0], path_file_two[1], outdir=tmp_path
    )
    dir_nucmer = tmp_path / "nucmer_output" / path_file_two[0].stem
    expected = (
        (
            "nucmer --mum -p "
            f"{dir_nucmer / str(path_file_two[0].stem + '_vs_' + path_file_two[1].stem)} "
            f"{path_file_two[0]} {path_file_two[1]}"
        ),
        (
            "delta_filter_wrapper.py delta-filter -1 "
            f"{dir_nucmer / str(path_file_two[0].stem + '_vs_' + path_file_two[1].stem + '.delta')} "
            f"{dir_nucmer / str(path_file_two[0].stem + '_vs_' + path_file_two[1].stem + '.filter')}"
        ),
    )
    assert cmds == expected


def test_mummer_job_generation(mummer_cmds_four):
    """Generate dependency tree of NUCmer/delta-filter jobs.

    Tests that the correct dependency graph and naming scheme is produced.
    """
    joblist = anim.generate_nucmer_jobs(mummer_cmds_four.infiles, jobprefix="test")
    assert len(joblist) == 6

    for idx, job in enumerate(joblist):
        assert job.name == "test_%06d-f" % idx  # filter job name
        assert len(job.dependencies) == 1  # has NUCmer job
        assert job.dependencies[0].name == "test_%06d-n" % idx
