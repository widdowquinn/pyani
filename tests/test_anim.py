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

import os

from pathlib import Path

import pandas as pd
import pytest

from pandas.util.testing import assert_frame_equal

from pyani import anim, pyani_files, pyani_tools


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
    dir_nucmer = tmp_path / "nucmer_output"
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
    print(f"\n{cmds}")
    print((mummer_cmds_four.ncmds, mummer_cmds_four.fcmds))
    assert cmds == (mummer_cmds_four.ncmds, mummer_cmds_four.fcmds)


def test_mummer_single(tmp_path, path_file_two):
    """Generate single NUCmer/delta-filter command-line."""
    cmds = anim.construct_nucmer_cmdline(
        path_file_two[0], path_file_two[1], outdir=tmp_path
    )
    dir_nucmer = tmp_path / "nucmer_output"
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
