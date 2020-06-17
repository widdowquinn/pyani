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
import unittest

from pathlib import Path

import pandas as pd
import pytest

from pandas.util.testing import assert_frame_equal

from pyani import anim, pyani_files, pyani_tools


# Test MUMmer command generation
def test_mummer_single(tmp_path):
    """Generate single NUCmer/delta-filter command-line."""
    cmds = anim.construct_nucmer_cmdline(
        Path("file1.fna"), Path("file2.fna"), outdir=tmp_path
    )
    expected = (
        f"nucmer --mum -p {tmp_path / 'nucmer_output' / 'file1_vs_file2'} file1.fna file2.fna",
        (
            f"delta_filter_wrapper.py delta-filter -1 {tmp_path / 'nucmer_output' / 'file1_vs_file2.delta'} "
            f"{tmp_path / 'nucmer_output' / 'file1_vs_file2.filter'}"
        ),
    )
    assert cmds == expected


@pytest.mark.skipif(
    not pyani_tools.has_dependencies().mummer, reason="nucmer executable not available"
)
class TestNUCmerCmdline(unittest.TestCase):

    """Class defining tests of NUCmer command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        # Basic NUCmer and delta-filter command-line targets
        testdir = Path("tests")
        self.ntgt = " ".join(
            [
                "nucmer --mum -p",
                "tests/test_output/anim/nucmer_output/file1_vs_file2",
                "file1.fna file2.fna",
            ]
        )
        self.ntgtmax = " ".join(
            [
                "nucmer --maxmatch -p",
                "tests/test_output/anim/nucmer_output/file1_vs_file2",
                "file1.fna file2.fna",
            ]
        )
        self.ftgt = " ".join(
            [
                "delta_filter_wrapper.py delta-filter -1",
                "tests/test_output/anim/nucmer_output/file1_vs_file2.delta",
                "tests/test_output/anim/nucmer_output/file1_vs_file2.filter",
            ]
        )
        self.files = [Path("file1"), Path("file2"), Path("file3"), Path("file4")]
        self.ncmdlist = [
            "nucmer --mum -p nucmer_output/file1_vs_file2 file1 file2",
            "nucmer --mum -p nucmer_output/file1_vs_file3 file1 file3",
            "nucmer --mum -p nucmer_output/file1_vs_file4 file1 file4",
            "nucmer --mum -p nucmer_output/file2_vs_file3 file2 file3",
            "nucmer --mum -p nucmer_output/file2_vs_file4 file2 file4",
            "nucmer --mum -p nucmer_output/file3_vs_file4 file3 file4",
        ]
        self.fcmdlist = [
            " ".join(
                [
                    "delta_filter_wrapper.py delta-filter -1",
                    "nucmer_output/file1_vs_file2.delta",
                    "nucmer_output/file1_vs_file2.filter",
                ]
            ),
            " ".join(
                [
                    "delta_filter_wrapper.py delta-filter -1",
                    "nucmer_output/file1_vs_file3.delta",
                    "nucmer_output/file1_vs_file3.filter",
                ]
            ),
            " ".join(
                [
                    "delta_filter_wrapper.py delta-filter -1",
                    "nucmer_output/file1_vs_file4.delta",
                    "nucmer_output/file1_vs_file4.filter",
                ]
            ),
            " ".join(
                [
                    "delta_filter_wrapper.py delta-filter -1",
                    "nucmer_output/file2_vs_file3.delta",
                    "nucmer_output/file2_vs_file3.filter",
                ]
            ),
            " ".join(
                [
                    "delta_filter_wrapper.py delta-filter -1",
                    "nucmer_output/file2_vs_file4.delta",
                    "nucmer_output/file2_vs_file4.filter",
                ]
            ),
            " ".join(
                [
                    "delta_filter_wrapper.py delta-filter -1",
                    "nucmer_output/file3_vs_file4.delta",
                    "nucmer_output/file3_vs_file4.filter",
                ]
            ),
        ]
        self.outdir = testdir / "test_output" / "anim"
        self.indir = testdir / "test_input" / "anim"

    def test_maxmatch_cmd_generation(self):
        """Generate NUCmer command line with maxmatch."""
        ncmd, _ = anim.construct_nucmer_cmdline(
            Path("file1.fna"), Path("file2.fna"), outdir=self.outdir, maxmatch=True
        )
        self.assertEqual(ncmd, self.ntgtmax)

    def test_multi_cmd_generation(self):
        """Generate multiple abstract NUCmer/delta-filter command-lines.

        Tests that all the input files are correctly-paired
        """
        cmds = anim.generate_nucmer_commands(self.files)
        self.assertEqual(cmds, (self.ncmdlist, self.fcmdlist))

    def test_nucmer_job_generation(self):
        """Generate dependency tree of NUCmer/delta-filter jobs.

        Tests that the correct dependency graph and naming scheme is produced.
        """
        joblist = anim.generate_nucmer_jobs(self.files, jobprefix="test")
        self.assertEqual(len(joblist), 6)
        for idx, job in enumerate(joblist):
            self.assertEqual(job.name, "test_%06d-f" % idx)  # filter job name
            self.assertEqual(len(job.dependencies), 1)  # has NUCmer job
            self.assertEqual(
                job.dependencies[0].name, "test_%06d-n" % idx
            )  # NUCmer job name


class TestDeltafileProcessing(unittest.TestCase):

    """Class defining tests for .delta/.filter file parsing."""

    def setUp(self):
        """Set parameters for tests."""
        testdir = Path("tests")
        self.indir = testdir / "test_input" / "anim"
        self.seqdir = testdir / "test_input" / "sequences"
        self.deltafile = self.indir / "test.delta"
        self.deltadir = self.indir / "deltadir"
        self.df_pid = pd.DataFrame(
            [
                [1.000000, 0.851120, 0.999974, 0.868010],
                [0.851120, 1.000000, 0.851200, 0.852965],
                [0.999974, 0.851200, 1.000000, 0.868061],
                [0.868010, 0.852965, 0.868061, 1.000000],
            ],
            columns=["NC_002696", "NC_010338", "NC_011916", "NC_014100"],
            index=["NC_002696", "NC_010338", "NC_011916", "NC_014100"],
        )

    def test_deltafile_import(self):
        """Parse test NUCmer .delta/.filter file."""
        result = anim.parse_delta(self.deltafile)
        self.assertEqual(result, (4074001, 2191))

    def test_process_deltadir(self):
        """Process test directory of .delta files into ANIResults."""
        seqfiles = pyani_files.get_fasta_files(self.seqdir)
        orglengths = pyani_files.get_sequence_lengths(seqfiles)
        result = anim.process_deltadir(self.deltadir, orglengths)
        assert_frame_equal(
            result.percentage_identity.sort_index(1).sort_index(),
            self.df_pid.sort_index(1).sort_index(),
        )
