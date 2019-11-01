#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for pyani package command-line generation.

These tests are intended to be run from the repository root using:

pytest -v

(c) The James Hutton Institute 2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2018 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import unittest

from pathlib import Path

from pyani import anim


class TestNUCmerCmdline(unittest.TestCase):

    """Class defining tests of NUCmer command-line generation."""

    def setUp(self):
        """Set up test parameters."""
        pass

    # One pairwise comparison
    def test_anim_pairwise_basic(self):
        """Test generation of basic NUCmer pairwise comparison command."""
        cmd_nucmer, cmd_filter = anim.construct_nucmer_cmdline(
            Path("file1.fna"), Path("file2.fna")
        )
        tgt_nucmer = "nucmer --mum -p nucmer_output/file1_vs_file2 file1.fna file2.fna"
        tgt_filter = "delta_filter_wrapper.py delta-filter -1 nucmer_output/file1_vs_file2.delta nucmer_output/file1_vs_file2.filter"
        self.assertEqual(cmd_nucmer, tgt_nucmer)
        self.assertEqual(cmd_filter, tgt_filter)

    def test_anim_pairwise_maxmatch(self):
        """Test generation of NUCmer pairwise comparison command with maxmatch."""
        cmd_nucmer, cmd_filter = anim.construct_nucmer_cmdline(
            Path("file1.fna"), Path("file2.fna"), maxmatch=True
        )
        tgt_nucmer = (
            "nucmer --maxmatch -p nucmer_output/file1_vs_file2 file1.fna file2.fna"
        )
        tgt_filter = "delta_filter_wrapper.py delta-filter -1 nucmer_output/file1_vs_file2.delta nucmer_output/file1_vs_file2.filter"
        self.assertEqual(cmd_nucmer, tgt_nucmer)
        self.assertEqual(cmd_filter, tgt_filter)

    # List of pairwise comparisons
    def test_anim_collection(self):
        """Test generation of list of NUCmer comparison commands."""
        files = [Path("file1"), Path("file2"), Path("file3"), Path("file4")]
        cmds_nucmer, cmds_filter = anim.generate_nucmer_commands(files)
        tgts_nucmer = [
            "nucmer --mum -p nucmer_output/file1_vs_file2 file1 file2",
            "nucmer --mum -p nucmer_output/file1_vs_file3 file1 file3",
            "nucmer --mum -p nucmer_output/file1_vs_file4 file1 file4",
            "nucmer --mum -p nucmer_output/file2_vs_file3 file2 file3",
            "nucmer --mum -p nucmer_output/file2_vs_file4 file2 file4",
            "nucmer --mum -p nucmer_output/file3_vs_file4 file3 file4",
        ]
        tgts_filter = [
            "delta_filter_wrapper.py delta-filter -1 nucmer_output/file1_vs_file2.delta nucmer_output/file1_vs_file2.filter",
            "delta_filter_wrapper.py delta-filter -1 nucmer_output/file1_vs_file3.delta nucmer_output/file1_vs_file3.filter",
            "delta_filter_wrapper.py delta-filter -1 nucmer_output/file1_vs_file4.delta nucmer_output/file1_vs_file4.filter",
            "delta_filter_wrapper.py delta-filter -1 nucmer_output/file2_vs_file3.delta nucmer_output/file2_vs_file3.filter",
            "delta_filter_wrapper.py delta-filter -1 nucmer_output/file2_vs_file4.delta nucmer_output/file2_vs_file4.filter",
            "delta_filter_wrapper.py delta-filter -1 nucmer_output/file3_vs_file4.delta nucmer_output/file3_vs_file4.filter",
        ]
        self.assertEqual(cmds_nucmer, tgts_nucmer)
        self.assertEqual(cmds_filter, tgts_filter)
