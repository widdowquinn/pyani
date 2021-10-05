#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2013-2019
# (c) The University of Strathclude 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2013-2019 The James Hutton Institute
# (c) The University of Strathclude 2019-2020
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
"""Test anib subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v

Each command CMD available at the command line as pyani <CMD> is
tested in its own class as a subclass of unittest.TestCase, where
setUp() defines input/output files, a null logger (which is also
picked up by nosetests), and a dictionary of command lines, keyed
by test name, with values representing command-line options.

For each test, command line options are defined in a Namespace and
passed as the sole argument to the appropriate subcommand.
"""

import logging
import os
import unittest

from argparse import Namespace
from collections import namedtuple
from pathlib import Path

import pytest

from pyani.scripts import subcommands


# Convenience struct with paths to third-party executables
ThirdPartyExes = namedtuple("ThirdPartyExes", "blastn_exe format_exe")

# Convenience struct with paths to working directories
DirPaths = namedtuple("DirPaths", "indir outdir")

# Convenience struct for label/class files
LabelPaths = namedtuple("LabelPaths", "classes labels")


@pytest.mark.xfail(reason="ANIb is not currently fully implemented")
class TestANIbsubcommand(unittest.TestCase):

    """Class defining tests of the pyani anib subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.dirpaths = DirPaths(
            Path("tests/test_input/subcmd_anib"), Path("tests/test_output/subcmd_anib")
        )
        os.makedirs(self.dirpaths.outdir, exist_ok=True)
        self.dbpath = Path("tests/test_output/subcmd_createdb/pyanidb")
        self.lblfiles = LabelPaths(
            self.dirpaths.indir / "classes.txt", self.dirpaths.indir / "labels.txt"
        )
        self.exes = ThirdPartyExes("blastn", "makeblastdb")
        self.scheduler = "multiprocessing"

        # Null logger instance
        self.logger = logging.getLogger("TestIndexSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "anib": Namespace(
                indir=self.dirpaths.indir,
                outdir=self.dirpaths.outdir,
                dbpath=self.dbpath,
                force=False,
                name="test_subcmd_anib",
                classes=self.lblfiles.classes,
                labels=self.lblfiles.labels,
                recovery=False,
                cmdline="ANIb test suite",
                blastn_exe=self.exes.blastn_exe,
                format_exe=self.exes.format_exe,
                fragsize=1020,
                scheduler=self.scheduler,
                workers=None,
                disable_tqdm=True,
                jobprefix="ANIbTest",
            )
        }

    def test_anib(self):
        """Test anib run."""
        subcommands.subcmd_anib(self.argsdict["anib"])
