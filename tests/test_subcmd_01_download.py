#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
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
# Copyright (c) 2016-2019 The James Hutton Institute
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
"""Test download subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v

Each command CMD available at the command line as pyani <CMD> is
tested in its own class as a subclass of unittest.TestCase, where
setUp() defines input/output files, a null logger (which is also
picked up by pytest), and a dictionary of command lines, keyed
by test name, with values representing command-line options.

For each test, command line options are defined in a Namespace and
passed as the sole argument to the appropriate subcommand.

As the download operations are slow, and subject to network issues,
especially with CI, we mock the download operations.
"""

import logging

from argparse import Namespace
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from pyani.scripts import subcommands


class TestDownloadSubcommand(TestCase):

    """Class defining tests of the pyani download subcommand."""

    @classmethod
    def setup_class(cls):
        """Set up mocking for class."""
        # Mock patcher for downloads
        cls.mock_subcmd_download_patcher = patch(
            "pyani.scripts.subcommands.subcmd_download"
        )
        cls.mock_subcmd_download = cls.mock_subcmd_download_patcher.start()

    @classmethod
    def teardown_class(cls):
        """Close down mocking for class."""
        cls.mock_subcmd_download_patcher.stop()

    def setUp(self):
        """Configure parameters for tests."""
        testdir = Path("tests")
        self.outdir = testdir / "test_output" / "subcmd_download"
        self.krakendir = testdir / "test_output" / "subcmd_download_kraken"
        self.outdir.mkdir(exist_ok=True, parents=True)
        self.krakendir.mkdir(exist_ok=True, parents=True)
        self.api_keypath = Path("~/.ncbi/api_key")

        # Null logger instance
        self.logger = logging.getLogger("TestDownloadSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "dry_run": Namespace(
                outdir=self.outdir,
                taxon="203804",
                email="my.email@my.domain",
                dryrun=True,
                retries=20,
                batchsize=10000,
                timeout=10,
                classfname="classes.txt",
                labelfname="labels.txt",
                noclobber=False,
                force=False,
                kraken=False,
                disable_tqdm=True,
                api_keypath=self.api_keypath,
            ),
            "C_blochmannia": Namespace(
                outdir=self.outdir,
                taxon="203804",
                email="my.email@my.domain",
                retries=20,
                batchsize=10000,
                timeout=10,
                classfname="classes.txt",
                labelfname="labels.txt",
                kraken=False,
                force=True,
                noclobber=False,
                dryrun=False,
                disable_tqdm=True,
                api_keypath=self.api_keypath,
            ),
            "kraken": Namespace(
                outdir=self.krakendir,
                taxon="203804",
                email="my.email@my.domain",
                kraken=True,
                retries=20,
                batchsize=10000,
                timeout=10,
                classfname="classes.txt",
                labelfname="labels.txt",
                force=True,
                noclobber=False,
                dryrun=False,
                disable_tqdm=True,
                api_keypath=self.api_keypath,
            ),
        }

    def test_download_dry_run(self):
        """Dry run of C. blochmannia download."""
        subcommands.subcmd_download(self.argsdict["dry_run"])

    def test_download_c_blochmannia(self):
        """Test C. blochmannia download."""
        subcommands.subcmd_download(self.argsdict["C_blochmannia"])

    def test_download_kraken(self):
        """C. blochmannia download in Kraken format."""
        subcommands.subcmd_download(self.argsdict["kraken"])
