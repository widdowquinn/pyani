#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2018-2019
# (c) University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
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
# Copyright (c) 2018-2019 The James Hutton Institute
# Copyright (c) 2019 University of Strathclyde
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
"""Test command-line parsing for pyani.

This test suite is intended to be run from the repository root using:

pytest -v
"""

import logging

from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from pyani.scripts import pyani_script


class TestCLIParsing(TestCase):

    """Class defining tests of pyani CLI parsing."""

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
        """Set attributes for tests."""
        testdir = Path("tests")
        self.indir = testdir / "test_input" / "sequences"
        self.outdir = testdir / "test_output" / "parsertests"
        self.downloadpath = self.outdir / "downloads"

        self.email = "pyani@pyani.org"
        self.testdbpath = testdir / "test_output" / "parsertests" / "testdb"

        # Lists of command-line arguments for each tests
        self.argsdict = {
            "createdb": ["createdb", "--dbpath", self.testdbpath, "--force"],
            "download": [
                "download",
                "-t",
                "218491",
                "--email",
                self.email,
                self.downloadpath,
                "--force",
            ],
        }

        # Null logger for testing
        self.logger = logging.getLogger("TestCLIParsing logger")
        self.logger.addHandler(logging.NullHandler())

    def test_createdb(self):
        """Create empty test database."""
        pyani_script.run_main(self.argsdict["createdb"], logger=self.logger)

    def test_download(self):
        """Download a single genome."""
        pyani_script.run_main(self.argsdict["download"], logger=self.logger)
