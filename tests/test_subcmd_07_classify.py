#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_07_classify.py

Test classify subcommand for pyani

The test suite is intended to be run from the repository root using:

pytest -v

Individual test classes can be run using, e.g.:

pytest -v tests/test_subcmd_01_download.py::TestDownloadSubcommand
pytest -v tests/test_subcmd_01_download.py::TestDownloadSubcommand::test_download_dry_run

Each command CMD available at the command line as pyani <CMD> is
tested in its own class as a subclass of unittest.TestCase, where
setUp() defines input/output files, a null logger (which is also
picked up by nosetests), and a dictionary of command lines, keyed
by test name, with values representing command-line options.

For each test, command line options are defined in a Namespace and
passed as the sole argument to the appropriate subcommand.

(c) The James Hutton Institute 2018-2019

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2018-2019 The James Hutton Institute

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


import logging
import os
import unittest

import pytest

from argparse import Namespace
from pathlib import Path

from pyani.scripts import subcommands


class TestClassifySubcommand(unittest.TestCase):
    """Class defining tests of the pyani classify subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.dbpath = Path("tests/test_output/subcmd_createdb/pyanidb")
        self.outdir = Path("tests/test_output/subcmd_classify")
        os.makedirs(self.outdir, exist_ok=True)
        self.run_id = "1"

        # Null logger instance
        self.logger = logging.getLogger("TestClassifySubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command-line namespaces
        self.argsdict = {
            "default": Namespace(
                outdir=self.outdir / "default",
                run_id=self.run_id,
                dbpath=self.dbpath,
                cov_min=0.5,
                id_min=0.8,
                min_id=None,
                max_id=None,
                resolution=0.001,
                show_all=False,
                disable_tqdm=True,
            ),
            "no_threshold": Namespace(
                outdir=self.outdir / "no_threshold",
                run_id=self.run_id,
                dbpath=self.dbpath,
                cov_min=0,
                id_min=0,
                min_id=None,
                max_id=None,
                resolution=0.001,
                show_all=False,
                disable_tqdm=True,
            ),
        }

    @pytest.mark.skip_if_exe_missing("nucmer")
    def test_classify_defaults(self):
        """test classification with default settings"""
        subcommands.subcmd_classify(self.argsdict["default"])

    @pytest.mark.skip_if_exe_missing("nucmer")
    def test_classify_no_thresh(self):
        """test classification with no thresholds"""
        subcommands.subcmd_classify(self.argsdict["no_threshold"])
