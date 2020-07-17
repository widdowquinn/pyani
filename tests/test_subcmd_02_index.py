#!/usr/bin/env python3
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
"""Test index subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v
"""

import logging
import shutil
import unittest

from argparse import Namespace
from pathlib import Path

from pyani.scripts import subcommands


class TestIndexSubcommand(unittest.TestCase):
    """Class defining tests of the pyani index subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        testdir = Path("tests")
        self.indir = testdir / "test_input" / "subcmd_index"
        self.outdir = testdir / "test_output" / "subcmd_index"
        self.tgtdir = testdir / "test_targets" / "subcmd_index"
        # Remove output directory before copying input data across
        if self.outdir.exists():
            shutil.rmtree(self.outdir)
        shutil.copytree(self.indir, self.outdir)

        # Null logger instance
        self.logger = logging.getLogger("TestIndexSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "index": Namespace(
                indir=self.outdir, labelfname="labels.txt", classfname="classes.txt"
            )
        }

    def test_index(self):
        """Test indexing of downloaded files."""
        # Create index
        print(self.argsdict["index"])
        subcommands.subcmd_index(self.argsdict["index"])

        # Check class and label files
        for fname in ("labels.txt", "classes.txt"):
            with (self.outdir / fname).open() as ofh:
                with (self.tgtdir / fname).open() as tfh:
                    print(ofh.name, tfh.name)
                    self.assertEqual(
                        sorted(ofh.readlines()),
                        sorted(tfh.readlines()),
                        msg="{} files differ".format(fname),
                    )

        #  Check hashes
        for fname in (_ for _ in self.outdir.iterdir() if _.suffix == ".md5"):
            with fname.open() as ofh:
                with (self.tgtdir / fname.name).open() as tfh:
                    self.assertEqual(
                        ofh.read(), tfh.read(), msg="{} files differ".format(fname)
                    )
