#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_06_anim.py

Test anib subcommand for pyani

The test suite is intended to be run from the repository root using:

pytest -v

Each command CMD available at the command line as pyani <CMD> is
tested in its own class as a subclass of unittest.TestCase, where
setUp() defines input/output files, a null logger (which is also
picked up by nosetests), and a dictionary of command lines, keyed
by test name, with values representing command-line options.

For each test, command line options are defined in a Namespace and
passed as the sole argument to the appropriate subcommand.

(c) The James Hutton Institute 2018

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

import logging
import os
import unittest

from argparse import Namespace

from pyani.scripts import subcommands


class TestANIbSubcommand(unittest.TestCase):
    """Class defining tests of the pyani anib subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.indir = os.path.join("tests", "test_input", "subcmd_anib")
        self.outdir = os.path.join("tests", "test_output", "subcmd_anib")
        os.makedirs(self.outdir, exist_ok=True)
        self.dbpath = os.path.join("tests", "test_output", "subcmd_createdb", "pyanidb")
        self.classes = os.path.join(self.indir, "classes.txt")
        self.labels = os.path.join(self.indir, "labels.txt")
        self.blastn_exe = "blastn"
        self.makeblastdb_exe = "makeblastdb"
        self.scheduler = "multiprocessing"

        # Null logger instance
        self.logger = logging.getLogger("TestANIbSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "anib": Namespace(
                indir=self.indir,
                outdir=self.outdir,
                dbpath=self.dbpath,
                force=False,
                name="test_anib",
                classes=self.classes,
                labels=self.labels,
                recovery=False,
                cmdline="ANIb test suite",
                blastn_exe=self.blastn_exe,
                scheduler=self.scheduler,
                workers=None,
                disable_tqdm=True,
                jobprefix="ANIbTest",
                fragsize=1000,
                makeblastdb_exe=self.makeblastdb_exe,
            )
        }

    def test_anib(self):
        """test ANIb run"""
        subcommands.subcmd_anib(self.argsdict["anib"], self.logger)
