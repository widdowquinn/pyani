#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_04_anim.py

Test anim subcommand for pyani

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

Copyright (c) 2017 The James Hutton Institute

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


class TestANImSubcommand(unittest.TestCase):
    """Class defining tests of the pyani anim subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.indir = os.path.join("tests", "test_input", "subcmd_anim")
        self.outdir = os.path.join("tests", "test_output", "subcmd_anim")
        os.makedirs(self.outdir, exist_ok=True)
        self.dbpath = os.path.join("tests", "test_output", "subcmd_createdb", "pyanidb")
        self.classes = os.path.join(self.indir, "classes.txt")
        self.labels = os.path.join(self.indir, "labels.txt")
        self.nucmer_exe = "nucmer"
        self.filter_exe = "delta-filter"
        self.scheduler = "multiprocessing"

        # Null logger instance
        self.logger = logging.getLogger("TestIndexSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "anim": Namespace(
                indir=self.indir,
                outdir=self.outdir,
                dbpath=self.dbpath,
                force=False,
                name="test_anim",
                classes=self.classes,
                labels=self.labels,
                recovery=False,
                cmdline="ANIm test suite",
                nucmer_exe=self.nucmer_exe,
                filter_exe=self.filter_exe,
                maxmatch=False,
                nofilter=False,
                scheduler=self.scheduler,
                workers=None,
                disable_tqdm=True,
                jobprefix="ANImTest",
            )
        }

    def test_anim(self):
        """test anim run"""
        subcommands.subcmd_anim(self.argsdict["anim"], self.logger)
