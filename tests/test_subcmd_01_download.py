#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_01_download.py

Test download subcommand for pyani

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

import pytest

from pyani.scripts import subcommands


class TestDownloadSubcommand(unittest.TestCase):
    """Class defining tests of the pyani download subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.outdir = os.path.join("tests", "test_output", "subcmd_download")
        os.makedirs(self.outdir, exist_ok=True)
        self.krakendir = os.path.join("tests", "test_output", "subcmd_download_kraken")
        os.makedirs(self.krakendir, exist_ok=True)

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
            ),
        }

    def test_download_dry_run(self):
        """dry run of C. blochmannia download"""
        subcommands.subcmd_download(self.argsdict["dry_run"], self.logger)

    def test_download_c_blochmannia(self):
        """test C. blochmannia download"""
        subcommands.subcmd_download(self.argsdict["C_blochmannia"], self.logger)

    @pytest.mark.skip(reason="Investigating Travis-CI test failures")
    def test_download_kraken(self):
        """C. blochmannia download in Kraken format"""
        subcommands.subcmd_download(self.argsdict["kraken"], self.logger)
