#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""test_subcmd_02_index.py

Test index subcommand for pyani

The test suite is intended to be run from the repository root using:

nosetests -v

Individual test classes can be run using, e.g.:

nosetests -v tests/test_subcmd_01_download.py:TestDownloadSubcommand

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
import shutil
import unittest

from argparse import Namespace

from nose.tools import assert_equal

from pyani.scripts import subcommands


class TestIndexSubcommand(unittest.TestCase):
    """Class defining tests of the pyani index subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.indir = os.path.join("tests", "test_input", "subcmd_index")
        self.outdir = os.path.join("tests", "test_output", "subcmd_index")
        self.tgtdir = os.path.join("tests", "test_targets", "subcmd_index")
        # Remove output directory before copying input data across
        shutil.rmtree(self.outdir)
        shutil.copytree(self.indir, self.outdir)

        # Null logger instance
        self.logger = logging.getLogger("TestIndexSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "index": Namespace(
                indir=self.outdir,
                labelfname=os.path.join("labels.txt"),
                classfname=os.path.join("classes.txt"),
            )
        }

    def test_index(self):
        """test indexing of downloaded files"""
        # Create index
        subcommands.subcmd_index(self.argsdict["index"], self.logger)

        # Check class and label files
        for fname in ("labels.txt", "classes.txt"):
            with open(os.path.join(self.outdir, fname)) as ofh:
                with open(os.path.join(self.tgtdir, fname)) as tfh:
                    assert_equal(
                        ofh.read(), tfh.read(), msg="{} files differ".format(fname)
                    )

        #  Check hashes
        for fname in (
            _ for _ in os.listdir(self.outdir) if os.path.splitext(_)[-1] == ".md5"
        ):
            with open(os.path.join(self.outdir, fname)) as ofh:
                with open(os.path.join(self.tgtdir, fname)) as tfh:
                    assert_equal(
                        ofh.read(), tfh.read(), msg="{} files differ".format(fname)
                    )
