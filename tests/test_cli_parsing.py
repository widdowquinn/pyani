#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test_cli_parsing.py


Test command-line parsing for pyani

This test suite is intended to be run from the repository root using:

nosetests -v

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

import logging
import os
import unittest

from pyani.scripts import pyani_script


class TestCLIParsing(unittest.TestCase):
    """Class defining tests of pyani CLI parsing."""

    def setUp(self):
        """Set attributes for tests."""
        self.indir = os.path.join("tests", "test_input", "sequences")
        self.outdir = os.path.join("tests", "test_output", "parsertests")
        self.downloadpath = os.path.join(self.outdir, "downloads")

        self.email = "pyani@pyani.org"
        self.testdbpath = os.path.join("tests", "test_output", "parsertests", "testdb")

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
