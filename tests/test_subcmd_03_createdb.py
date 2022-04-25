#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# 161 Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2019-2020 University of Strathclyde
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
"""Test createdb subcommand for pyani.

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
import unittest

from argparse import Namespace
from pathlib import Path

from pyani.scripts import subcommands


class TestCreatedbSubcommand(unittest.TestCase):

    """Class defining tests of the pyani createdb subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.dbpath = Path("tests/test_output/subcmd_createdb/pyanidb")

        # Null logger instance
        self.logger = logging.getLogger("TestIndexSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "create_newdb": Namespace(dbpath=self.dbpath, force=False),
            "create_newdb_force": Namespace(dbpath=self.dbpath, force=True),
        }

    def test_create_newdb(self):
        """Test creation of new empty pyani database."""
        # Remove existing dbpath first
        self.dbpath.unlink(missing_ok=True)

        # Create new database
        subcommands.subcmd_createdb(self.argsdict["create_newdb"])

    def test_create_newdb_force(self):
        """Test creation of new pyani database, overwriting old."""
        # Ensure dbpath exists
        assert self.dbpath.exists()

        # Create new database
        subcommands.subcmd_createdb(self.argsdict["create_newdb_force"])
