#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) The University of Strathclude 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# Cathedral Street
# Glasgow
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2018 The James Hutton Institute
# (c) The University of Strathclude 2019
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
"""Test report subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v
"""

import logging
import os
import unittest

from argparse import Namespace

from pyani.scripts import subcommands


class TestReportSubcommand(unittest.TestCase):
    """Class defining tests of the pyani report subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.dbpath = os.path.join("tests", "test_output", "subcmd_createdb", "pyanidb")
        self.outdir = os.path.join("tests", "test_output", "subcmd_report")
        os.makedirs(self.outdir, exist_ok=True)

        # Null logger instance
        self.logger = logging.getLogger("TestReportSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "runs": Namespace(
                outdir=self.outdir,
                dbpath=self.dbpath,
                show_runs=True,
                show_genomes=False,
                show_runs_genomes=False,
                show_genomes_runs=False,
                run_results=False,
                run_matrices=False,
                force=True,
                formats="html,excel,stdout",
            ),
            "genomes": Namespace(
                outdir=self.outdir,
                dbpath=self.dbpath,
                show_runs=False,
                show_genomes=True,
                show_runs_genomes=False,
                show_genomes_runs=False,
                run_results=False,
                run_matrices=False,
                force=True,
                formats="html,excel,stdout",
            ),
            "runs_genomes": Namespace(
                outdir=self.outdir,
                dbpath=self.dbpath,
                show_runs=False,
                show_genomes=False,
                show_runs_genomes=True,
                show_genomes_runs=False,
                run_results=False,
                run_matrices=False,
                force=True,
                formats="html,excel,stdout",
            ),
            "genomes_runs": Namespace(
                outdir=self.outdir,
                dbpath=self.dbpath,
                show_runs=False,
                show_genomes=True,
                show_runs_genomes=False,
                show_genomes_runs=True,
                run_results=False,
                run_matrices=False,
                force=True,
                formats="html,excel,stdout",
            ),
            "run_results": Namespace(
                outdir=self.outdir,
                dbpath=self.dbpath,
                show_runs=False,
                show_genomes=False,
                show_runs_genomes=False,
                show_genomes_runs=False,
                run_results="1",
                run_matrices=False,
                force=True,
                formats="html,excel,stdout",
            ),
            "run_matrices": Namespace(
                outdir=self.outdir,
                dbpath=self.dbpath,
                show_runs=False,
                show_genomes=True,
                show_runs_genomes=False,
                show_genomes_runs=False,
                run_results=False,
                run_matrices="1",
                force=True,
                formats="html,excel",
            ),
        }

    def test_runs(self):
        """Test reporting of runs in the database."""
        subcommands.subcmd_report(self.argsdict["runs"], self.logger)

    def test_genomes(self):
        """Test reporting of genomes in the database."""
        subcommands.subcmd_report(self.argsdict["genomes"], self.logger)

    def test_runs_genomes(self):
        """Test reporting of runs/genomes in the database."""
        subcommands.subcmd_report(self.argsdict["runs_genomes"], self.logger)

    def test_genomes_runs(self):
        """Test reporting of genomes_runs in the database."""
        subcommands.subcmd_report(self.argsdict["genomes_runs"], self.logger)

    def test_results(self):
        """Test reporting of run results in the database."""
        subcommands.subcmd_report(self.argsdict["run_results"], self.logger)

    def test_matrices(self):
        """Test reporting of run matrices in the database."""
        subcommands.subcmd_report(self.argsdict["run_matrices"], self.logger)
