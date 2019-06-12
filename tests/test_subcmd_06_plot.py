#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_subcmd_06_plot.py

Test plot subcommand for pyani

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

from argparse import Namespace
from pathlib import Path

from pyani.scripts import subcommands


class TestPlotSubcommand(unittest.TestCase):
    """Class defining tests of the pyani plot subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.dbpath = Path("tests/test_output/subcmd_createdb/pyanidb")
        self.outdir = Path("tests/test_output/subcmd_plot")
        os.makedirs(self.outdir, exist_ok=True)
        self.run_id = "1"

        # Null logger instance
        self.logger = logging.getLogger("TestPlotSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command-line namespaces
        self.argsdict = {
            "mpl_pdf": Namespace(
                outdir=self.outdir / "mpl",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="pdf",
                method="mpl",
            ),
            "mpl_png": Namespace(
                outdir=self.outdir / "mpl",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="png",
                method="mpl",
            ),
            "mpl_svg": Namespace(
                outdir=self.outdir / "mpl",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="svg",
                method="mpl",
            ),
            "mpl_jpg": Namespace(
                outdir=self.outdir / "mpl",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="jpg",
                method="mpl",
            ),
            "seaborn_pdf": Namespace(
                outdir=self.outdir / "seaborn",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="pdf",
                method="seaborn",
            ),
            "seaborn_png": Namespace(
                outdir=self.outdir / "seaborn",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="png",
                method="seaborn",
            ),
            "seaborn_svg": Namespace(
                outdir=self.outdir / "seaborn",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="svg",
                method="seaborn",
            ),
            "seaborn_jpg": Namespace(
                outdir=self.outdir / "seaborn",
                run_id=self.run_id,
                dbpath=self.dbpath,
                formats="jpg",
                method="seaborn",
            ),
        }

    def test_plot_mpl_pdf(self):
        """test PDF output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["mpl_pdf"], self.logger)

    def test_plot_mpl_png(self):
        """test PNG output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["mpl_png"], self.logger)

    def test_plot_mpl_svg(self):
        """test SVG output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["mpl_svg"], self.logger)

    def test_plot_mpl_jpg(self):
        """test JPG output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["mpl_jpg"], self.logger)

    def test_plot_seaborn_pdf(self):
        """test PDF output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["seaborn_pdf"], self.logger)

    def test_plot_seaborn_png(self):
        """test PNG output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["seaborn_png"], self.logger)

    def test_plot_seaborn_svg(self):
        """test SVG output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["seaborn_svg"], self.logger)

    def test_plot_seaborn_jpg(self):
        """test JPG output from Matplotlib"""
        subcommands.subcmd_plot(self.argsdict["seaborn_jpg"], self.logger)
