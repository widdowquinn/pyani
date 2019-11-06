#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2018-2019
# (c) The University of Strathclyde 2019
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
# (c) The James Hutton Institute 2018-2019
# (c) The University of Strathclyde 2019
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
"""Test plot subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v
"""

import logging
import os
import unittest

from argparse import Namespace
from pathlib import Path

import matplotlib.pyplot as plt

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
        """Test PDF output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["mpl_pdf"], self.logger)
        plt.close("all")

    def test_plot_mpl_png(self):
        """Test PNG output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["mpl_png"], self.logger)
        plt.close("all")

    def test_plot_mpl_svg(self):
        """Test SVG output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["mpl_svg"], self.logger)
        plt.close("all")

    def test_plot_mpl_jpg(self):
        """Test JPG output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["mpl_jpg"], self.logger)
        plt.close("all")

    def test_plot_seaborn_pdf(self):
        """Test PDF output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["seaborn_pdf"], self.logger)
        plt.close("all")

    def test_plot_seaborn_png(self):
        """Test PNG output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["seaborn_png"], self.logger)
        plt.close("all")

    def test_plot_seaborn_svg(self):
        """Test SVG output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["seaborn_svg"], self.logger)
        plt.close("all")

    def test_plot_seaborn_jpg(self):
        """Test JPG output from Matplotlib."""
        subcommands.subcmd_plot(self.argsdict["seaborn_jpg"], self.logger)
        plt.close("all")
