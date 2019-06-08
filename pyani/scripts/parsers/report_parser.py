# -*- coding: utf-8 -*-
"""report_parser.py

Provides parser for report subcommand

(c) The James Hutton Institute 2016-2019
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

Copyright (c) 2016-2019 The James Hutton Institute

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

from argparse import ArgumentDefaultsHelpFormatter

from pyani.scripts import subcommands


def build(subps, parents=None):
    """Return a command-line parser for the report subcommand."""
    parser = subps.add_parser(
        "report", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional argument: output directory
    parser.add_argument(
        action="store",
        dest="outdir",
        default=None,
        help="output analysis results directory",
    )
    # Optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    parser.add_argument(
        "--runs",
        action="store_true",
        dest="show_runs",
        default=False,
        help="Report table of analysis runs in database",
    )
    parser.add_argument(
        "--genomes",
        action="store_true",
        dest="show_genomes",
        default=False,
        help="Report table of genomes in database",
    )
    parser.add_argument(
        "--runs_genomes",
        action="store_true",
        dest="show_runs_genomes",
        default=False,
        help="Report table of all genomes for each run in " + "database",
    )
    parser.add_argument(
        "--genomes_runs",
        action="store_true",
        dest="show_genomes_runs",
        default=False,
        help="Report table of all runs in which each genome "
        + "in the database participates",
    )
    parser.add_argument(
        "--run_results",
        action="store",
        dest="run_results",
        default=False,
        help="Report table of results for a pyani run",
    )
    parser.add_argument(
        "--run_matrices",
        action="store",
        dest="run_matrices",
        default=False,
        help="Report matrices of results for a pyani run",
    )
    parser.add_argument(
        "--formats",
        dest="formats",
        action="store",
        default=None,
        help="Output formats (in addition to .tab)",
    )
    parser.set_defaults(func=subcommands.subcmd_report)
