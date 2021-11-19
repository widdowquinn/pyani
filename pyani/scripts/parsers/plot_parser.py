# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
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
"""Provides parser for plot subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani.scripts import subcommands


def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the plot subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    The plot subcommand takes specific arguments:

    --method        (graphics method to use)
    """
    parser = subps.add_parser(
        "plot", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required arguments: output directory and run ID
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        dest="outdir",
        default=None,
        type=Path,
        help="output directory",
        required=True,
    )
    parser.add_argument(
        "--run_id",
        action="store",
        dest="run_id",
        default=None,
        help="run ID to plot",
        required=True,
    )
    # Other optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=Path(".pyani/pyanidb"),
        type=Path,
        help="path to pyani database",
    )
    # Graphics methods and formats
    parser.add_argument(
        "--formats",
        dest="formats",
        action="store",
        default="png",
        help="graphics output format (pdf/png/svg/jpg)",
    )
    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        default="seaborn",
        help="graphics method to use for plotting",
        choices=["seaborn", "mpl", "plotly"],
    )
    parser.set_defaults(func=subcommands.subcmd_plot)
