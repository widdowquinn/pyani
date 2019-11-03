# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019 University of Strathclyde
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
"""Provides parser for index subcommand."""

from argparse import ArgumentDefaultsHelpFormatter
from pathlib import Path

from pyani.scripts import subcommands


def build(subps, parents=None):
    """Return a command-line parser for the index subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    The index subcommand takes a single positional argument:

    - indir           (directory containing input genome sequence files)
    """
    parser = subps.add_parser(
        "index", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional argument: input directory
    parser.add_argument(
        action="store", dest="indir", default=None, type=Path, help="input directory"
    )
    # Names for output files
    parser.add_argument(
        "--labels",
        dest="labelfname",
        action="store",
        default=Path("labels.txt"),
        type=Path,
        help="Filename for labels file",
    )
    parser.add_argument(
        "--classes",
        dest="classfname",
        action="store",
        default=Path("classes.txt"),
        type=Path,
        help="Filename for classes file",
    )
    parser.set_defaults(func=subcommands.subcmd_index)
