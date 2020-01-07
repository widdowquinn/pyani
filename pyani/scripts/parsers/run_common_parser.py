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
"""Provides parser for arguments common to analysis run subcommands."""

from argparse import ArgumentParser
from pathlib import Path


# Common parser for all subcommands
def build() -> ArgumentParser:
    """Return the common argument parser for analysis run subcommands.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    Common arguments are:

    --name     (human-readable name for the run)
    --labels   (genome labels for this run)
    --classes  (genome classes for this run)
    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "--name",
        dest="name",
        action="store",
        default=None,
        help="Provide a name to help identify the analysis run",
    )
    parser.add_argument(
        "--classes",
        dest="classes",
        action="store",
        default=None,
        type=Path,
        help="Path to genome classes for the run",
    )
    parser.add_argument(
        "--labels",
        dest="labels",
        action="store",
        default=None,
        type=Path,
        help="Path to genome labels for the run",
    )
    parser.add_argument(
        "--recovery",
        dest="recovery",
        action="store_true",
        default=False,
        help="Recover comparison results for this run",
    )
    return parser
