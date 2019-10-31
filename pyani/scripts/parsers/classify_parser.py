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
"""Provides parser for classify subcommand."""

from argparse import ArgumentDefaultsHelpFormatter

from pyani.scripts import subcommands


def build(subps, parents=None):
    """Return a command-line parser for the classify subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    The classify subcommand takes specific arguments:

    --cov_min       (minimum coverage threshold for an edge)
    --id_min        (minimum identity threshold for an edge)
    --resolution    (number of identity thresholds to test)
    """
    parser = subps.add_parser(
        "classify", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional arguments: output directory and run ID
    parser.add_argument(
        action="store", dest="outdir", default=None, help="output directory"
    )
    parser.add_argument(
        action="store", dest="run_id", default=None, help="run ID to classify"
    )
    # Other optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    # Parameters for classification: minimum %coverage, %identity,
    # and the resolution of thresholds to test
    parser.add_argument(
        "--cov_min",
        dest="cov_min",
        action="store",
        type=float,
        default=0.5,
        help="minimum %%coverage for an edge",
    )
    parser.add_argument(
        "--id_min",
        dest="id_min",
        action="store",
        type=float,
        default=0.8,
        help="minimum %%identity for an edge",
    )
    parser.add_argument(
        "--min_id",
        dest="min_id",
        action="store",
        type=float,
        default=None,
        help="minimum identity threshold to test",
    )
    parser.add_argument(
        "--max_id",
        dest="max_id",
        action="store",
        type=float,
        default=None,
        help="maximum identity threshold to test",
    )
    parser.add_argument(
        "--resolution",
        dest="resolution",
        action="store",
        type=float,
        default=1e-4,
        help="number of identity thresholds to test",
    )
    parser.add_argument(
        "--show_all",
        dest="show_all",
        action="store_true",
        default=False,
        help=(
            "report all intervals in log (default: "
            "only intervals where all subgraphs are k-complete"
        ),
    )
    parser.set_defaults(func=subcommands.subcmd_classify)
