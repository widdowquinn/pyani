# -*- coding: utf-8 -*-
"""classify_parser.py

Provides parser for classify subcommand

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
    """Return a command-line parser for the classify subcommand.

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
        "--resolution",
        dest="resolution",
        action="store",
        type=int,
        default=1500,
        help="number of identity thresholds to test",
    )
    parser.set_defaults(func=subcommands.subcmd_classify)
