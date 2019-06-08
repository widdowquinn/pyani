# -*- coding: utf-8 -*-
"""db_parser.py

Provides parser for db subcommand

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
    """Return a command-line parser for the db subcommand."""
    parser = subps.add_parser(
        "db", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        action="store",
        type=int,
        dest="run_id",
        default=None,
        help="work with this runID",
    )
    # Path to database (default: .pyani/pyanidb)
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    parser.add_argument(
        "--relabel",
        action="store",
        dest="relabelfname",
        default=None,
        help="relabel genomes with hashes in passed file",
    )
    parser.add_argument(
        "--reclass",
        action="store",
        dest="reclassfname",
        default=None,
        help="reclass genomes with hashes in passed file",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        dest="force",
        default=False,
        help="force database operation",
    )
    parser.set_defaults(func=subcommands.subcmd_db)
