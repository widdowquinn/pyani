# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD2 5DA,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
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
"""Provides parser for anib subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani import pyani_config
from pyani.scripts import subcommands


def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the anib subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    The terminology may be confusing, but in practice the main parser collects
    command-line arguments that are then available to this parser, which inherits
    options from the parsers in `parents` in addition to those defined below.
    """
    parser = subps.add_parser(
        "anib", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional arguments: input and output directories
    parser.add_argument(
        action="store",
        dest="indir",
        default=None,
        type=Path,
        help="input genome directory",
    )
    parser.add_argument(
        action="store",
        dest="outdir",
        default=None,
        type=Path,
        help="output analysis results directory",
    )
    # Optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=Path(".pyani/pyanidb"),
        type=Path,
        help="path to pyani database",
    )
    parser.add_argument(
        "--blastn_exe",
        dest="blastn_exe",
        action="store",
        default=pyani_config.BLASTN_DEFAULT,
        type=Path,
        help="path to blastn executable",
    )
    parser.add_argument(
        "--format_exe",
        dest="format_exe",
        action="store",
        default=pyani_config.MAKEBLASTDB_DEFAULT,
        type=Path,
        help="path to makeblastdb executable",
    )
    parser.add_argument(
        "--fragsize",
        dest="fragsize",
        action="store",
        type=int,
        default=pyani_config.FRAGSIZE,
        help="blastn query fragment size",
    )
    parser.set_defaults(func=subcommands.subcmd_anib)
