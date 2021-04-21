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
"""Provides parser for anim subcommand."""
import sys

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani import pyani_config
from pyani import anim
from pyani.scripts import subcommands

from pyani import __version__


def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the anim subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited
    """
    parser = subps.add_parser(
        "anim", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional arguments: input and output directories
    parser.add_argument(
        "-i",
        "--indir",
        action="store",
        dest="indir",
        default=None,
        type=Path,
        help="input genome directory",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        dest="outdir",
        default=None,
        type=Path,
        help="output analysis results directory",
    )
    # Optional arguments
    # parser.add_argument(
    #    "--version",
    #    action="version",
    #    version="pyani.py version: {0}§ {1} version: {2}".format(
    #        __version__,
    #        pyani_config.NUCMER_DEFAULT,
    #        anim.get_version(pyani_config.NUCMER_DEFAULT).split("_")[-1],
    #    ).replace('§', '\n\n    '),
    # )
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=Path(".pyani/pyanidb"),
        type=Path,
        help="path to pyani database",
    )
    parser.add_argument(
        "--nucmer_exe",
        dest="nucmer_exe",
        action="store",
        default=pyani_config.NUCMER_DEFAULT,
        type=Path,
        help="path to NUCmer executable",
    )
    parser.add_argument(
        "--filter_exe",
        dest="filter_exe",
        action="store",
        default=pyani_config.FILTER_DEFAULT,
        type=Path,
        help="path to delta-filter executable",
    )
    parser.add_argument(
        "--maxmatch",
        dest="maxmatch",
        action="store_true",
        default=False,
        help="override MUMmer to allow all NUCmer matches",
    )
    parser.add_argument(
        "--nofilter",
        dest="nofilter",
        action="store_true",
        default=False,
        help="do not use delta-filter for 1:1 NUCmer matches",
    )
    parser.set_defaults(func=subcommands.subcmd_anim)
