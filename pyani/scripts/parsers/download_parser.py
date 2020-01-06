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
"""Provides parser for download subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani.scripts import subcommands


# Subcommand parsers
def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the download subcommand.

    :param subps:  ArgumentParser.subparser
    :param parents:  additional Parser objects

    The download subcommand takes specific arguments:

    -t, --taxon (NCBI taxonomy IDs - comma-separated list, or one ID)
    --email     (email for providing to Entrez services)
    --api_key   (path to file containing personal API key for Entrez)
    --retries   (number of Entrez retry attempts to make)
    --batchsize (number of Entrez records to download in a batch)
    --timeout   (how long to wait for Entrez query timeout)
    -f, --force (allow existing directory overwrite)
    --noclobber (don't replace existing files)
    --labels    (path to write labels file)
    --classes   (path to write classes file)
    """
    parser = subps.add_parser(
        "download", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional argument: output directory
    parser.add_argument(
        action="store", dest="outdir", default=None, type=Path, help="output directory"
    )
    # Required arguments for NCBI query
    parser.add_argument(
        "-t",
        "--taxon",
        dest="taxon",
        action="store",
        default=None,
        help="NCBI taxonomy IDs (required, " + "comma-separated list)",
        required=True,
    )
    parser.add_argument(
        "--email",
        dest="email",
        action="store",
        default=None,
        help="Email associated with NCBI queries (required)",
        required=True,
    )
    # Optional argument for faster/more robust NCBI query
    parser.add_argument(
        "--api_key",
        dest="api_keypath",
        action="store",
        default="~/.ncbi/api_key",
        type=Path,
        help="Path to NCBI Entrez API key file",
    )
    # Arguments controlling connection to NCBI for download
    parser.add_argument(
        "--retries",
        dest="retries",
        action="store",
        default=20,
        help="Number of Entrez retry attempts per request",
    )
    parser.add_argument(
        "--batchsize",
        dest="batchsize",
        action="store",
        default=10000,
        help="Entrez record return batch size",
    )
    parser.add_argument(
        "--timeout",
        dest="timeout",
        action="store",
        default=10,
        help="Timeout for URL connection (s)",
    )
    # Arguments controlling local filehandling
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Allow download to existing directory",
    )
    parser.add_argument(
        "--noclobber",
        dest="noclobber",
        action="store_true",
        default=False,
        help="Don't replace existing files",
    )
    # Names for output files
    parser.add_argument(
        "--labels",
        dest="labelfname",
        action="store",
        default="labels.txt",
        help="Filename for labels file",
    )
    parser.add_argument(
        "--classes",
        dest="classfname",
        action="store",
        default="classes.txt",
        help="Filename for classes file",
    )
    # Output for Kraken
    parser.add_argument(
        "--kraken",
        dest="kraken",
        action="store_true",
        default=False,
        help="Modify downloaded sequence ID for Kraken",
    )
    # Dry-run: do everything except download
    parser.add_argument(
        "--dry-run",
        dest="dryrun",
        action="store_true",
        default=False,
        help="Dry run only, do not download or overwrite.",
    )
    parser.set_defaults(func=subcommands.subcmd_download)
