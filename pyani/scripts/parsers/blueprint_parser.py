# -*- coding: utf-8 -*-
# (c) The University of Strathclyde 2021–Present
# Author: Bailey Harrington
#
# Contact: bailey.harrington@strath.ac.uk

"""Parser for example subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani import pyani_config
from pyani.scripts import subcommands


# Function that adds all of the arguments accepted by the parser
# along with default values, help information, and other bits
# specific to this submodule parser
def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the exaple subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    The terminology may be confusing, but in practice the main parser collects
    command-line arguments that are then available to this parser, which inherits
    options from the parsers in `parents` in addition to those defined below.
    """
    parser = subps.add_parser(
        "blueprint", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required arguments: input and output directories
    parser.add_argument(
        "-i",
        "--indir",
        required=True,
        action="store",
        dest="indir",
        default=None,
        type=Path,
        help="input directory",
    )
    parser.add_argument(  # I think this should be handled with pyani's overall output
        "-o",
        "--outdir",
        required=True,
        action="store",
        dest="outdir",
        default=None,
        type=Path,
        help="output analysis results directory",
    )
    # Optional arguments
    parser.add_argument(  # keep ;)
        "--dbpath",
        action="store",
        dest="dbpath",
        default=Path(".pyani/pyanidb"),
        type=Path,
        help="path to pyani database",
    )
    parser.add_argument(
        "--blueprint_exe",
        dest="blueprint_exe",
        action="store",
        default=pyani_config.BLUEPRINT_DEFAULT,
        type=Path,
        help="path to blueprint executable",
    )
    parser.set_defaults(func=subcommands.subcmd_blueprint)
