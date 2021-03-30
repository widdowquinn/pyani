"""Provides parser for fastani subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani import pyani_config
from pyani.scripts import subcommands


def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the fastani subcommand.

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
        "--fastani_exe",
        dest="fastani_exe",
        action="store",
        default=pyani_config.FASTANI_DEFAULT,
        type=Path,
        help="path to fastani executable",
    )
    '''parser.add_argument(
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
    )'''
    parser.set_defaults(func=subcommands.subcmd_fastani)
