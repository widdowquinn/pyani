"""Provides parser for the compare subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani.scripts import subcommands


def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the plot subcommand.

    :param subps:  collectin of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    The compare subcommand takes specific arguments:

    """
    parser = subps.add_parser(
        "compare", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required arguments: output directory and run IDs
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        dest="outdir",
        type=Path,
        required=True,
        help="output directory",
    )
    parser.add_argument(
        "--run_a",
        action="store",
        dest="run_a",
        required=True,
        help="ID of first run to compare",
    )
    parser.add_argument(
        "--run_b",
        action="store",
        dest="run_b",
        required=True,
        help="ID of second run to compare",
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
    # Graphics methods and formats
    parser.add_argument(
        "--formats",
        dest="formats",
        action="store",
        default="png",
        help="graphics output format (pdf/png/svg/jpg)",
    )
    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        default="seaborn",
        help="graphics method to use for plotting",
        choices=["seaborn", "mpl", "plotly"],
    )
    parser.set_defaults(func=subcommands.subcmd_compare)
