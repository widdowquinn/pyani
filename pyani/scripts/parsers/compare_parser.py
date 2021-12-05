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
    # parser.add_argument(
    #     "--run_a",
    #     action="store",
    #     dest="run_a",
    #     required=True,
    #     help="ID of first run to compare",
    # )
    # parser.add_argument(
    #     "--run_b",
    #     action="store",
    #     dest="run_b",
    #     required=True,
    #     help="ID of second run to compare",
    # )
    parser.add_argument(
        "--ref_ids",
        action="store",
        dest="ref_ids",
        default=None,
        metavar="RUN_ID",
        nargs="+",
        required=True,
        type=int,
        help="Space-separated list of run_ids to use as reference(s) for comparisons",
    )
    parser.add_argument(
        "--run_ids",
        action="store",
        dest="run_ids",
        default=None,
        metavar="RUN_ID",
        nargs="+",
        required=True,
        type=int,
        help="Space-separated list of run_ids to compare to reference(s)",
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
        default=["png"],
        metavar="FORMAT",
        nargs="+",
        choices=["pdf", "png", "svg", "jpg"],
        help="graphics output format; options: (pdf, png, svg, jpg)",
    )
    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        default="seaborn",
        metavar="METHOD",
        nargs="+",
        choices=["seaborn", "mpl", "plotly"],
        help="graphics method to use for plotting; options (seaborn, mpl, plotly)",
    )
    parser.add_argument(
        "--workers",
        dest="workers",
        action="store",
        default=None,
        type=int,
        help="Number of worker processes for multiprocessing "
        "(default zero, meaning use all available cores)",
    )
    parser.set_defaults(func=subcommands.subcmd_compare)
