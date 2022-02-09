"""Provides parser for versiondb subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani import pyani_config

from pyani.scripts import subcommands


def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the versiondb subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    """
    parser = subps.add_parser(
        "versiondb", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Path to database (default: .pyani/pyanidb)
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=Path(".pyani/pyanidb"),
        type=Path,
        help="path to pyani database",
    )
    direction = parser.add_mutually_exclusive_group(required=False)
    direction.add_argument(
        "-u",
        "--upgrade",
        action="store",
        dest="upgrade",
        default="head",
        metavar="VERSION",
        help="update an existing database to a newer schema; default is to upgrade to the newest version",
    )
    direction.add_argument(
        "-d",
        "--downgrade",
        action="store",
        dest="downgrade",
        default=None,
        metavar="VERSION",
        help="revert an existing database to a older schema",
    )
    parser.add_argument(
        "--alembic_exe",
        action="store",
        dest="alembic_exe",
        default=pyani_config.ALEMBIC_DEFAULT,
        type=Path,
        help="path to alembic executable",
    )
    parser.add_argument(
        "-n",
        "--name",
        action="store",
        dest="dbname",
        default=None,
        metavar="NAME",
        help="used to specify an individual database in a multidb setup",
    )
    parser.add_argument(
        "-c",
        "--config",
        action="store",
        dest="alembic_config",
        default=None,
        metavar="FILE",
        help="used to specify a config file for alembic",
    )
    parser.set_defaults(func=subcommands.subcmd_versiondb)
