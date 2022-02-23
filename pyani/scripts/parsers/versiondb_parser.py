"""Provides parser for versiondb subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani import pyani_config

from pyani.scripts import subcommands


def get_dry_run_args():
    pass


def build(
    subps: _SubParsersAction, parents: Optional[List[ArgumentParser]] = None
) -> None:
    """Return a command-line parser for the versiondb subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    """
    parser = subps.add_parser(
        "versiondb",
        parents=parents,
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="One of --upgrade, --downgrade, or --dry-run must be specified.",
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
    direction = parser.add_mutually_exclusive_group(required=True)
    direction.add_argument(
        "--upgrade",
        action="store",
        dest="upgrade",
        nargs="?",
        default=None,
        const="head",
        metavar="VERSION",
        help="update an existing database to a newer schema; if no argument is given, 'head' will be used",
    )
    direction.add_argument(
        "--downgrade",
        action="store",
        dest="downgrade",
        nargs="?",
        default=None,
        const="base",
        metavar="VERSION",
        help="revert an existing database to a older schema; if no argument is given, 'base' will be used",
    )
    direction.add_argument(
        "--dry-run",
        action="store",
        dest="dry_run",
        nargs="?",
        metavar="START:END",
        default=None,
        help="produce the SQL that would be run in migrations, without altering the database; start and end versions must be specified",
    )
    parser.add_argument(
        "--alembic_exe",
        action="store",
        dest="alembic_exe",
        required=False,
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
        required=False,
        metavar="NAME",
        help="used to specify an individual database in a multidb setup",
    )
    parser.add_argument(
        "-c",
        "--config",
        action="store",
        dest="alembic_config",
        default=None,
        required=False,
        metavar="FILE",
        help="used to specify a config file for alembic",
    )
    parser.set_defaults(func=subcommands.subcmd_versiondb)
