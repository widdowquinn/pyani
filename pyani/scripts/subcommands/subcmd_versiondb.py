import logging
import os
import platform
import re
import shutil
import subprocess
import sys

from logging import Logger
from pathlib import Path
from typing import List

from pyani import pyani_config

from argparse import Namespace

from pyani import (
    pyani_orm,
    versiondb,
)

from pyani.pyani_tools import termcolor


def subcmd_versiondb(args: Namespace) -> int:
    """Up/downgrade a pyani database.

    :param args:  Namespace, command-line arguments
    :param logger:  logging object
    """
    # Create logger
    logger = logging.getLogger(__name__)

    # Announce what's happening
    if args.upgrade:
        logger.info(termcolor("Downgrading database to %s", bold=True), args.downgrade)
    else:
        logger.info(termcolor("Upgrading database to %s", bold=True), args.upgrade)

    # Get current alembic version
    alembic_version = versiondb.get_version(args.alembic_exe)
    logger.info(termcolor("Alembic version: %s", "cyan"), alembic_version)

    # If the database doesn't exist, raise an error
    if not args.dbpath.is_file():
        logger.error("Database %s does not exist (exiting)", args.dbpath)
        raise SystemError(1)

    # Up/downgrade database
    if args.downgrade:
        logger.info("Downgrading database schema to: %s", args.downgrade)
        versiondb.downgrade_database(args)
    elif args.upgrade:
        logger.info("Upgrading database schema to: %s", args.upgrade)
        versiondb.upgrade_database(args)

    # Do some stuff in alembic

    return 0


# Valid SQLite URL forms are:
# sqlite:///:memory: (or, sqlite://)
# sqlite:///relative/path/to/file.db
# sqlite:////absolute/path/to/file.db

# alembic init --package <dir>  # alembic
# need to change location of sqlalchemy.url

# - alembic.ini
# - alembic
#   - __init__.py
#   - env.py
#   - README
#   - script.py.mako
#   - versions
#       - __init__.py
