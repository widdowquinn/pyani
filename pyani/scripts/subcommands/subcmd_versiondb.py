import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import datetime

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
    logger.info(termcolor("Database: %s", "cyan"), str(args.dbpath.resolve()))
    if args.downgrade:
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

    # Create environment variables for alembic to access
    abs_path = args.dbpath.resolve()
    os.environ["PYANI_DATABASE"] = str(abs_path)
    os.environ["ALEMBIC_MIGRATIONS_DIR"] = "alembic"

    # Create a backup of the database
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    # Up/downgrade database
    if args.dry_run:
        logger.info(
            "(Dry run): Migrating database from %s to %s", *args.dry_run.split(":")
        )
        versiondb.migrate_database(args.direction, args, timestamp=timestamp)
    elif args.downgrade:
        logger.info("Downgrading database schema to: %s", args.downgrade)
        shutil.copy(args.dbpath.resolve(), f"{abs_path}.{timestamp}.bak")
        versiondb.migrate_database("downgrade", args)
    elif args.upgrade:
        logger.info("Upgrading database schema to: %s", args.upgrade)
        shutil.copy(args.dbpath.resolve(), f"{abs_path}.{timestamp}.bak")
        versiondb.migrate_database("upgrade", args)

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
