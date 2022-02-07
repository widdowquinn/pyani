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

from pyani import pyani_orm


def subcmd_versiondb(args: Namespace) -> int:
    """Create an empty pyani database.

    :param args:  Namespace, command-line arguments
    :param logger:  logging object
    """
    # Create logger
    logger = logging.getLogger(__name__)

    # If the database exists, raise an error rather than overwrite
    if not args.dbpath.is_file():
        logger.error("Database %s does not exist (exiting)", args.dbpath)
        raise SystemError(1)

    # If the path to the database doesn't exist, create it
    if args.upgrade:
        logger.info("Upgrading database schema to: %s", args.upgrade)
    elif args.downgrade:
        logger.info("Downgrading database schema to: %s", args.downgrade)

    # Do some stuff in alembic

    return 0


def get_version(alembic_exe: Path = pyani_config.ALEMBIC_DEFAULT) -> str:
    """Return ALembic package version as a string.

    :param alembic_exe:  path to Alembic executable

    We expect Alembic to return a string on STDOUT as

    .. code-block:: bash

        $ alembic --version
        alembic 1.7.5

    we concatenate this with the OS name.

    The following circumstances are explicitly reported as strings:

    - no executable at passed path
    - non-executable file at passed path (this includes cases where the user doesn't have execute permissions on the file)
    - no version info returned
    """

    try:
        alembic_path = Path(shutil.which(alembic_exe))  # type:ignore
    except TypeError:
        return f"{alembic_exe} is not found in $PATH"

    if not alembic_path.is_file():  # no executable
        return f"No alembic at {alembic_path}"

    # This should catch cases when the file can't be executed by the user
    if not os.access(alembic_path, os.X_OK):  # file exists but not executable
        return f"alembic exists at {alembic_path} but not executable"

    cmdline = [alembic_exe, "--version"]  # type: List
    result = subprocess.run(
        cmdline,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )

    if result.stdout:
        match = re.search(r"(?<=alembic\s)[0-9\.]*", str(result.stdout, "utf-8"))

    version = match.group()  # type: ignore

    if 0 == len(version.strip()):
        return f"alembic exists at {alembic_path} but could not retrieve version"

    return f"{platform.system()}_{version} ({alembic_path})"
