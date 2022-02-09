import logging
import os
import platform
import re
import shutil
import subprocess
import sys

from typing import List

from pathlib import Path

from pyani import pyani_config

from argparse import Namespace


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


def get_optional_args(args: Namespace):
    opts = []
    if args.dbname:
        opts.extend(["-n", args.dbname])
    if args.alembic_config:
        opts.extend(["-c", args.config])
    return opts


def construct_alembic_cmdline(
    direction,
    args: Namespace,
    alembic_exe=pyani_config.ALEMBIC_DEFAULT,
):
    if direction == "upgrade":
        return [alembic_exe, direction, args.upgrade, *get_optional_args(args)]
    elif direction == "downgrade":
        return [alembic_exe, direction, args.downgrade, *get_optional_args(args)]


def upgrade_database(args: Namespace):
    logger = logging.getLogger(__name__)
    cmdline = construct_alembic_cmdline("upgrade", args)
    result = subprocess.run(
        cmdline,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    # with str(result.stdout, "utf-8") as pipe:
    # for line in str(result.stdout, "utf-8"):
    # logger.info('Alembic: %s', str(result.stderr, "utf-8"))
    for line in str(result.stderr, "utf-8").split("\n"):
        if line:
            logger.info("Alembic: %s", line)


def downgrade_database(args: Namespace):
    logger = logging.getLogger(__name__)
    cmdline = construct_alembic_cmdline("downgrade", args)
    result = subprocess.run(
        cmdline,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    # logger.info('A: %s', str(result.stderr, "utf-8"))
    for line in str(result.stderr, "utf-8").split("\n"):
        if line:
            logger.info("Alembic: %s", line)
