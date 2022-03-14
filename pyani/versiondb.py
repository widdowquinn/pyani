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
        opts.extend(["-c", args.alembic_config])
    return opts


def construct_alembic_cmdline(
    direction,
    args: Namespace,
    alembic_exe=pyani_config.ALEMBIC_DEFAULT,
):
    if args.dry_run:
        return [
            str(alembic_exe),
            direction,
            args.dry_run,
            "--sql",
            *get_optional_args(args),
        ]  # FAILED: downgrade with --sql requires <fromrev>:<torev>
    elif direction == "upgrade":
        return [str(alembic_exe), direction, args.upgrade, *get_optional_args(args)]
    elif direction == "downgrade":
        return [str(alembic_exe), direction, args.downgrade, *get_optional_args(args)]


def log_output_and_errors(result, direction, args: Namespace, timestamp=None):
    logger = logging.getLogger(__name__)
    if result.stdout:
        logger.info("Alembic stdout:")
        for line in str(result.stdout, "utf-8").split("\n"):
            if line:
                logger.info(line)
        if args.dry_run:
            abs_path = args.dbpath.resolve()
            with open(f"{abs_path}.{direction}.{timestamp}.sql", "w") as sqlfile:
                for line in str(result.stdout, "utf-8").split("\n"):
                    if line:
                        sqlfile.write(f"{line}\n")

    if result.stderr:
        logger.info("Alembic stderr:")
        for line in str(result.stderr, "utf-8").split("\n"):
            if line:
                logger.info(line)


def migrate_database(direction, args: Namespace, timestamp=None):
    cmdline = construct_alembic_cmdline(direction, args)

    result = subprocess.run(
        cmdline,
        shell=False,
        capture_output=True,
    )

    log_output_and_errors(result, direction, args, timestamp)
