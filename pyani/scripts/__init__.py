# -*- coding: utf-8 -*-
"""Modules to support use of pyani as a script."""

import logging
import shutil

from pathlib import Path

from pyani.pyani_tools import termcolor


class PyaniScriptException(Exception):

    """General exception for pyani.py script."""

    def __init__(self, msg="Error in pyani.py script"):
        """Insantiate exception."""
        Exception.__init__(self, msg)


def make_outdir(outdir: Path, force: bool, noclobber: bool) -> None:
    """Create output directory (allows for force and noclobber).

    :param outdir:  Path, path to output directory
    :param force:  bool, True if an existing directory will be reused
    :param noclobber:  bool, True if existing files are not overwritten

    The intended outcomes are:
    outdir doesn't exist: create outdir
    outdir exists: raise exception
    outdir exists, --force only: remove the directory tree
    outdir exists, --force --noclobber: continue with existing directory tree

    So long as the outdir is created with this function, we need only check
    for args.noclobber elsewhere to see how to proceed when a file exists.
    """
    # Create logger
    logger = logging.getLogger(__name__)
    logger.info("Creating output directory %s", outdir)

    if force:
        logger.warning(termcolor("Output directory overwrite forced", "red"))
        if outdir.is_dir() and noclobber is False:
            logger.warning(termcolor("Clobbering existing directory %s", "red"), outdir)
            shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=force)
