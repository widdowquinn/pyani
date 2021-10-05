"""Test aniblastall subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v

Each command CMD available at the command line as pyani <CMD> is
tested in its own class as a subclass of unittest.TestCase, where
setUp() defines input/output files, a null logger (which is also
picked up by nosetests), and a dictionary of command lines, keyed
by test name, with values representing command-line options.

For each test, command line options are defined in a Namespace and
passed as the sole argument to the appropriate subcommand.
"""

import logging
import os
import unittest

from argparse import Namespace
from typing import NamedTuple
from pathlib import Path

import pytest

from pyani.scripts import subcommands


# Convenience struct with paths to third-party executables
class ThirdPartyExes(NamedTuple):
    blastall_exe: Path
    format_exe: Path


# Convenience struct with paths to working directories
class DirPaths(NamedTuple):
    indir: Path
    outdir: Path


# Convenience struct for label/class files
class LabelPaths(NamedTuple):
    classes: Path
    labels: Path


class TestANIblastallSubcommand(unittest.TestCase):

    """Class defining tests of the pyani aniblastall subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.dirpaths = DirPaths(
            Path("tests/test_input/subcmd_anib"),
            Path("tests/test_output/subcmd_aniblastall"),
        )
        os.makedirs(self.dirpaths.outdir, exist_ok=True)
        self.dbpath = Path("tests/test_output/subcmd_createdb/pyanidb")
        self.lblfiles = LabelPaths(
            self.dirpaths.indir / "classes.txt", self.dirpaths.indir / "labels.txt"
        )
        self.exes = ThirdPartyExes("blastall", "formatdb")
        self.scheduler = "multiprocessing"

        # Null logger instance
        self.logger = logging.getLogger("TestIndexSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "aniblastall": Namespace(
                indir=self.dirpaths.indir,
                outdir=self.dirpaths.outdir,
                dbpath=self.dbpath,
                force=False,
                name="test_subcmd_aniblastall",
                classes=self.lblfiles.classes,
                labels=self.lblfiles.labels,
                recovery=False,
                cmdline="ANIblastall test suite",
                blastall_exe=self.exes.blastall_exe,
                formatdb_exe=self.exes.format_exe,
                fragsize=1020,
                scheduler=self.scheduler,
                workers=None,
                disable_tqdm=True,
                jobprefix="ANIblastallTest",
            )
        }

    def test_anib(self):
        """Test aniblastall run."""
        subcommands.subcmd_aniblastall(self.argsdict["aniblastall"])
