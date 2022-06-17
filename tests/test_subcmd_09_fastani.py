"""Test fastANI subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v
"""

import logging
import unittest

from argparse import Namespace
from typing import NamedTuple
from pathlib import Path

from pyani.scripts import subcommands


# Convenience struct with paths to third-party executables
class ThirdPartyExes(NamedTuple):
    fastani_exe: Path


# Convenience struct with paths to working directories
class DirPaths(NamedTuple):
    indir: Path
    outdir: Path


# Convenience struct for label/class files
class LabelPaths(NamedTuple):
    classes: Path
    labels: Path


class TestfastANISubcommand(unittest.TestCase):

    """Class defining tests of the pyani fastANI subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        testdir = Path("tests")
        self.dirpaths = DirPaths(
            testdir / "test_input" / "subcmd_fastani",
            testdir / "test_output" / "subcmd_fastani",
        )
        self.dirpaths.outdir.mkdir(exist_ok=True)
        self.dbpath = testdir / "test_output" / "subcmd_createdb" / "pyanidb"
        self.lblfiles = LabelPaths(
            self.dirpaths.indir / "classes.txt", self.dirpaths.indir / "labels.txt"
        )
        self.exes = ThirdPartyExes("fastANI")
        self.scheduler = "multiprocessing"

        # Null logger instance
        self.logger = logging.getLogger("TestIndexSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "fastani": Namespace(
                indir=self.dirpaths.indir,
                outdir=self.dirpaths.outdir,
                dbpath=self.dbpath,
                force=False,
                name="test_fastani",
                classes=self.lblfiles.classes,
                labels=self.lblfiles.labels,
                recovery=False,
                cmdline="fastANI test suite",
                fastani_exe=self.exes.fastani_exe,
                fragLen=3000,
                kmerSize=16,
                minFraction=0.2,
                scheduler=self.scheduler,
                workers=None,
                disable_tqdm=True,
                jobprefix="fastANITest",
            )
        }

    # @unittest.skip(
    #     "This test currently fails for reasons unknown. \
    #     Something about the comparisons.kmersize column not existing."
    # )
    def test_fastani(self):
        """Test fastani run."""
        print(self.argsdict["fastani"])
        subcommands.subcmd_fastani(self.argsdict["fastani"])
