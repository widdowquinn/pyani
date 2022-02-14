import logging
import unittest

from argparse import Namespace
from typing import NamedTuple
from pathlib import Path

from pyani.scripts import subcommands


# Convenience struct with paths to third-party executables
class ThirdPartyExes(NamedTuple):
    alembic_exe: Path


# Convenience struct with paths to working directories
class DirPaths(NamedTuple):
    indir: Path
    outdir: Path


class TestVersiondbSubcommand(unittest.TestCase):

    """Class defining tests of the pyani versiondb subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        testdir = Path("tests")
        self.dirpaths = DirPaths(
            testdir / "test_input" / "subcmd_versiondb",
            testdir / "test_output" / "subcmd_versiondb",
        )
        self.dirpaths.outdir.Mkdir(exist_ok=True)
        self.dbpath = testdir / "test_input" / "subcmd_versiondb" / "pyanidb"
        self.exes = ThirdPartyExes("alembic")

        # Null logger instance
        self.logger = logging.getLogger("TestVersiondbSubcommand logger")
        self.logger.addHandler(logging.NullHandler())

        # Command line namespaces
        self.argsdict = {
            "versiondb": Namespace(
                indir=self.dirpaths.indir,
                outdir=self.dirpaths.outdir,
                dbpath=self.dbpath,
                name="test_versiondb",
                cmdline="alembic test suite",
                alembic_exe=self.exes.alembic_exe,
                workers=None,
                jobprefix="alembicTest",
            )
        }

        def test_versiondb(self):
            """Test versiondb run."""
            print(self.argsdict["versiondb"])
            subcommands.subcmd_versiondb(self.argsdict["versiondb"])
