#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019 University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""Test legacy scripts.

The test suite is intended to be run from the repository root using:

pytest -v

The two legacy scripts download genomes, then carry out ANI analysis. The
pytest ordering plug-in is used to guarantee that the download script
tests are conducted first.
"""

import logging
import shutil

from argparse import Namespace
from collections import namedtuple
from pathlib import Path

import pytest

from pyani import anim, pyani_tools
from pyani.scripts import average_nucleotide_identity, genbank_get_genomes_by_taxon

from tools import modify_namespace, PyaniTestCase

# Convenience struct for executables
Executables = namedtuple(
    "Executables", "nucmer filter blastn blastall makeblastdb formatdb"
)

# Convenience struct for test directories
TestDirs = namedtuple("TestDirs", "outdir tgtdir dldir")


@pytest.mark.xfail(reason="legacy scripts not yet ready for test")
class TestLegacyScripts(PyaniTestCase):
    """Class defining tests of the pyani download subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        testdir = Path("tests")
        self.testdirs = TestDirs(
            testdir / "test_output" / "legacy_scripts",
            testdir / "test_targets" / "legacy_scripts",
            testdir / "test_output" / "legacy_scripts" / "C_blochmannia",
        )
        self.testdirs.outdir.mkdir(exist_ok=True)
        self.exes = Executables(
            "nucmer", "delta-filter", "blastn", "blastall", "makeblastdb", "formatdb"
        )

        # Base namespaces for each script
        self.base_download = Namespace(
            outdirname=self.testdirs.dldir,
            taxon="203804",
            verbose=False,
            force=True,
            noclobber=False,
            logfile=None,
            format="fasta",
            email="emailme@my.email.domain",
            retries=20,
            batchsize=10000,
            timeout=10,
        )
        self.base_ani = Namespace(
            outdirname=self.testdirs.outdir / f"ANIm_seaborn_{anim.get_version()}",
            indirname=self.testdirs.dldir,
            verbose=False,
            force=True,
            fragsize=1020,
            logfile="test_ANIm.log",
            skip_nucmer=False,
            skip_blastn=False,
            noclobber=False,
            nocompress=False,
            graphics=True,
            gformat="pdf,png",
            gmethod="seaborn",
            labels=None,
            classes=None,
            method="ANIm",
            scheduler="multiprocessing",
            workers=None,
            sgeargs=None,
            sgegroupsize=10000,
            maxmatch=False,
            nucmer_exe=self.exes.nucmer,
            filter_exe=self.exes.filter,
            blastn_exe=self.exes.blastn,
            blastall_exe=self.exes.blastall,
            makeblastdb_exe=self.exes.makeblastdb,
            formatdb_exe=self.exes.formatdb,
            write_excel=False,
            rerender=False,
            subsample=None,
            seed=None,
            jobprefix="ANI",
        )

        # Null logger instance
        self.logger = logging.getLogger("TestLegacyScripts logger")
        self.logger.addHandler(logging.NullHandler())

        # Command-line namespaces
        self.argsdict = {
            "download": self.base_download,
            "anim_seaborn": self.base_ani,  #  Baseline namespace
            "anim_mpl": modify_namespace(
                self.base_ani,
                {
                    "outdirname": self.testdirs.outdir
                    / f"ANIm_mpl_{anim.get_version()}",
                    "gmethod": "mpl",
                },
            ),
            "anib_seaborn": modify_namespace(
                self.base_ani,
                {"outdirname": self.testdirs.outdir / "ANIb_seaborn", "method": "ANIb"},
            ),
            "anib_mpl": modify_namespace(
                self.base_ani,
                {
                    "outdirname": self.testdirs.outdir / "ANIb_mpl",
                    "gmethod": "mpl",
                    "method": "ANIb",
                },
            ),
            "tetra_seaborn": modify_namespace(
                self.base_ani,
                {
                    "outdirname": self.testdirs.outdir / "TETRA_seaborn",
                    "method": "TETRA",
                },
            ),
            "tetra_mpl": modify_namespace(
                self.base_ani,
                {
                    "outdirname": self.testdirs.outdir / "TETRA_mpl",
                    "gmethod": "mpl",
                    "method": "TETRA",
                },
            ),
        }

    @pytest.mark.run(order=1)
    def test_legacy_genome_downloads(self):
        r"""Use legacy script to download genomes.

        genbank_get_genomes_by_taxon.py \
            -o tests/test_output/legacy_scripts/C_blochmannia \
            --email emailme@my.email.domain \
            -t 203804 -f
        """
        if self.testdirs.dldir.is_dir():
            shutil.rmtree(self.testdirs.dldir)  # Clean before running
        genbank_get_genomes_by_taxon.run_main(self.argsdict["download"], self.logger)
        self.assertDirsEqual(
            self.testdirs.dldir, self.testdirs.tgtdir / "C_blochmannia"
        )

    @pytest.mark.run(order=2)
    @pytest.mark.skipif(
        not pyani_tools.has_dependencies().mummer,
        reason="nucmer executable not available",
    )
    def test_legacy_anim_seaborn(self):
        r"""Use legacy script to run ANIm (seaborn output).

        average_nucleotide_identity.py \
            -l test_ANIm.log \
            -i tests/test_output/legacy_scripts/C_blochmannia \
            -o tests/test_output/legacy_scripts/ANIm_seaborn \
            -g --gmethod seaborn --gformat pdf,png \
            -f --jobprefix ANI
        """
        args = self.argsdict["anim_seaborn"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname,  # pylint: disable=no-member
            self.testdirs.tgtdir / args.outdirname.name,  # pylint: disable=no-member
        )

    @pytest.mark.run(order=2)
    @pytest.mark.skipif(
        not pyani_tools.has_dependencies().mummer,
        reason="nucmer executable not available",
    )
    def test_legacy_anim_mpl(self):
        r"""Use legacy script to run ANIm (mpl output).

        average_nucleotide_identity.py \
            -l test_ANIm.log \
            -i tests/test_output/legacy_scripts/C_blochmannia \
            -o tests/test_output/legacy_scripts/ANIm_mpl \
            -g --gmethod mpl --gformat pdf,png \
            -f --jobprefix ANI
        """
        args = self.argsdict["anim_mpl"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )

    @pytest.mark.run(order=2)
    @pytest.mark.skipif(
        not pyani_tools.has_dependencies().blast,
        reason="blastn executable not available",
    )
    def test_legacy_anib_seaborn(self):
        r"""Use legacy script to run ANIb (seaborn output).

        average_nucleotide_identity.py \
            -l test_ANIb.log \
            -i tests/test_output/legacy_scripts/C_blochmannia \
            -o tests/test_output/legacy_scripts/ANIb_seaborn \
            -g --gmethod seaborn --gformat pdf,png \
            -f --jobprefix ANI
        """
        args = self.argsdict["anib_seaborn"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )

    @pytest.mark.run(order=2)
    @pytest.mark.skipif(
        not pyani_tools.has_dependencies().blast,
        reason="blastn executable not available",
    )
    def test_legacy_anib_mpl(self):
        r"""Use legacy script to run ANIb (mpl output).

        average_nucleotide_identity.py \
            -l test_ANIb.log \
            -i tests/test_output/legacy_scripts/C_blochmannia \
            -o tests/test_output/legacy_scripts/ANIb_mpl \
            -g --gmethod mpl --gformat pdf,png \
            -f --jobprefix ANI
        """
        args = self.argsdict["anib_mpl"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )

    @pytest.mark.run(order=2)
    def test_legacy_tetra_seaborn(self):
        r"""Use legacy script to run TETRA (seaborn output)."""
        args = self.argsdict["tetra_seaborn"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )

    @pytest.mark.run(order=2)
    def test_legacy_tetra_mpl(self):
        r"""Use legacy script to run TETRA (mpl output)."""
        args = self.argsdict["tetra_mpl"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )
