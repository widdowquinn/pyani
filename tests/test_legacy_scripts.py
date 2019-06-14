#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_legacy_scripts.py

Test legacy average_nucleotide_identity.py and
genbank_download_genomes_by_taxon.py scripts

The test suite is intended to be run from the repository root using:

pytest -v

The two legacy scripts download genomes, then carry out ANI analysis. The
pytest ordering plug-in is used to guarantee that the download script
tests are conducted first.

(c) The James Hutton Institute 2019

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2019 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import logging
import os
import shutil

from argparse import Namespace
from collections import namedtuple
from pathlib import Path

import pytest

from pyani.scripts import average_nucleotide_identity, genbank_get_genomes_by_taxon

from tools import modify_namespace, PyaniTestCase

# Convenience struct for executables
Executables = namedtuple(
    "Executables", "nucmer filter blastn blastall makeblastdb formatdb"
)

# Convenience struct for test directories
TestDirs = namedtuple("TestDirs", "outdir tgtdir dldir")


class TestLegacyScripts(PyaniTestCase):
    """Class defining tests of the pyani download subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.testdirs = TestDirs(
            Path("tests/test_output/legacy_scripts"),
            Path("tests/test_targets/legacy_scripts"),
            Path("tests/test_output/legacy_scripts/C_blochmannia"),
        )
        os.makedirs(self.testdirs.outdir, exist_ok=True)
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
            outdirname=self.testdirs.outdir / "ANIm_seaborn",
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
            "anim_seaborn": self.base_ani,
            "anim_mpl": modify_namespace(
                self.base_ani,
                {"outdirname": self.testdirs.outdir / "ANIm_mpl", "gmethod": "mpl"},
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
        """Use legacy script to download genomes"""
        shutil.rmtree(self.testdirs.outdir)  # Clean before running
        genbank_get_genomes_by_taxon.run_main(self.argsdict["download"], self.logger)
        self.assertDirsEqual(
            self.testdirs.dldir, self.testdirs.tgtdir / "C_blochmannia"
        )

    @pytest.mark.run(order=2)
    def test_legacy_anim_seaborn(self):
        """Use legacy script to run ANIm (seaborn output)"""
        args = self.argsdict["anim_seaborn"]
        average_nucleotide_identity.run_main(args, self.logger)
        # self.assertDirsEqual(
        #     args.outdirname,  # pylint: disable=no-member
        #     self.testdirs.tgtdir / args.outdirname.name,  # pylint: disable=no-member
        # )

    @pytest.mark.run(order=2)
    def test_legacy_anim_mpl(self):
        """Use legacy script to run ANIm (mpl output)"""
        args = self.argsdict["anim_mpl"]
        average_nucleotide_identity.run_main(args, self.logger)
        # self.assertDirsEqual(
        #     args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        # )

    @pytest.mark.run(order=2)
    def test_legacy_anib_seaborn(self):
        """Use legacy script to run ANIb (seaborn output)"""
        args = self.argsdict["anib_seaborn"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )

    @pytest.mark.run(order=2)
    def test_legacy_anib_mpl(self):
        """Use legacy script to run ANIb (mpl output)"""
        args = self.argsdict["anib_mpl"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )

    @pytest.mark.run(order=2)
    def test_legacy_tetra_seaborn(self):
        """Use legacy script to run TETRA (seaborn output)"""
        args = self.argsdict["tetra_seaborn"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )

    @pytest.mark.run(order=2)
    def test_legacy_tetra_mpl(self):
        """Use legacy script to run TETRA (mpl output)"""
        args = self.argsdict["tetra_mpl"]
        average_nucleotide_identity.run_main(args, self.logger)
        self.assertDirsEqual(
            args.outdirname, self.testdirs.tgtdir / args.outdirname.name
        )
