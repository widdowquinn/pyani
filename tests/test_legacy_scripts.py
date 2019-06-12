#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_legacy_scripts.py

Test legacy average_nucleotide_identity.py and
genbank_download_genomes_by_taxon.py scripts

The test suite is intended to be run from the repository root using:

pytest -v

Each command CMD available at the command line as pyani <CMD> is
tested in its own class as a subclass of unittest.TestCase, where
setUp() defines input/output files, a null logger (which is also
picked up by nosetests), and a dictionary of command lines, keyed
by test name, with values representing command-line options.

For each test, command line options are defined in a Namespace and
passed as the sole argument to the appropriate subcommand.

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
import unittest

from argparse import Namespace
from collections import namedtuple
from pathlib import Path

import pytest

from pyani.scripts import average_nucleotide_identity, genbank_get_genomes_by_taxon


# Convenience struct for executables
Executables = namedtuple(
    "Executables", "nucmer filter blastn blastall makeblastdb formatdb"
)


class TestLegacyScripts(unittest.TestCase):
    """Class defining tests of the pyani download subcommand."""

    def setUp(self):
        """Configure parameters for tests."""
        self.outdir = Path("tests/test_output/legacy_scripts")
        self.dldir = self.outdir / "C_blochmannia"
        os.makedirs(self.outdir, exist_ok=True)
        self.exes = Executables(
            "nucmer", "delta-filter", "blastn", "blastall", "makeblastdb", "formatdb"
        )

        # Null logger instance
        self.logger = logging.getLogger("TestLegacyScripts logger")
        self.logger.addHandler(logging.NullHandler())

        # Command-line namespaces
        self.argsdict = {
            "download": Namespace(
                outdirname=self.dldir,
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
            ),
            "anim_seaborn": Namespace(
                outdirname=self.outdir / "ANIm_output",
                indirname=self.dldir,
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
            ),
        }

    @pytest.mark.run(order=1)
    def test_legacy_genome_downloads(self):
        """Uue legacy script to download genomes"""
        genbank_get_genomes_by_taxon.run_main(self.argsdict["download"], self.logger)

    @pytest.mark.run(order=2)
    def test_legacy_anim_seaborn(self):
        """Uue legacy script to run ANIm"""
        average_nucleotide_identity.run_main(self.argsdict["anim_seaborn"], self.logger)
