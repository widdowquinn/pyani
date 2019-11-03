#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact: leighton.pritchard@strath.ac.uk
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
"""Test anib.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import os
import unittest

from pathlib import Path

import pandas as pd
import pytest

from pandas.util.testing import assert_frame_equal

from pyani import anib, pyani_files, pyani_tools


@pytest.mark.skipif(
    not pyani_tools.has_dependencies().blast, reason="BLASTN executable not found"
)
class TestBLASTCmdline(unittest.TestCase):

    """Class defining tests of BLAST command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        testdir = Path("tests")
        self.indir = testdir / "test_input" / "anib"
        self.outdir = testdir / "test_output" / "anib"
        self.seqdir = testdir / "test_input" / "sequences"
        self.infiles = list(self.seqdir.iterdir())
        self.fraglen = 1000
        self.fmtdboutdir = self.outdir / "formatdb"
        self.fmtdbcmd = (
            f"formatdb -p F -i {self.fmtdboutdir / 'NC_002696.fna'} -t NC_002696"
        )
        self.makeblastdbdir = self.outdir / "makeblastdb"
        self.makeblastdbcmd = f"makeblastdb -dbtype nucl -in {self.seqdir / 'NC_002696.fna'} -title NC_002696 -out {self.makeblastdbdir / 'NC_002696.fna'}"
        self.blastdbfnames = [
            self.seqdir / fname for fname in ("NC_002696.fna", "NC_010338.fna")
        ]
        self.blastdbtgt = [
            (
                f"makeblastdb -dbtype nucl -in {self.seqdir / 'NC_002696.fna'} -title NC_002696 -out {self.outdir / 'NC_002696.fna'}",
                self.outdir / "NC_002696.fna",
            ),
            (
                f"makeblastdb -dbtype nucl -in {self.seqdir / 'NC_010338.fna'} -title NC_010338 -out {self.outdir / 'NC_010338.fna'}",
                self.outdir / "NC_010338.fna",
            ),
        ]
        self.blastdbtgtlegacy = [
            (
                f"formatdb -p F -i {self.outdir / 'NC_002696.fna'} -t NC_002696",
                self.outdir / "NC_002696.fna",
            ),
            (
                f"formatdb -p F -i {self.outdir / 'NC_010338.fna'} -t NC_010338",
                self.outdir / "NC_010338.fna",
            ),
        ]
        self.blastncmd = (
            f"blastn -out {self.outdir / 'NC_002696_vs_NC_010338.blast_tab'} "
            f"-query {self.seqdir / 'NC_002696.fna'} "
            f"-db {self.seqdir / 'NC_010338.fna'} "
            "-xdrop_gap_final 150 -dust no -evalue 1e-15 -max_target_seqs 1 "
            "-outfmt '6 qseqid sseqid length mismatch pident nident qlen slen qstart qend sstart send positive ppos gaps' "
            "-task blastn"
        )
        self.blastallcmd = (
            f"blastall -p blastn -o {self.outdir / 'NC_002696_vs_NC_010338.blast_tab'} "
            f"-i {self.seqdir / 'NC_002696.fna'} "
            f"-d {self.seqdir / 'NC_010338.fna'} "
            "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
        )
        self.blastntgt = [
            (
                f"blastn -out {self.outdir / 'NC_002696_vs_NC_010338.blast_tab'} "
                f"-query {self.seqdir / 'NC_002696.fna'} "
                f"-db {self.seqdir / 'NC_010338.fna'} "
                "-xdrop_gap_final 150 -dust no -evalue 1e-15 "
                "-max_target_seqs 1 -outfmt '6 qseqid sseqid "
                "length mismatch pident nident qlen slen qstart "
                "qend sstart send positive ppos gaps' -task blastn"
            ),
            (
                f"blastn -out {self.outdir / 'NC_010338_vs_NC_002696.blast_tab'} "
                f"-query {self.seqdir / 'NC_010338.fna'} "
                f"-db {self.seqdir / 'NC_002696.fna'} "
                "-xdrop_gap_final 150 -dust no -evalue 1e-15 "
                "-max_target_seqs 1 -outfmt '6 qseqid sseqid length "
                "mismatch pident nident qlen slen qstart qend "
                "sstart send positive ppos gaps' -task blastn"
            ),
        ]
        self.blastalltgt = [
            (
                f"blastall -p blastn -o {self.outdir / 'NC_002696_vs_NC_010338.blast_tab'} "
                f"-i {self.seqdir / 'NC_002696.fna'} "
                f"-d {self.seqdir / 'NC_010338.fna'} "
                "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
            ),
            (
                f"blastall -p blastn -o {self.outdir / 'NC_010338_vs_NC_002696.blast_tab'} "
                f"-i {self.seqdir / 'NC_010338.fna'} "
                f"-d {self.seqdir / 'NC_002696.fna'} "
                "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
            ),
        ]
        self.blastnjobdict = sorted(
            [
                (
                    self.outdir / "NC_002696.fna",
                    f"makeblastdb -dbtype nucl -in {self.seqdir / 'NC_002696.fna'} -title NC_002696 -out {self.outdir / 'NC_002696.fna'}",
                ),
                (
                    self.outdir / "NC_010338.fna",
                    f"makeblastdb -dbtype nucl -in {self.seqdir / 'NC_010338.fna'} -title NC_010338 -out {self.outdir / 'NC_010338.fna'}",
                ),
                (
                    self.outdir / "NC_011916.fna",
                    f"makeblastdb -dbtype nucl -in {self.seqdir / 'NC_011916.fna'} -title NC_011916 -out {self.outdir / 'NC_011916.fna'}",
                ),
                (
                    self.outdir / "NC_014100.fna",
                    f"makeblastdb -dbtype nucl -in {self.seqdir / 'NC_014100.fna'} -title NC_014100 -out {self.outdir / 'NC_014100.fna'}",
                ),
            ]
        )
        self.blastalljobdict = sorted(
            [
                (
                    self.outdir / "NC_002696.fna",
                    f"formatdb -p F -i {self.outdir / 'NC_002696.fna'} -t NC_002696",
                ),
                (
                    self.outdir / "NC_010338.fna",
                    f"formatdb -p F -i {self.outdir / 'NC_010338.fna'} -t NC_010338",
                ),
                (
                    self.outdir / "NC_011916.fna",
                    f"formatdb -p F -i {self.outdir / 'NC_011916.fna'} -t NC_011916",
                ),
                (
                    self.outdir / "NC_014100.fna",
                    f"formatdb -p F -i {self.outdir / 'NC_014100.fna '}-t NC_014100",
                ),
            ]
        )
        self.outdir.mkdir(exist_ok=True)
        self.fmtdboutdir.mkdir(exist_ok=True)
        self.makeblastdbdir.mkdir(exist_ok=True)

    def test_formatdb_generation(self):
        """Generate formatdb command-line."""
        cmd = anib.construct_formatdb_cmd(
            self.seqdir / "NC_002696.fna", self.fmtdboutdir
        )
        self.assertEqual(cmd[0], self.fmtdbcmd)  # correct command
        if not self.fmtdboutdir.is_dir():  # creates new file
            raise AssertionError()

    def test_makeblastdb_generation(self):
        """Generate makeblastdb command-line."""
        cmd = anib.construct_makeblastdb_cmd(
            self.seqdir / "NC_002696.fna", self.makeblastdbdir
        )
        self.assertEqual(cmd[0], self.makeblastdbcmd)  # correct command

    def test_blastdb_commands(self):
        """Generate BLAST+ db commands."""
        # BLAST+
        cmds = anib.generate_blastdb_commands(
            self.blastdbfnames, self.outdir, mode="ANIb"
        )
        self.assertEqual(cmds, self.blastdbtgt)

    def test_legacy_blastdb_commands(self):
        """Generate legacy BLAST db creation commands."""
        cmds = anib.generate_blastdb_commands(
            self.blastdbfnames, self.outdir, mode="ANIblastall"
        )
        self.assertEqual(cmds, self.blastdbtgtlegacy)

    def test_blastn_generation(self):
        """Generate BLASTN+ command-line."""
        cmd = anib.construct_blastn_cmdline(
            self.blastdbfnames[0], self.blastdbfnames[1], self.outdir
        )
        self.assertEqual(cmd, self.blastncmd)

    def test_blastall_generation(self):
        """Generate legacy BLASTN command-line."""
        cmd = anib.construct_blastall_cmdline(
            self.blastdbfnames[0], self.blastdbfnames[1], self.outdir
        )
        self.assertEqual(cmd, self.blastallcmd)

    def test_blastn_commands(self):
        """Generate BLASTN+ commands."""
        # BLAST+
        cmds = anib.generate_blastn_commands(
            self.blastdbfnames, self.outdir, mode="ANIb"
        )
        print(cmds)
        print(self.blastntgt)
        self.assertEqual(cmds, self.blastntgt)

    def test_legacy_blastn_commands(self):
        """Generate legacy BLASTN commands."""
        cmds = anib.generate_blastn_commands(
            self.blastdbfnames, self.outdir, mode="ANIblastall"
        )
        self.assertEqual(cmds, self.blastalltgt)

    def test_blastall_dbjobdict(self):
        """Generate dictionary of legacy BLASTN database jobs."""
        blastcmds = anib.make_blastcmd_builder("ANIblastall", self.outdir)
        jobdict = anib.build_db_jobs(self.infiles, blastcmds)
        self.assertEqual(
            sorted([(k, v.script) for (k, v) in jobdict.items()]), self.blastalljobdict
        )

    def test_blastn_dbjobdict(self):
        """Generate dictionary of BLASTN+ database jobs."""
        blastcmds = anib.make_blastcmd_builder("ANIb", self.outdir)
        jobdict = anib.build_db_jobs(self.infiles, blastcmds)
        print(sorted([(k, v.script) for (k, v) in jobdict.items()]))
        print(self.blastnjobdict)
        self.assertEqual(
            sorted([(k, v.script) for (k, v) in jobdict.items()]), self.blastnjobdict
        )

    def test_blastn_graph(self):
        """Create jobgraph for BLASTN+ jobs."""
        fragresult = anib.fragment_fasta_files(self.infiles, self.outdir, self.fraglen)
        blastcmds = anib.make_blastcmd_builder("ANIb", self.outdir)
        jobgraph = anib.make_job_graph(self.infiles, fragresult[0], blastcmds)
        # We check that the main script job is a blastn job, and that there
        # is a single dependency, which is a makeblastdb job
        for job in jobgraph:
            if not job.script.startswith("blastn"):
                raise AssertionError()
            self.assertEqual(1, len(job.dependencies))
            dep = job.dependencies[0]
            if not dep.script.startswith("makeblastdb"):
                raise AssertionError()

    def test_blastall_graph(self):
        """Create jobgraph for legacy BLASTN jobs."""
        fragresult = anib.fragment_fasta_files(self.infiles, self.outdir, self.fraglen)
        blastcmds = anib.make_blastcmd_builder("ANIblastall", self.outdir)
        jobgraph = anib.make_job_graph(self.infiles, fragresult[0], blastcmds)
        # We check that the main script job is a blastn job, and that there
        # is a single dependency, which is a makeblastdb job
        for job in jobgraph:
            if not job.script.startswith("blastall -p blastn"):
                raise AssertionError()
            self.assertEqual(1, len(job.dependencies))
            dep = job.dependencies[0]
            if not dep.script.startswith("formatdb"):
                raise AssertionError()


class TestFragments(unittest.TestCase):

    """Class defining tests of ANIb FASTA fragmentation."""

    def setUp(self):
        """Initialise parameters for tests."""
        testdir = Path("tests")
        self.outdir = testdir / "test_output" / "anib"
        self.seqdir = testdir / "test_input" / "sequences"
        self.infnames = [
            self.seqdir / fname
            for fname in (
                "NC_002696.fna",
                "NC_010338.fna",
                "NC_011916.fna",
                "NC_014100.fna",
            )
        ]
        self.outfnames = [
            self.outdir / fname
            for fname in (
                "NC_002696-fragments.fna",
                "NC_010338-fragments.fna",
                "NC_011916-fragments.fna",
                "NC_014100-fragments.fna",
            )
        ]
        self.fraglen = 1000
        self.outdir.mkdir(exist_ok=True)

    def test_fragment_files(self):
        """Fragment files for ANIb/ANIblastall."""
        result = anib.fragment_fasta_files(self.infnames, self.outdir, self.fraglen)
        # Are files created?
        for outfname in self.outfnames:
            if not outfname.is_file():
                raise AssertionError()

        # Test fragment lengths
        for _, fragdict in result[-1].items():
            for _, fraglen in fragdict.items():
                if not fraglen <= self.fraglen:
                    raise AssertionError()


class TestParsing(unittest.TestCase):

    """Class defining tests of BLAST output parsing."""

    def setUp(self):
        """Set up test parameters."""
        testdir = Path("tests")
        self.indir = testdir / "test_input" / "anib"
        self.seqdir = testdir / "test_input" / "sequences"
        self.fragdir = testdir / "test_input/anib" / "fragfiles"
        self.anibdir = testdir / "test_input/anib" / "blastn"
        self.aniblastalldir = testdir / "test_input" / "anib" / "blastall"
        self.fname_legacy = self.indir / "NC_002696_vs_NC_010338.blast_tab"
        self.fname = self.indir / "NC_002696_vs_NC_011916.blast_tab"
        self.fragfname = self.indir / "NC_002696-fragments.fna"
        self.fraglens = 1000
        self.infnames = [
            self.seqdir / fname
            for fname in (
                "NC_002696.fna",
                "NC_010338.fna",
                "NC_011916.fna",
                "NC_014100.fna",
            )
        ]
        self.fragfiles = [
            self.fragdir / fname
            for fname in (
                "NC_002696-fragments.fna",
                "NC_010338-fragments.fna",
                "NC_011916-fragments.fna",
                "NC_014100-fragments.fna",
            )
        ]
        self.anibtgt = pd.DataFrame(
            [
                [1.000_000, 0.796_974, 0.999_977, 0.837_285],
                [0.795_958, 1.000_000, 0.795_917, 0.798_250],
                [0.999_922, 0.795_392, 1.000_000, 0.837_633],
                [0.836_780, 0.798_704, 0.836_823, 1.000_000],
            ],
            columns=["NC_002696", "NC_010338", "NC_011916", "NC_014100"],
            index=["NC_002696", "NC_010338", "NC_011916", "NC_014100"],
        )
        self.aniblastalltgt = pd.DataFrame(
            [
                [1.000_000, 0.785_790, 0.999_977, 0.830_641],
                [0.781_319, 1.000_000, 0.781_281, 0.782_723],
                [0.999_937, 0.782_968, 1.000_000, 0.830_431],
                [0.828_919, 0.784_533, 0.828_853, 1.000_000],
            ],
            columns=["NC_002696", "NC_010338", "NC_011916", "NC_014100"],
            index=["NC_002696", "NC_010338", "NC_011916", "NC_014100"],
        )

    def test_parse_blasttab(self):
        """Parses ANIb .blast_tab output."""
        fragdata = anib.get_fraglength_dict([self.fragfname])
        result = anib.parse_blast_tab(self.fname, fragdata, mode="ANIb")
        self.assertEqual(result, (4_016_551, 93, 99.997_693_577_050_029))

    def test_parse_legacy_blasttab(self):
        """Parses legacy .blast_tab output."""
        # ANIblastall output
        fragdata = anib.get_fraglength_dict([self.fragfname])
        result = anib.parse_blast_tab(self.fname_legacy, fragdata, mode="ANIblastall")
        self.assertEqual(result, (1_966_922, 406_104, 78.578_978_313_253_018))

    def test_blastdir_processing(self):
        """Parses directory of .blast_tab output."""
        orglengths = pyani_files.get_sequence_lengths(self.infnames)
        fraglengths = anib.get_fraglength_dict(self.fragfiles)
        # ANIb
        result = anib.process_blast(self.anibdir, orglengths, fraglengths, mode="ANIb")
        assert_frame_equal(
            result.percentage_identity.sort_index(1).sort_index(),
            self.anibtgt.sort_index(1).sort_index(),
        )

    def test_legacy_blastdir_processing(self):
        """Parses directory of legacy .blast_tab output."""
        orglengths = pyani_files.get_sequence_lengths(self.infnames)
        fraglengths = anib.get_fraglength_dict(self.fragfiles)
        result = anib.process_blast(
            self.aniblastalldir, orglengths, fraglengths, mode="ANIblastall"
        )
        assert_frame_equal(
            result.percentage_identity.sort_index(1).sort_index(),
            self.aniblastalltgt.sort_index(1).sort_index(),
        )
