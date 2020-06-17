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
# Copyright (c) 2017-2019 The James Hutton Institute
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

import unittest

from pathlib import Path

import pandas as pd
import pytest

from pandas.util.testing import assert_frame_equal

from pyani import anib, pyani_files, pyani_tools


# Test legacy BLAST (blastall) command generation
def test_blastall_dbjobdict(path_fna_all, tmp_path):
    """Generate dictionary of legacy BLASTN database jobs."""
    blastcmds = anib.make_blastcmd_builder("ANIblastall", tmp_path)
    jobdict = anib.build_db_jobs(path_fna_all, blastcmds)
    expected = [
        (tmp_path / _.name, f"formatdb -p F -i {tmp_path / _.name} -t {_.stem}")
        for _ in path_fna_all
    ]
    assert sorted([(k, v.script) for (k, v) in jobdict.items()]) == sorted(expected)


def test_blastall_graph(path_fna_all, tmp_path, fragment_length):
    """Create jobgraph for legacy BLASTN jobs."""
    fragresult = anib.fragment_fasta_files(path_fna_all, tmp_path, fragment_length)
    blastcmds = anib.make_blastcmd_builder("ANIblastall", tmp_path)
    jobgraph = anib.make_job_graph(path_fna_all, fragresult[0], blastcmds)
    # We check that the main script job is a blastn job, and that there
    # is a single dependency, which is a makeblastdb job
    for job in jobgraph:
        assert job.script.startswith("blastall -p blastn")
        assert len(job.dependencies) == 1
        assert job.dependencies[0].script.startswith("formatdb")


def test_blastall_multiple(path_fna_two, tmp_path):
    """Generate legacy BLASTN commands."""
    cmds = anib.generate_blastn_commands(path_fna_two, tmp_path, mode="ANIblastall")
    expected = [
        (
            f"blastall -p blastn -o {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
            f"-i {path_fna_two[0]} "
            f"-d {path_fna_two[1]} "
            "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
        ),
        (
            f"blastall -p blastn -o {tmp_path / str(path_fna_two[1].stem + '_vs_' + path_fna_two[0].stem + '.blast_tab')} "
            f"-i {path_fna_two[1]} "
            f"-d {path_fna_two[0]} "
            "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
        ),
    ]
    assert cmds == expected


def test_blastall_single(path_fna_two, tmp_path):
    """Generate legacy BLASTN command-line."""
    cmd = anib.construct_blastall_cmdline(path_fna_two[0], path_fna_two[1], tmp_path)
    expected = (
        f"blastall -p blastn -o {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
        f"-i {path_fna_two[0]} "
        f"-d {path_fna_two[1]} "
        "-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8"
    )
    assert cmd == expected


# Test BLAST+ (blastn) command generation
def test_blastn_dbjobdict(path_fna_all, tmp_path):
    """Generate dictionary of BLASTN+ database jobs."""
    blastcmds = anib.make_blastcmd_builder("ANIb", tmp_path)
    jobdict = anib.build_db_jobs(path_fna_all, blastcmds)
    expected = [
        (
            tmp_path / _.name,
            f"makeblastdb -dbtype nucl -in {_} -title {_.stem} -out {tmp_path / _.name}",
        )
        for _ in path_fna_all
    ]
    assert sorted([(k, v.script) for (k, v) in jobdict.items()]) == sorted(expected)


def test_blastn_graph(path_fna_all, tmp_path, fragment_length):
    """Create jobgraph for BLASTN+ jobs."""
    fragresult = anib.fragment_fasta_files(path_fna_all, tmp_path, fragment_length)
    blastcmds = anib.make_blastcmd_builder("ANIb", tmp_path)
    jobgraph = anib.make_job_graph(path_fna_all, fragresult[0], blastcmds)
    # We check that the main script job is a blastn job, and that there
    # is a single dependency, which is a makeblastdb job
    for job in jobgraph:
        assert job.script.startswith("blastn")
        assert len(job.dependencies) == 1
        assert job.dependencies[0].script.startswith("makeblastdb")


def test_blastn_multiple(path_fna_two, tmp_path):
    """Generate BLASTN+ commands."""
    # BLAST+
    cmds = anib.generate_blastn_commands(path_fna_two, tmp_path, mode="ANIb")
    expected = [
        (
            f"blastn -out {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
            f"-query {path_fna_two[0]} "
            f"-db {path_fna_two[1]} "
            "-xdrop_gap_final 150 -dust no -evalue 1e-15 "
            "-max_target_seqs 1 -outfmt '6 qseqid sseqid "
            "length mismatch pident nident qlen slen qstart "
            "qend sstart send positive ppos gaps' -task blastn"
        ),
        (
            f"blastn -out {tmp_path / str(path_fna_two[1].stem + '_vs_' + path_fna_two[0].stem + '.blast_tab')} "
            f"-query {path_fna_two[1]} "
            f"-db {path_fna_two[0]} "
            "-xdrop_gap_final 150 -dust no -evalue 1e-15 "
            "-max_target_seqs 1 -outfmt '6 qseqid sseqid length "
            "mismatch pident nident qlen slen qstart qend "
            "sstart send positive ppos gaps' -task blastn"
        ),
    ]
    assert cmds == expected


def test_blastn_single(path_fna_two, tmp_path):
    """Generate BLASTN+ command-line."""
    cmd = anib.construct_blastn_cmdline(path_fna_two[0], path_fna_two[1], tmp_path)
    expected = (
        f"blastn -out {tmp_path / str(path_fna_two[0].stem + '_vs_' + path_fna_two[1].stem + '.blast_tab')} "
        f"-query {path_fna_two[0]} "
        f"-db {path_fna_two[1]} "
        "-xdrop_gap_final 150 -dust no -evalue 1e-15 -max_target_seqs 1 "
        "-outfmt '6 qseqid sseqid length mismatch pident nident qlen slen "
        "qstart qend sstart send positive ppos gaps' "
        "-task blastn"
    )
    assert cmd == expected


# Test legacy BLAST database formatting (formatdb) command generation
def test_formatdb_multiple(path_fna_two, tmp_path):
    """Generate legacy BLAST db creation commands."""
    cmds = anib.generate_blastdb_commands(path_fna_two, tmp_path, mode="ANIblastall")
    expected = [
        (
            f"formatdb -p F -i {tmp_path / path_fna_two[0].name} -t {path_fna_two[0].stem}",
            tmp_path / path_fna_two[0].name,
        ),
        (
            f"formatdb -p F -i {tmp_path / path_fna_two[1].name} -t {path_fna_two[1].stem}",
            tmp_path / path_fna_two[1].name,
        ),
    ]
    assert cmds == expected


def test_formatdb_single(path_fna, tmp_path):
    """Generate legacy BLAST formatdb command-line."""
    cmd = anib.construct_formatdb_cmd(path_fna, tmp_path)
    expected = f"formatdb -p F -i {tmp_path / path_fna.name} -t {path_fna.stem}"
    assert cmd[0] == expected


# Test BLAST+ database formatting (makeblastdb) command generation
def test_makeblastdb_multiple(path_fna_two, tmp_path):
    """Generate multiple BLAST+ makeblastdb command-lines."""
    cmds = anib.generate_blastdb_commands(path_fna_two, tmp_path, mode="ANIb")
    expected = [
        (
            (
                f"makeblastdb -dbtype nucl -in {path_fna_two[0]} "
                f"-title {path_fna_two[0].stem} -out {tmp_path / path_fna_two[0].name}"
            ),
            tmp_path / path_fna_two[0].name,
        ),
        (
            (
                f"makeblastdb -dbtype nucl -in {path_fna_two[1]} "
                f"-title {path_fna_two[1].stem} -out {tmp_path / path_fna_two[1].name}"
            ),
            tmp_path / path_fna_two[1].name,
        ),
    ]
    assert cmds == expected


def test_makeblastdb_single(path_fna, tmp_path):
    """Generate single BLAST+ makeblastdb command-line."""
    cmd = anib.construct_makeblastdb_cmd(path_fna, tmp_path)
    expected = (
        f"makeblastdb -dbtype nucl -in {path_fna} "
        f"-title {path_fna.stem} -out {tmp_path / path_fna.name}"
    )
    assert cmd[0] == expected


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
