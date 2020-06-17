#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2013-2019
# (c) The University of Strathclude 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2013-2019 The James Hutton Institute
# (c) The University of Strathclude 2019-2020
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
"""Pytest configuration file."""

from pathlib import Path
from typing import List, NamedTuple, Tuple

import pandas as pd
import pytest

# Path to tests, contains tests and data subdirectories
TESTSPATH = Path(__file__).parents[0]
FIXTUREPATH = TESTSPATH / "fixtures"


class ANIbOutput(NamedTuple):

    """Convenience struct for ANIb output."""

    fragfile: Path
    tabfile: Path
    legacytabfile: Path


class ANIbOutputDir(NamedTuple):

    """Convenience struct for ANIb output."""

    infiles: List[Path]
    fragfiles: List[Path]
    blastdir: Path
    legacyblastdir: Path
    blastresult: pd.DataFrame
    legacyblastresult: pd.DataFrame


class DeltaDir(NamedTuple):

    """Convenience struct for MUMmer .delta file and associated parsed output."""

    seqdir: Path
    deltadir: Path
    deltaresult: pd.DataFrame


class DeltaParsed(NamedTuple):

    """Convenience struct for MUMmer .delta file and associated parsed output."""

    filename: Path
    data: Tuple[int]


class MUMmerExample(NamedTuple):

    """Convenience struct for MUMmer command-line examples."""

    infiles: List[Path]
    ncmds: List[str]
    fcmds: List[str]


@pytest.fixture
def anib_output(dir_anib_in):
    """Namedtuple of example ANIb output.

    fragfile - fragmented FASTA query file
    tabfile  - BLAST+ tabular output
    legacytabfile - blastall tabular output
    """
    return ANIbOutput(
        dir_anib_in / "NC_002696-fragments.fna",
        dir_anib_in / "NC_002696_vs_NC_011916.blast_tab",
        dir_anib_in / "NC_002696_vs_NC_010338.blast_tab",
    )


@pytest.fixture
def anib_output_dir(dir_anib_in):
    """Namedtuple of example ANIb output - full directory.

    infiles - list of FASTA query files
    fragfiles - list of fragmented FASTA query files
    blastdir - path to BLAST+ output data
    legacyblastdir - path to blastall output data
    blastresult - pd.DataFrame result for BLAST+
    legacyblastresult - pd.DataFrame result for blastall
    """
    return ANIbOutputDir(
        [
            _
            for _ in (dir_anib_in / "sequences").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        [
            _
            for _ in (dir_anib_in / "fragfiles").iterdir()
            if _.is_file() and _.suffix == ".fna"
        ],
        dir_anib_in / "blastn",
        dir_anib_in / "blastall",
        pd.read_csv(dir_anib_in / "dataframes" / "blastn_result.csv", index_col=0),
        pd.read_csv(dir_anib_in / "dataframes" / "blastall_result.csv", index_col=0),
    )


@pytest.fixture
def delta_output_dir(dir_anim_in):
    """Namedtuple of example MUMmer .delta file output."""
    return DeltaDir(
        dir_anim_in / "sequences",
        dir_anim_in / "deltadir",
        pd.read_csv(dir_anim_in / "dataframes" / "deltadir_result.csv", index_col=0),
    )


@pytest.fixture
def deltafile_parsed(dir_anim_in):
    """Example parsed deltafile data."""
    return DeltaParsed(dir_anim_in / "test.delta", (4074001, 2191))


@pytest.fixture
def dir_anib_in():
    """Input files for ANIb tests."""
    return FIXTUREPATH / "anib"


@pytest.fixture
def dir_anim_in():
    """Input files for ANIm tests."""
    return FIXTUREPATH / "anim"


@pytest.fixture
def dir_seq():
    """Sequence files for tests."""
    return FIXTUREPATH / "sequences"


@pytest.fixture
def dir_targets():
    """Target files for output comparisons."""
    return FIXTUREPATH / "targets"


@pytest.fixture
def dir_tgt_fragments(dir_targets):
    """Target files for FASTA file fragmentation."""
    return dir_targets / "fragments"


@pytest.fixture
def fragment_length():
    """Fragment size for ANIb-related analyses."""
    return 1000


@pytest.fixture
def mummer_cmds_four(path_file_four):
    """Example MUMmer commands (four files)."""
    return MUMmerExample(
        path_file_four,
        [
            "nucmer --mum -p nucmer_output/file1_vs_file2 file1.fna file2.fna",
            "nucmer --mum -p nucmer_output/file1_vs_file3 file1.fna file3.fna",
            "nucmer --mum -p nucmer_output/file1_vs_file4 file1.fna file4.fna",
            "nucmer --mum -p nucmer_output/file2_vs_file3 file2.fna file3.fna",
            "nucmer --mum -p nucmer_output/file2_vs_file4 file2.fna file4.fna",
            "nucmer --mum -p nucmer_output/file3_vs_file4 file3.fna file4.fna",
        ],
        [
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file1_vs_file2.delta "
                "nucmer_output/file1_vs_file2.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file1_vs_file3.delta "
                "nucmer_output/file1_vs_file3.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file1_vs_file4.delta "
                "nucmer_output/file1_vs_file4.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file2_vs_file3.delta "
                "nucmer_output/file2_vs_file3.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file2_vs_file4.delta "
                "nucmer_output/file2_vs_file4.filter"
            ),
            (
                "delta_filter_wrapper.py delta-filter -1 "
                "nucmer_output/file3_vs_file4.delta "
                "nucmer_output/file3_vs_file4.filter"
            ),
        ],
    )


@pytest.fixture
def path_file_two():
    """Path to two arbitrary filenames."""
    return [Path(f"file{_:d}.fna") for _ in range(1, 3)]


@pytest.fixture
def path_file_four():
    """Path to four arbitrary filenames."""
    return [Path(f"file{_:d}.fna") for _ in range(1, 5)]


@pytest.fixture
def path_fna(dir_seq):
    """Path to one .fna sequence file from dir_seq."""
    fnapaths = [_ for _ in dir_seq.iterdir() if _.is_file() and _.suffix == ".fna"]
    return fnapaths[0]


@pytest.fixture
def path_fna_two(dir_seq):
    """Paths to two .fna sequence file in dir_seq."""
    fnapaths = [_ for _ in dir_seq.iterdir() if _.is_file() and _.suffix == ".fna"]
    return fnapaths[:2]


@pytest.fixture
def path_fna_all(dir_seq):
    """Paths to all .fna sequence file in dir_seq."""
    return [_ for _ in dir_seq.iterdir() if _.is_file() and _.suffix == ".fna"]
