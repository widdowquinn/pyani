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

import pytest

# Path to tests, contains tests and data subdirectories
TESTSPATH = Path(__file__).parents[0]
FIXTUREPATH = TESTSPATH / "fixtures"


@pytest.fixture
def dir_anib_in():
    """Input files for ANIb tests."""
    return FIXTUREPATH / "anib"


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
