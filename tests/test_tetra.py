#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) The University of Strathclude 2019-2022
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
# Copyright (c) 2017-2019 The James Hutton Institute
# (c) The University of Strathclude 2019-2022
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
"""Test tetra.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import json
import unittest

from pathlib import Path

import pandas as pd

from pandas.testing import assert_frame_equal

from pyani.tetra import (
    calculate_correlations,
    calculate_tetra_zscore,
    calculate_tetra_zscores,
    tetra_clean,
)


def ordered(obj):
    """Return consistently ordered version of the passed object."""
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    return obj


def test_tetraclean():
    """Detect unambiguous IUPAC symbols correctly."""
    assert (
        tetra_clean("ACGTYACGTACNGTACGWTACGT"),
        tetra_clean("ACGTACGTACGTACGTACGTAC"),
    ) == (False, True)


def test_zscore(dir_seq, dir_targets):
    """Test that TETRA Z-score calculated correctly."""
    tetra_z = calculate_tetra_zscore(dir_seq / "NC_002696.fna")
    with (dir_targets / "tetra" / "zscore.json").open("r") as ifh:
        target = json.load(ifh)
    assert ordered(tetra_z) == ordered(target)


def test_correlations(path_fna_all, dir_targets):
    """Test that TETRA correlation calculated correctly."""
    infiles = ordered(path_fna_all)[:2]  # only test a single correlation
    corr = calculate_correlations(calculate_tetra_zscores(infiles))
    target = pd.read_csv(
        dir_targets / "tetra" / "correlation.tab", sep="\t", index_col=0
    )
    assert_frame_equal(corr, target)
