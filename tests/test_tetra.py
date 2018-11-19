#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test_tetra.py

Test tetra.py module.

These tests are intended to be run from the repository root using:

pytest -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.

(c) The James Hutton Institute 2017-2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2017-2018 The James Hutton Institute

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

import json
import os
import unittest

import pandas as pd

from pandas.util.testing import assert_frame_equal

from pyani import tetra


def ordered(obj):
    if isinstance(obj, dict):
        return sorted((k, ordered(v)) for k, v in obj.items())
    elif isinstance(obj, list):
        return sorted(ordered(x) for x in obj)
    else:
        return obj


class TestTETRA(unittest.TestCase):

    """Class defining tests of TETRA algorithm."""

    def setUp(self):
        """Define parameters and values for tests."""
        self.indir = os.path.join("tests", "test_input", "tetra")
        self.tgtdir = os.path.join("tests", "test_targets", "tetra")
        self.seqdir = os.path.join("tests", "test_input", "sequences")
        self.infile = os.path.join(self.seqdir, "NC_002696.fna")
        self.infiles = [
            os.path.join(self.seqdir, fname) for fname in os.listdir(self.seqdir)
        ]

    def test_tetraclean(self):
        """detects unambiguous IUPAC symbols correctly."""
        self.assertFalse(tetra.tetra_clean("ACGTYACGTACNGTACGWTACGT"))
        self.assertTrue(tetra.tetra_clean("ACGTACGTACGTACGTACGTAC"))

    def test_zscore(self):
        """TETRA Z-score calculated correctly."""
        tetra_z = tetra.calculate_tetra_zscore(self.infile)
        with open(os.path.join(self.tgtdir, "zscore.json"), "r") as ifh:
            target = json.load(ifh)
        self.assertEqual(ordered(tetra_z), ordered(target))

    def test_correlations(self):
        """TETRA correlation calculated correctly."""
        infiles = ordered(self.infiles)[:2]  # only test a single correlation
        corr = tetra.calculate_correlations(tetra.calculate_tetra_zscores(infiles))
        target = pd.read_csv(
            os.path.join(self.tgtdir, "correlation.tab"), sep="\t", index_col=0
        )
        assert_frame_equal(corr, target)
