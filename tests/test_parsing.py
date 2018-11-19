#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test_parsing.py

Tests for pyani package intermediate file parsing

These tests are intended to be run from the repository root using:

pytest -v

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

import os
import unittest

from pyani import anim

# Work out where we are. We need to do this to find related data files
# for testing
curdir = os.path.dirname(os.path.abspath(__file__))

# Path to test .delta file
DELTAFILE = os.path.join(curdir, "test_ani_data", "NC_002696_vs_NC_011916.delta")


class TestIntermediateParsing(unittest.TestCase):

    """Test parsing of intermediate files"""

    def test_anim_delta(self):
        """Test parsing of NUCmer delta file."""
        aln, sim = anim.parse_delta(DELTAFILE)
        self.assertEqual(aln, 4074001)
        self.assertEqual(sim, 2191)
        print("Alignment length: {0}\nSimilarity Errors: {1}".format(aln, sim))
