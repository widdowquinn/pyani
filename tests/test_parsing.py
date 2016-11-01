#!/usr/bin/env python

"""Tests for pyani package intermediate file parsing

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).
"""

import os

from nose.tools import assert_equal
from pyani import anim

# Work out where we are. We need to do this to find related data files
# for testing
curdir = os.path.dirname(os.path.abspath(__file__))

# Path to test .delta file
DELTAFILE = os.path.join(curdir, 'test_ani_data',
                         'NC_002696_vs_NC_011916.delta')

# Test ANIm command-lines
# One pairwise comparison
def test_anim_delta():
    """Test parsing of NUCmer delta file."""
    aln, sim = anim.parse_delta(DELTAFILE)
    assert_equal(aln, 4073917)
    assert_equal(sim, 2191)
    print("Alignment length: {0}\nSimilarity Errors: {1}".format(aln, sim))
