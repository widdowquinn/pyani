#!/usr/bin/env python

"""Tests for pyani package intermediate file parsing

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
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


# Run as script
if __name__ == '__main__':
    import inspect
    import test_parsing
    functions = [o[0] for o in inspect.getmembers(test_parsing) if
                 inspect.isfunction(o[1])]
    for fn in functions:
        print("\nFunction called: {}()".format(fn))
        locals()[fn]()
