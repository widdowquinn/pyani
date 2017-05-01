#!/usr/bin/env python

"""Tests for pyani package multiprocessing runs

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).
"""


from nose.tools import assert_equal
from pyani import run_multiprocessing


# Test ANIm command-lines
# One pairwise comparison
def test_multiprocessing_run():
    """Test basic multiprocessing function
    """
    cmdlist = ['for i in %s; do echo "Thread %d: value ${i}"; done' %
               (' '.join([str(e) for e in range(v)]), v) for
               v in range(5)]
    run_multiprocessing.multiprocessing_run(cmdlist)
