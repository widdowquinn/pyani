#!/usr/bin/env python

"""Tests for availability of pyani dependencies

We only test for dependencies from non-standard libraries.

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
"""


def test_dependency_biopython():
    """Test Biopython import."""
    import Bio


def test_dependency_matplotlib():
    """Test matplotlib import."""
    import matplotlib


def test_dependency_numpy():
    """Test numpy import."""
    import numpy


def test_dependency_pandas():
    """Test pandas import."""
    import pandas


def test_dependency_scipy():
    """Test scipy import."""
    import scipy


# Run as script
if __name__ == '__main__':
    import inspect
    import test_cmdlines
    functions = [o[0] for o in inspect.getmembers(test_cmdlines) if
                 inspect.isfunction(o[1])]
    for fn in functions:
        print("\nFunction called: {}()".format(fn))
        locals()[fn]()
