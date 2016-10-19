#!/usr/bin/env python

"""Tests for availability of pyani dependencies

We only test for dependencies from non-standard libraries.

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
"""

import subprocess
import sys

def test_import_biopython():
    """Test Biopython import."""
    import Bio


def test_import_matplotlib():
    """Test matplotlib import."""
    import matplotlib


def test_import_numpy():
    """Test numpy import."""
    import numpy


def test_import_pandas():
    """Test pandas import."""
    import pandas


def test_import_scipy():
    """Test scipy import."""
    import scipy


def test_run_blast():
    """Test that BLAST+ is runnable."""
    cmd = "blastn -version"
    subprocess.run(cmd, shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE)


# Run as script
if __name__ == '__main__':
    import inspect
    import test_cmdlines
    functions = [o[0] for o in inspect.getmembers(test_cmdlines) if
                 inspect.isfunction(o[1])]
    for fn in functions:
        print("\nFunction called: {}()".format(fn))
        locals()[fn]()
