#!/usr/bin/env python

"""Tests for pyani package multiprocessing runs

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
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


# Run as script
if __name__ == '__main__':
    import inspect
    import test_multiprocessing
    functions = [o[0] for o in inspect.getmembers(test_multiprocessing) if
                 inspect.isfunction(o[1])]
    for fn in functions:
        print("\nFunction called: {}()".format(fn))
        locals()[fn]()
