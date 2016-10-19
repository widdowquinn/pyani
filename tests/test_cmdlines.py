#!/usr/bin/env python

"""Tests for pyani package command-line generation

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/), from the repository root
directory.

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
"""


from nose.tools import assert_equal
from pyani import anim


# Test ANIm command-lines
# One pairwise comparison
def test_anim_pairwise_basic():
    """Test generation of basic NUCmer pairwise comparison command.
    """
    cmd = anim.construct_nucmer_cmdline("file1.fna", "file2.fna")
    assert_equal(cmd, "nucmer -mum -p ./nucmer_output/file1_vs_file2 " +
                 "file1.fna file2.fna")
    print(cmd)

def test_anim_pairwise_maxmatch():
    """Test generation of NUCmer pairwise comparison command with maxmatch.
    """
    cmd = anim.construct_nucmer_cmdline("file1.fna", "file2.fna",
                                        maxmatch=True)
    assert_equal(cmd, "nucmer -maxmatch -p ./nucmer_output/file1_vs_file2 " +
                 "file1.fna file2.fna")
    print(cmd)


# List of pairwise comparisons
def test_anim_collection():
    """Test generation of list of NUCmer comparison commands.
    """
    files = ["file1", "file2", "file3", "file4"]
    cmdlist = anim.generate_nucmer_commands(files)
    assert_equal(cmdlist, ['nucmer -mum -p ./nucmer_output/file1_vs_file2 ' +
                           'file1 file2',
                           'nucmer -mum -p ./nucmer_output/file1_vs_file3 ' +
                           'file1 file3',
                           'nucmer -mum -p ./nucmer_output/file1_vs_file4 ' +
                           'file1 file4',
                           'nucmer -mum -p ./nucmer_output/file2_vs_file3 ' +
                           'file2 file3',
                           'nucmer -mum -p ./nucmer_output/file2_vs_file4 ' +
                           'file2 file4',
                           'nucmer -mum -p ./nucmer_output/file3_vs_file4 ' +
                           'file3 file4'])
    print(cmdlist)


# Run as script
if __name__ == '__main__':
    import inspect
    import test_cmdlines
    functions = [o[0] for o in inspect.getmembers(test_cmdlines) if
                 inspect.isfunction(o[1])]
    for fn in functions:
        print("\nFunction called: {}()".format(fn))
        locals()[fn]()
