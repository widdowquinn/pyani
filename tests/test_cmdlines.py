#!/usr/bin/env python

"""Tests for pyani package command-line generation

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/), from the repository root
directory.
"""


from nose.tools import assert_equal
from pyani import anim


# Test ANIm command-lines
# One pairwise comparison
def test_anim_pairwise_basic():
    """Test generation of basic NUCmer pairwise comparison command."""
    ncmd, fcmd = anim.construct_nucmer_cmdline("file1.fna", "file2.fna")
    print("\n{}\n{}\n".format(ncmd, fcmd))
    assert_equal(ncmd, "nucmer --mum -p ./nucmer_output/file1_vs_file2 " +
                 "file1.fna file2.fna")
    assert_equal(fcmd, "delta_filter_wrapper.py delta-filter -1 " +
                 "./nucmer_output/file1_vs_file2.delta " +
                 "./nucmer_output/file1_vs_file2.filter")


def test_anim_pairwise_maxmatch():
    """Test generation of NUCmer pairwise comparison command with maxmatch."""
    ncmd, fcmd = anim.construct_nucmer_cmdline("file1.fna", "file2.fna",
                                               maxmatch=True)
    print("\n{}\n{}\n".format(ncmd, fcmd))
    assert_equal(ncmd, "nucmer --maxmatch -p ./nucmer_output/file1_vs_file2 " +
                 "file1.fna file2.fna")
    assert_equal(fcmd, "delta_filter_wrapper.py delta-filter -1 " +
                 "./nucmer_output/file1_vs_file2.delta " + 
                 "./nucmer_output/file1_vs_file2.filter")    
    print("\n{}\n{}\n".format(ncmd, fcmd))


# List of pairwise comparisons
def test_anim_collection():
    """Test generation of list of NUCmer comparison commands.
    """
    files = ["file1", "file2", "file3", "file4"]
    ncmds, fcmds = anim.generate_nucmer_commands(files)
    print('\n'.join([ncmd for ncmd in ncmds]))
    print('\n'.join([fcmd for fcmd in fcmds]))
    assert_equal(ncmds, ['nucmer --mum -p ./nucmer_output/file1_vs_file2 file1 file2',
                         'nucmer --mum -p ./nucmer_output/file1_vs_file3 file1 file3',
                         'nucmer --mum -p ./nucmer_output/file1_vs_file4 file1 file4',
                         'nucmer --mum -p ./nucmer_output/file2_vs_file3 file2 file3',
                         'nucmer --mum -p ./nucmer_output/file2_vs_file4 file2 file4',
                         'nucmer --mum -p ./nucmer_output/file3_vs_file4 file3 file4'])
    assert_equal(fcmds, ['delta_filter_wrapper.py delta-filter -1 ' +
                         './nucmer_output/file1_vs_file2.delta ' +
                         './nucmer_output/file1_vs_file2.filter',
                         'delta_filter_wrapper.py delta-filter -1 ' +
                         './nucmer_output/file1_vs_file3.delta ' +
                         './nucmer_output/file1_vs_file3.filter',
                         'delta_filter_wrapper.py delta-filter -1 ' +
                         './nucmer_output/file1_vs_file4.delta ' +
                         './nucmer_output/file1_vs_file4.filter',
                         'delta_filter_wrapper.py delta-filter -1 ' +
                         './nucmer_output/file2_vs_file3.delta ' +
                         './nucmer_output/file2_vs_file3.filter',
                         'delta_filter_wrapper.py delta-filter -1 ' +
                         './nucmer_output/file2_vs_file4.delta ' +
                         './nucmer_output/file2_vs_file4.filter',
                         'delta_filter_wrapper.py delta-filter -1 ' +
                         './nucmer_output/file3_vs_file4.delta ' +
                         './nucmer_output/file3_vs_file4.filter'])
