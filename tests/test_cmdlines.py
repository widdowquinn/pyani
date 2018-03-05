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
    cmd_nucmer, cmd_filter = anim.construct_nucmer_cmdline(
        "file1.fna", "file2.fna")
    tgt_nucmer = ' '.join(["nucmer --mum -p ./nucmer_output/file1_vs_file2",
                           "file1.fna file2.fna"])
    tgt_filter = ' '.join(["delta_filter_wrapper.py delta-filter -1",
                           "./nucmer_output/file1_vs_file2.delta",
                           "./nucmer_output/file1_vs_file2.filter"])
    assert_equal(cmd_nucmer, tgt_nucmer)
    assert_equal(cmd_filter, tgt_filter)


def test_anim_pairwise_maxmatch():
    """Test generation of NUCmer pairwise comparison command with maxmatch."""
    cmd_nucmer, cmd_filter = anim.construct_nucmer_cmdline("file1.fna", "file2.fna",
                                                           maxmatch=True)
    tgt_nucmer = ' '.join(["nucmer --maxmatch -p ./nucmer_output/file1_vs_file2",
                           "file1.fna file2.fna"])
    tgt_filter = ' '.join(["delta_filter_wrapper.py delta-filter -1",
                           "./nucmer_output/file1_vs_file2.delta",
                           "./nucmer_output/file1_vs_file2.filter"])
    assert_equal(cmd_nucmer, tgt_nucmer)
    assert_equal(cmd_filter, tgt_filter)


# List of pairwise comparisons
def test_anim_collection():
    """Test generation of list of NUCmer comparison commands."""
    files = ["file1", "file2", "file3", "file4"]
    cmds_nucmer, cmds_filter = anim.generate_nucmer_commands(files)
    tgts_nucmer = ['nucmer --mum -p ./nucmer_output/file1_vs_file2 file1 file2',
                   'nucmer --mum -p ./nucmer_output/file1_vs_file3 file1 file3',
                   'nucmer --mum -p ./nucmer_output/file1_vs_file4 file1 file4',
                   'nucmer --mum -p ./nucmer_output/file2_vs_file3 file2 file3',
                   'nucmer --mum -p ./nucmer_output/file2_vs_file4 file2 file4',
                   'nucmer --mum -p ./nucmer_output/file3_vs_file4 file3 file4']
    tgts_filter = ['delta_filter_wrapper.py delta-filter -1 ' +
                   './nucmer_output/file1_vs_file2.delta ' +
                   './nucmer_output/file1_vs_file2.filter',
                   'delta_filter_wrapper.py delta-filter -1 ' +
                   './nucmer_output/file1_vs_file3.delta ' +
                   './nucmer_output/file1_vs_file3.filter',
                   'delta_filter_wrapper.py delta-filter -1 ' +
                   './nucmer_output/file1_vs_file4.delta ' +
                   './nucmer_output/file1_vs_file4.filter'
                   'delta_filter_wrapper.py delta-filter -1 ' +
                   './nucmer_output/file2_vs_file3.delta ' +
                   './nucmer_output/file2_vs_file3.filter',
                   'delta_filter_wrapper.py delta-filter -1 ' +
                   './nucmer_output/file2_vs_file4.delta ' +
                   './nucmer_output/file2_vs_file4.filter',
                   'delta_filter_wrapper.py delta-filter -1 ' +
                   './nucmer_output/file3_vs_file4.delta ' +
                   './nucmer_output/file3_vs_file4.filter']
    assert_equal(cmds_nucmer, tgts_nucmer)
    assert_equal(cmds_filter, cmds_filter)
