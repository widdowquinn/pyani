#!/usr/bin/env python

"""Tests for pyani package command-line generation
"""

from pyani import anim

# Test ANIm command-lines
# One pairwise comparison
def test_anim_pairwise():
    """Test generation of pairwise comparison command.
    """
    cmd = anim.construct_nucmer_cmdline("file1", "file2")
    print cmd

# List of pairwise comparisons
def test_anim_pairwise():
    """Test generation of pairwise comparison command.
    """
    files = ["file1", "file2", "file3", "file4"]
    cmdlist = anim.generate_nucmer_commands(files)
    print cmdlist
    
    
