#!/usr/bin/env python

"""test_anim.py

Test anim.py module.

These tests are intended to be run from the repository root using:

nosetests -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.

(c) The James Hutton Institute 2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2017 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
import unittest

from nose.tools import assert_equal

from pyani import (anim,)


class TestNUCmerCmdline(unittest.TestCase):

    """Class defining tests of NUCmer command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        # Basic NUCmer and delta-filter command-line targets
        self.ntgt = ' '.join(["nucmer --mum -p",
                              "tests/test_output/nucmer_output/file1_vs_file2",
                              "file1.fna file2.fna"])
        self.ntgtmax = ' '.join(["nucmer --maxmatch -p",
                                 "tests/test_output/nucmer_output/file1_vs_file2",
                                 "file1.fna file2.fna"])
        self.ftgt = ' '.join(["delta_filter_wrapper.py delta-filter -1",
                              "tests/test_output/nucmer_output/file1_vs_file2.delta",
                              "tests/test_output/nucmer_output/file1_vs_file2.filter"])
        self.ncmdlist = ['nucmer --mum -p ./nucmer_output/file1_vs_file2 file1 file2',
                         'nucmer --mum -p ./nucmer_output/file1_vs_file3 file1 file3',
                         'nucmer --mum -p ./nucmer_output/file1_vs_file4 file1 file4',
                         'nucmer --mum -p ./nucmer_output/file2_vs_file3 file2 file3',
                         'nucmer --mum -p ./nucmer_output/file2_vs_file4 file2 file4',
                         'nucmer --mum -p ./nucmer_output/file3_vs_file4 file3 file4']
        self.fcmdlist = [' '.join(['delta_filter_wrapper.py delta-filter -1',
                                   './nucmer_output/file1_vs_file2.delta',
                                   './nucmer_output/file1_vs_file2.filter']),
                         ' '.join(['delta_filter_wrapper.py delta-filter -1',
                                   './nucmer_output/file1_vs_file3.delta',
                                   './nucmer_output/file1_vs_file3.filter']),
                         ' '.join(['delta_filter_wrapper.py delta-filter -1',
                                   './nucmer_output/file1_vs_file4.delta',
                                   './nucmer_output/file1_vs_file4.filter']),
                         ' '.join(['delta_filter_wrapper.py delta-filter -1',
                                   './nucmer_output/file2_vs_file3.delta',
                                   './nucmer_output/file2_vs_file3.filter']),
                         ' '.join(['delta_filter_wrapper.py delta-filter -1',
                                   './nucmer_output/file2_vs_file4.delta',
                                   './nucmer_output/file2_vs_file4.filter']),
                         ' '.join(['delta_filter_wrapper.py delta-filter -1',
                                  './nucmer_output/file3_vs_file4.delta',
                                  './nucmer_output/file3_vs_file4.filter'])]
        self.outdir = os.path.join('tests', 'test_output')
        

    def test_single_cmd_generation(self):
        """generate single abstract NUCmer/delta-filter command-line."""
        cmds = anim.construct_nucmer_cmdline("file1.fna", "file2.fna",
                                             outdir=self.outdir)
        assert_equal(cmds, (self.ntgt, self.ftgt))

    def test_maxmatch_cmd_generation(self):
        """generate NUCmer command line with maxmatch."""
        ncmd, fcmd = anim.construct_nucmer_cmdline("file1.fna", "file2.fna",
                                                   outdir=self.outdir,
                                                   maxmatch=True)
        assert_equal(ncmd, self.ntgtmax)
        
    def test_multi_cmd_generation(self):
        """generate multiple abstract NUCmer/delta-filter command-lines."""
        files = ["file1", "file2", "file3", "file4"]
        cmds = anim.generate_nucmer_commands(files)
        assert_equal(cmds, (self.ncmdlist, self.fcmdlist))
