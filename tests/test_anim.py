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

import pandas as pd

from nose.tools import (assert_equal,)
from pandas.util.testing import (assert_frame_equal,)

from pyani import (anim, pyani_files)


class TestNUCmerCmdline(unittest.TestCase):

    """Class defining tests of NUCmer command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        # Basic NUCmer and delta-filter command-line targets
        self.ntgt = ' '.join(["nucmer --mum -p",
                              "tests/test_output/anim/nucmer_output/file1_vs_file2",
                              "file1.fna file2.fna"])
        self.ntgtmax = ' '.join(["nucmer --maxmatch -p",
                                 "tests/test_output/anim/nucmer_output/file1_vs_file2",
                                 "file1.fna file2.fna"])
        self.ftgt = ' '.join(["delta_filter_wrapper.py delta-filter -1",
                              "tests/test_output/anim/nucmer_output/file1_vs_file2.delta",
                              "tests/test_output/anim/nucmer_output/file1_vs_file2.filter"])
        self.files = ["file1", "file2", "file3", "file4"]
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
        self.seqdir = os.path.join('tests', 'test_input', 'sequences')
        self.outdir = os.path.join('tests', 'test_output', 'anim')
        self.indir = os.path.join('tests', 'test_input', 'anim')
        self.deltafile = os.path.join(self.indir, 'test.delta')
        self.deltadir = os.path.join(self.indir, 'deltadir')
        self.df_pid = pd.DataFrame([[1.000000, 0.850994, 0.999974, 0.867940],
                                    [0.850994, 1.000000, 0.851074, 0.852842],
                                    [0.999974, 0.851074, 1.000000, 0.867991],
                                    [0.867940, 0.852842, 0.867991, 1.000000]],
                                   columns=['NC_002696',  'NC_010338',
                                            'NC_011916',  'NC_014100'],
                                   index=['NC_002696',  'NC_010338',
                                          'NC_011916',  'NC_014100'])

    def test_single_cmd_generation(self):
        """generate single abstract NUCmer/delta-filter command-line.

        Tests that a single NUCmer/delta-filter command-line pair is
        produced correctly
        """
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
        """generate multiple abstract NUCmer/delta-filter command-lines.

        Tests that all the input files are correctly-paired
        """
        cmds = anim.generate_nucmer_commands(self.files)
        assert_equal(cmds, (self.ncmdlist, self.fcmdlist))

    def test_nucmer_job_generation(self):
        """generate dependency tree of NUCmer/delta-filter jobs.

        Tests that the correct dependency graph and naming scheme is produced.
        """
        joblist = anim.generate_nucmer_jobs(self.files,
                                            jobprefix="test")
        assert_equal(len(joblist), 6)
        for idx, job in enumerate(joblist):
            assert_equal(job.name, "test_%06d-f" % idx)  # filter job name
            assert_equal(len(job.dependencies), 1)       # has NUCmer job
            assert_equal(job.dependencies[0].name,
                         "test_%06d-n" % idx)            # NUCmer job name

    def test_deltafile_import(self):
        """parses NUCmer .delta/.filter file."""
        result = anim.parse_delta(self.deltafile)
        assert_equal(result, (4073917, 2191))

    def test_process_deltadir(self):
        """processes directory of .delta files into ANIResults."""
        seqfiles = pyani_files.get_fasta_files(self.seqdir)
        orglengths = pyani_files.get_sequence_lengths(seqfiles)
        result = anim.process_deltadir(self.deltadir, orglengths)
        assert_frame_equal(result.percentage_identity.sort_index(1).sort_index(),
                           self.df_pid.sort_index(1).sort_index())
