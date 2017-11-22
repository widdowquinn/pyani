#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""test_multiprocessing.py

Test run_multiprocessing.py module.

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

from nose.tools import (assert_equal, )

from pyani import (run_multiprocessing, pyani_jobs, anib)


class TestMultiprocessing(unittest.TestCase):

    """Class defining tests of pyani's multiprocessing module."""

    def setUp(self):
        """Define parameters and arguments for tests."""
        self.cmdlist = ['for i in %s; do echo "Thread %d: value ${i}"; done' %
                        (' '.join([str(e) for e in range(v)]), v) for
                        v in range(5)]
        self.cmds = ['ls -ltrh', 'echo ${PWD}']
        self.seqdir = os.path.join('tests', 'test_input', 'sequences')
        self.outdir = os.path.join('tests', 'test_output', 'multiprocessing')
        self.infiles = [os.path.join(self.seqdir, fname) for fname in
                        os.listdir(self.seqdir)][:2]
        self.fraglen = 1000
        os.makedirs(self.outdir, exist_ok=True)

    def test_multiprocessing_run(self):
        """multiprocessing() runs basic jobs."""
        result = run_multiprocessing.multiprocessing_run(self.cmdlist)
        assert_equal(0, result)

    def test_cmdsets(self):
        """module builds command sets."""
        job1 = pyani_jobs.Job('dummy_with_dependency', self.cmds[0])
        job2 = pyani_jobs.Job('dummy_dependency', self.cmds[1])
        job1.add_dependency(job2)
        cmdsets = run_multiprocessing.populate_cmdsets(job1, list(), depth=1)
        target = [{cmd} for cmd in self.cmds]
        assert_equal(cmdsets, target)

    def test_dependency_graph_run(self):
        """module runs dependency graph."""
        fragresult = anib.fragment_fasta_files(self.infiles, self.outdir,
                                               self.fraglen)
        blastcmds = anib.make_blastcmd_builder("ANIblastall", self.outdir)
        jobgraph = anib.make_job_graph(self.infiles, fragresult[0],
                                       blastcmds)
        result = run_multiprocessing.run_dependency_graph(jobgraph)
        assert_equal(0, result)
