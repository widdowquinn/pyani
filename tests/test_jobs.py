#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""test_jobs.py

Test pyani_jobs.py module.

These tests are intended to be run from the repository root using:

nosetests -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.

(c) The James Hutton Institute 2017-2018
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2017-2018 The James Hutton Institute

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

import unittest

from nose.tools import (assert_equal, )

from pyani import (pyani_jobs, )


class TestJob(unittest.TestCase):

    """Class defining tests of Job object generation."""

    def setUp(self):
        """Define parameters and values for tests."""
        self.cmds = ['ls -ltrh', 'echo ${PWD}']

    def test_create_job(self):
        """create a dummy job."""
        job = pyani_jobs.Job('empty', '')
        assert_equal(job.script, "")

    def test_create_job_with_command(self):
        """create dummy job with command."""
        job = pyani_jobs.Job('dummy', self.cmds[0])
        assert_equal(job.script, self.cmds[0])

    def test_add_dependency(self):
        """create dummy job with dependency."""
        job1 = pyani_jobs.Job('dummy_with_dependency', self.cmds[0])
        job2 = pyani_jobs.Job('dummy_dependency', self.cmds[1])
        job1.add_dependency(job2)
        assert_equal(self.cmds[0], job1.script)
        assert_equal(1, len(job1.dependencies))
        dep = job1.dependencies[0]
        assert_equal(self.cmds[1], dep.script)

    def test_remove_dependency(self):
        """create dummy job, add and remove dependency."""
        job1 = pyani_jobs.Job('dummy_with_dependency', self.cmds[0])
        job2 = pyani_jobs.Job('dummy_dependency', self.cmds[1])
        job1.add_dependency(job2)
        assert_equal(1, len(job1.dependencies))
        dep = job1.dependencies[0]
        assert_equal(self.cmds[1], dep.script)
        job1.remove_dependency(dep)
        assert_equal(0, len(job1.dependencies))


class TestJobGroup(unittest.TestCase):

    """Class defining tests of JobGroup object generation."""

    def setUp(self):
        """define parameters and values for tests."""
        self.emptyscript = 'let "TASK_ID=$SGE_TASK_ID - 1"\n\n\n\n'
        self.params1 = {'-f': ['file1', 'file2', 'file3']}
        self.p1script = "".join(['let "TASK_ID=$SGE_TASK_ID - 1"\n',
                                 '-f_ARRAY=( file1 file2 file3  )\n\n',
                                 'let "-f_INDEX=$TASK_ID % 3"\n',
                                 '-f=${-f_ARRAY[$-f_INDEX]}\n',
                                 'let "TASK_ID=$TASK_ID / 3"\n\n',
                                 'cat\n'])
        self.params2 = {'-f': ['file1', 'file2'],
                        '--format': ['fmtA', 'fmtB']}
        self.p2script = "".join(['let "TASK_ID=$SGE_TASK_ID - 1"\n',
                                 '-f_ARRAY=( file1 file2  )\n',
                                 '--format_ARRAY=( fmtA fmtB  )\n\n',
                                 'let "-f_INDEX=$TASK_ID % 2"\n',
                                 '-f=${-f_ARRAY[$-f_INDEX]}\n',
                                 'let "TASK_ID=$TASK_ID / 2"\n',
                                 'let "--format_INDEX=$TASK_ID % 2"\n',
                                 '--format=${--format_ARRAY[$--format_INDEX]}\n',
                                 'let "TASK_ID=$TASK_ID / 2"\n\n',
                                 'myprog\n'])

    def test_create_jobgroup(self):
        """create dummy jobgroup."""
        jobgroup = pyani_jobs.JobGroup('empty', '')
        assert_equal(jobgroup.script, self.emptyscript)

    def test_1d_jobgroup(self):
        """create dummy 1-parameter sweep jobgroup."""
        jobgroup = pyani_jobs.JobGroup(
            '1d-sweep', 'cat', arguments=self.params1)
        assert_equal(jobgroup.script, self.p1script)
        assert_equal(3, jobgroup.tasks)

    def test_2d_jobgroup(self):
        """create dummy 2-parameter sweep jobgroup."""
        jobgroup = pyani_jobs.JobGroup(
            '2d-sweep', 'myprog', arguments=self.params2)
        print(jobgroup.script)
        print(self.p2script)
        assert_equal(jobgroup.script, self.p2script)
        assert_equal(4, jobgroup.tasks)

    def test_add_dependency(self):
        """add jobgroup dependency."""
        jg1 = pyani_jobs.JobGroup('1d-sweep', 'cat', arguments=self.params1)
        jg2 = pyani_jobs.JobGroup('2d-sweep', 'myprog', arguments=self.params2)
        jg2.add_dependency(jg1)
        assert_equal(4, jg2.tasks)
        assert_equal(1, len(jg2.dependencies))
        dep = jg2.dependencies[0]
        assert_equal(3, dep.tasks)
        assert_equal('1d-sweep', dep.name)

    def test_remove_dependency(self):
        """add and remove jobgroup dependency."""
        jg1 = pyani_jobs.JobGroup('1d-sweep', 'cat', arguments=self.params1)
        jg2 = pyani_jobs.JobGroup('2d-sweep', 'myprog', arguments=self.params2)
        jg2.add_dependency(jg1)
        assert_equal(1, len(jg2.dependencies))
        dep = jg2.dependencies[0]
        assert_equal('1d-sweep', dep.name)
        jg2.remove_dependency(dep)
        assert_equal(0, len(jg2.dependencies))
