#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# 161 Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019-2020 University of Strathclyde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
"""test_jobs.py

Test pyani_jobs.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from typing import Dict, NamedTuple

import pytest

from pyani import pyani_jobs


class JobScript(NamedTuple):

    """Convenience struct for job script creation tests."""

    params: Dict[str, str]
    script: str


@pytest.fixture
def job_dummy_cmds():
    """Dummy commands for testing job creation."""
    return ["ls -ltrh", "echo ${PWD}"]


@pytest.fixture
def job_empty_script():
    """Empty script for testing job creation."""
    return 'let "TASK_ID=$SGE_TASK_ID - 1"\n\n\n\n'


@pytest.fixture
def job_scripts():
    """Return two JobScript namedtuples for testing job creation."""
    return (
        JobScript(
            {"-f": ["file1", "file2", "file3"]},
            "".join(
                [
                    'let "TASK_ID=$SGE_TASK_ID - 1"\n',
                    "-f_ARRAY=( file1 file2 file3  )\n\n",
                    'let "-f_INDEX=$TASK_ID % 3"\n',
                    "-f=${-f_ARRAY[$-f_INDEX]}\n",
                    'let "TASK_ID=$TASK_ID / 3"\n\n',
                    "cat\n",
                ]
            ),
        ),
        JobScript(
            {"-f": ["file1", "file2"], "--format": ["fmtA", "fmtB"]},
            "".join(
                [
                    'let "TASK_ID=$SGE_TASK_ID - 1"\n',
                    "--format_ARRAY=( fmtA fmtB  )\n",
                    "-f_ARRAY=( file1 file2  )\n\n",
                    'let "--format_INDEX=$TASK_ID % 2"\n',
                    "--format=${--format_ARRAY[$--format_INDEX]}\n",
                    'let "TASK_ID=$TASK_ID / 2"\n',
                    'let "-f_INDEX=$TASK_ID % 2"\n',
                    "-f=${-f_ARRAY[$-f_INDEX]}\n",
                    'let "TASK_ID=$TASK_ID / 2"\n\n',
                    "myprog\n",
                ]
            ),
        ),
    )


def test_create_job():
    """Create a dummy job."""
    job = pyani_jobs.Job("empty", "")
    assert job.script == ""


def test_create_job_with_command(job_dummy_cmds):
    """Create dummy job with command."""
    job = pyani_jobs.Job("dummy", job_dummy_cmds[0])
    assert job.script == job_dummy_cmds[0]


def test_add_dependency(job_dummy_cmds):
    """Create dummy job with dependency."""
    job1 = pyani_jobs.Job("dummy_with_dependency", job_dummy_cmds[0])
    job2 = pyani_jobs.Job("dummy_dependency", job_dummy_cmds[1])
    job1.add_dependency(job2)
    dep = job1.dependencies[0]

    assert (job_dummy_cmds[0], job_dummy_cmds[1], 1) == (
        job1.script,
        dep.script,
        len(job1.dependencies),
    )


def test_remove_dependency(job_dummy_cmds):
    """Create dummy job, add and remove dependency."""
    job1 = pyani_jobs.Job("dummy_with_dependency", job_dummy_cmds[0])
    job2 = pyani_jobs.Job("dummy_dependency", job_dummy_cmds[1])
    job1.add_dependency(job2)
    dep = job1.dependencies[0]
    job1.remove_dependency(dep)

    assert 0 == len(job1.dependencies)


def test_create_jobgroup(job_empty_script):
    """create dummy jobgroup."""
    jobgroup = pyani_jobs.JobGroup("empty", "")
    assert jobgroup.script == job_empty_script


def test_1d_jobgroup(job_scripts):
    """create dummy 1-parameter sweep jobgroup."""
    jobgroup = pyani_jobs.JobGroup("1d-sweep", "cat", arguments=job_scripts[0].params)
    assert (jobgroup.script, 3) == (job_scripts[0].script, jobgroup.tasks)


def test_2d_jobgroup(job_scripts):
    """create dummy 2-parameter sweep jobgroup."""
    jobgroup = pyani_jobs.JobGroup(
        "2d-sweep", "myprog", arguments=job_scripts[1].params
    )
    assert (jobgroup.script, 4) == (job_scripts[1].script, jobgroup.tasks)


def test_add_group_dependency(job_scripts):
    """add jobgroup dependency."""
    jg1 = pyani_jobs.JobGroup("1d-sweep", "cat", arguments=job_scripts[0].params)
    jg2 = pyani_jobs.JobGroup("2d-sweep", "myprog", arguments=job_scripts[1].params)
    jg2.add_dependency(jg1)
    dep = jg2.dependencies[0]

    assert (1, 3, 4, "1d-sweep") == (
        len(jg2.dependencies),
        dep.tasks,
        jg2.tasks,
        dep.name,
    )


def test_remove_group_dependency(job_scripts):
    """add and remove jobgroup dependency."""
    jg1 = pyani_jobs.JobGroup("1d-sweep", "cat", arguments=job_scripts[0].params)
    jg2 = pyani_jobs.JobGroup("2d-sweep", "myprog", arguments=job_scripts[1].params)
    jg2.add_dependency(jg1)
    dep = jg2.dependencies[0]
    jg2.remove_dependency(dep)

    assert 0 == len(jg2.dependencies)
