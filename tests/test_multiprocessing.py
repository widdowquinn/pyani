#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) The University of Strathclude 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# 161 Cathedral Street
# Glasgow
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
# (c) The University of Strathclude 2019-2020
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
"""Test run_multiprocessing.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

import unittest

from pathlib import Path

import pytest

from pyani.anib import fragment_fasta_files, make_blastcmd_builder, make_job_graph
from pyani.pyani_jobs import Job
from pyani.run_multiprocessing import (
    multiprocessing_run,
    populate_cmdsets,
    run_dependency_graph,
)


@pytest.fixture
def mp_cmdlist():
    """List of commands to be run using the multiprocessing interface."""
    return [
        'for i in %s; do echo "Thread %d: value ${i}"; done'
        % (" ".join([str(e) for e in range(v)]), v)
        for v in range(5)
    ]


@pytest.fixture
def mp_dummy_cmds():
    """Dummy commands for building command sets."""
    return ["ls -ltrh", "echo ${PWD}"]


def test_multiprocessing_run(mp_cmdlist):
    """Test that multiprocessing() runs basic jobs."""
    result = multiprocessing_run(mp_cmdlist)
    assert 0 == result


def test_cmdsets(mp_dummy_cmds):
    """Test that module builds command sets."""
    job1 = Job("dummy_with_dependency", mp_dummy_cmds[0])
    job2 = Job("dummy_dependency", mp_dummy_cmds[1])
    job1.add_dependency(job2)
    cmdsets = populate_cmdsets(job1, list(), depth=1)
    target = [{cmd} for cmd in mp_dummy_cmds]
    assert cmdsets == target


@pytest.mark.skip_if_exe_missing("blastn")
def test_dependency_graph_run(path_fna_two, fragment_length, tmp_path):
    """Test that module runs dependency graph."""
    fragresult = fragment_fasta_files(path_fna_two, tmp_path, fragment_length)
    blastcmds = make_blastcmd_builder(tmp_path)
    jobgraph = make_job_graph(path_fna_two, fragresult[0], blastcmds)
    result = run_dependency_graph(jobgraph)
    assert 0 == result
