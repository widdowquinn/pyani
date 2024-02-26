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
import logging

from pathlib import Path

import pytest
from typing import Dict, List, NamedTuple, Tuple
import subprocess
import pickle
import multiprocessing

from pyani.anib import fragment_fasta_files, make_blastcmd_builder, make_job_graph
from pyani.pyani_jobs import Job
from pyani.run_multiprocessing import (
    multiprocessing_run,
    populate_cmdsets,
    run_dependency_graph,
)

assertions = unittest.TestCase("__init__")
logger = logging.getLogger(__name__)


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


# Convenience struct wich mocks a multiprocessing.pool worker's task
# pool.apply_async returns a multiprocessing.pool.ApplyResult object,
# which is essentially a NamedTuple, but has a get method, hence
# the need to define one here that doesn't appear to do much
class MockProcess(NamedTuple):
    returncode: int
    args: List[str]
    stdout: str
    stderr: str

    def get(self):
        return self


# This must be a top-level function so that the pool workers can pickle it
def mock_run(*args, **kwargs):
    """Mock a call to `subprocess.run()`."""
    return MockProcess(1, ["ls", "-ltrh"], b"mock bytes", b"mock bytes")


@pytest.fixture
def bad_comparison(monkeypatch):
    """
    Mocks a call to `run_multiprocessing.multiprocessing_run()` that fails.
    """
    pool = multiprocessing.Pool(processes=8)

    def mock_apply_async(*args, **kwargs):
        """Mock a call to `multiprocessing.Pool().apply_async()`."""
        return multiprocessing.reduction.ForkingPickler.dumps(
            MockProcess(1, ["ls", "-ltrh"], b"mock bytes", b"mock bytes")
        )

    monkeypatch.setattr(pool, "apply_async", mock_apply_async)
    monkeypatch.setattr(subprocess, "run", mock_run)


def test_failed_comparison(mp_cmdlist, bad_comparison, monkeypatch):
    """Test that a warning is logged when an individual comparison fails."""

    with assertions.assertLogs(logger, level="WARNING") as scribe:
        multiprocessing_run(mp_cmdlist, logger=logger)
        assertions.assertEqual(
            scribe.output[0], "WARNING:test_multiprocessing:Comparison failed: ls -ltrh"
        )
