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
"""Test for availability of pyani dependencies.

We only test for dependencies from non-standard libraries.

These tests are intended to be run from the repository root using:

pytest -v
"""

import subprocess
import sys
import unittest

import pytest

from pyani import pyani_config


def test_import_biopython():
    """Test Biopython import."""
    import Bio


def test_import_matplotlib():
    """Test matplotlib import."""
    import matplotlib


def test_import_numpy():
    """Test numpy import."""
    import numpy


def test_import_pandas():
    """Test pandas import."""
    import pandas


def test_import_scipy():
    """Test scipy import."""
    import scipy


def test_blastn_available(blastn_available):
    """Test that BLAST+ is available."""
    assert blastn_available


@pytest.mark.xfail(reason="Optional third-party executable (blastall)")
def test_run_blastall(blastall_available):
    """Test that blastall is available."""
    assert blastall_available


def test_run_nucmer(nucmer_available):
    """Test that nucmer is available."""
    assert nucmer_available
