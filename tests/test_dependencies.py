#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test_dependencies.py

Tests for availability of pyani dependencies

We only test for dependencies from non-standard libraries.

These tests are intended to be run from the repository root using:

pytest -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.

(c) The James Hutton Institute 2017-2019
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

Copyright (c) 2017-2019 The James Hutton Institute

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


class TestDependencyExecutables(unittest.TestCase):

    """Class defining tests of third-party executables"""

    def test_run_blast(self):
        """Test that BLAST+ is runnable."""
        blastn_exe = pyani_config.BLASTN_DEFAULT
        cmd = [blastn_exe, "-version"]
        result = subprocess.run(
            cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
        print(result.stdout)
        self.assertEqual(result.stdout[:6], b"blastn")

    def test_run_blastall(self):
        """Test that legacy BLAST is runnable."""
        blastall_exe = pyani_config.BLASTALL_DEFAULT
        cmd = blastall_exe
        # Can't use check=True, as blastall without arguments returns 1!
        result = subprocess.run(
            cmd,
            shell=False,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print(result.stdout)
        self.assertEqual(result.stdout[1:9], b"blastall")

    def test_run_nucmer(self):
        """Test that NUCmer is runnable."""
        nucmer_exe = pyani_config.NUCMER_DEFAULT
        cmd = [nucmer_exe, "--version"]
        result = subprocess.run(
            cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
        )
        print(result.stderr)  # NUCmer puts output to STDERR!
        self.assertEqual(result.stderr[:6], b"nucmer")
