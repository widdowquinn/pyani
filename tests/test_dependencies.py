#!/usr/bin/env python

"""Tests for availability of pyani dependencies

We only test for dependencies from non-standard libraries.

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/).
"""

import subprocess
import sys
import unittest

import pytest


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
        cmd = "blastn -version"
        result = subprocess.run(
            cmd,
            shell=sys.platform != "win32",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        print(result.stdout)
        self.assertEqual(result.stdout[:6], b"blastn")

    @pytest.mark.skip(reason="Deprecate legacy BLAST")
    def test_run_blastall(self):
        """Test that legacy BLAST is runnable."""
        cmd = "blastall"
        # Can't use check=True, as blastall without arguments returns 1!
        result = subprocess.run(
            cmd,
            shell=sys.platform != "win32",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print(result.stdout)
        self.assertEqual(result.stdout[1:9], b"blastall")

    def test_run_nucmer(self):
        """Test that NUCmer is runnable."""
        cmd = "nucmer --version"
        result = subprocess.run(
            cmd,
            shell=sys.platform != "win32",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        print(result.stderr)  # NUCmer puts output to STDERR!
        self.assertEqual(result.stderr[:6], b"nucmer")
