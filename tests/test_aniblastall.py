#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) University of Strathclyde 2019-2021
# Author: Bailey Harrington
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
# Copyright (c) 2019-2021 University of Strathclyde
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
"""Test aniblastall.py module.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pathlib import Path
import unittest
from pyani import aniblastall

# Create object for accessing unittest assertions
assertions = unittest.TestCase("__init__")


# Test get_version()
# Test case 0: no executable location is specified
def test_get_version_nonetype():
    """Test behaviour when no location for the executable is given."""
    test_file_0 = None

    assert (
        aniblastall.get_version(test_file_0)
        == f"expected path to blastall executable; received {test_file_0}"
    )


# Test case 1: no such file exists
def test_get_version_random_string():
    """Test behaviour when the given 'file' is not one."""
    test_file_1 = "string"

    assert (
        aniblastall.get_version(test_file_1) == f"{test_file_1} is not found in $PATH"
    )


# Test case 2: there is no executable
def test_get_version_missing_exe(executable_missing):
    """Test behaviour when there is no file at the specified executable location."""
    test_file_2 = Path("/non/existent/blastall")
    assert aniblastall.get_version(test_file_2) == f"No blastall at {test_file_2}"


# Test case 3: there is a file, but it is not executable
def test_get_version_not_executable(executable_not_executable):
    """Test behaviour when the file at the executable location is not executable."""
    test_file_3 = Path("/non/executable/blastall")
    assert (
        aniblastall.get_version(test_file_3)
        == f"blastall exists at {test_file_3} but not executable"
    )


# Test case 4: there is an executable file, but the version can't be retrieved
def test_get_version_no_version(executable_without_version):
    """Test behaviour when the version for the executable can not be retrieved."""
    test_file_4 = Path("/missing/version/blastall")
    assert (
        aniblastall.get_version(test_file_4)
        == f"blastall exists at {test_file_4} but could not retrieve version"
    )


# Test case 4: there is an executable file, but it will not run on the OS
def test_get_version_os_incompatible(executable_incompatible_with_os):
    """Test behaviour when the program can't run on the operating system.
    This will happen with newer versions of MacOS."""
    test_file_4 = Path("/os/incompatible/blastall")
    assert (
        aniblastall.get_version(test_file_4)
        == f"blastall exists at {test_file_4} but could not be executed"
    )
