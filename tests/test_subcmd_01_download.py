#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) The University of Strathclude 2019-2024
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
# Copyright (c) 2016-2019 The James Hutton Institute
# (c) The University of Strathclude 2019-2024
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
"""Test download subcommand for pyani.

The test suite is intended to be run from the repository root using:

pytest -v

Each command CMD available at the command line as pyani <CMD> is
tested in its own class as a subclass of unittest.TestCase, where
setUp() defines input/output files, a null logger (which is also
picked up by pytest), and a dictionary of command lines, keyed
by test name, with values representing command-line options.

For each test, command line options are defined in a Namespace and
passed as the sole argument to the appropriate subcommand.

As the download operations are slow, and subject to network issues,
especially with CI, we mock the download operations.
"""

import logging

from argparse import Namespace
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import pytest

from pyani import download
from pyani.scripts import subcommands
from tools import modify_namespace


@pytest.fixture
def api_keypath():
    """NCBI API key path."""
    return Path("~/.ncbi/api_key")


@pytest.fixture
def base_download_namespace(api_keypath, path_fixtures_base, tmp_path):
    """Base namespace for the pyani download subcommand."""
    return Namespace(
        outdir=tmp_path / "C_blochmannia",
        taxon="203804",
        email="my.email@my.domain",
        retries=20,
        batchsize=10000,
        timeout=10,
        classfname="classes.txt",
        labelfname="labels.txt",
        kraken=False,
        force=True,
        noclobber=False,
        dryrun=False,
        disable_tqdm=True,
        api_keypath=api_keypath,
    )


@pytest.fixture
def dryrun_namespace(base_download_namespace):
    """Namespace for pyani download dry run."""
    return modify_namespace(base_download_namespace, dryrun=True)


@pytest.fixture
def kraken_namespace(base_download_namespace, tmp_path):
    """Namespace for downloading C. blochmannia with Kraken labels."""
    return modify_namespace(
        base_download_namespace, kraken=True, outdir=tmp_path / "kraken"
    )


# Create object for accessing unittest assertions
assertions = TestCase("__init__")


def test_create_hash():
    """Test that the expected exception is raised if the file doesn't exist."""
    test_file = "/this/is/not/a/file"
    with assertions.assertRaises(download.PyaniIndexException):
        download.create_hash(test_file)


def test_download_dry_run(dryrun_namespace):
    """Dry run of C. blochmannia download.

    Remote database access is mocked for this test as the number of genomes
    increased over time, and the test would be slow and subject to network
    issues.
    """
    subcommands.subcmd_download(dryrun_namespace)


def test_download_c_blochmannia(base_download_namespace, mock_blochmannia_dl):
    """Test C. blochmannia download."""
    subcommands.subcmd_download(base_download_namespace)


def test_download_kraken(kraken_namespace, mock_blochmannia_kraken_dl):
    """C. blochmannia download in Kraken format."""
    subcommands.subcmd_download(kraken_namespace)
