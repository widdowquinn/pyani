#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
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
# Copyright (c) 2017-2019 The James Hutton Institute
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
"""Test command-line parsing for pyani.

This test suite is intended to be run from the repository root using:

pytest -v
"""

import logging

from pathlib import Path

import pytest

from pyani import pyani_orm
from pyani.scripts import pyani_script


@pytest.fixture
def args_createdb(tmp_path):
    """Command-line arguments for database creation."""
    return ["createdb", "--dbpath", tmp_path / "pyanidb", "--force"]


@pytest.fixture
def args_single_genome_download(email_address, tmp_path):
    """Command-line arguments for single genome download."""
    return [
        "download",
        "-t",
        "218491",
        "--email",
        email_address,
        tmp_path,
        "--force",
    ]


def test_createdb(args_createdb, monkeypatch):
    """Create empty test database."""

    def mock_return_none(*args, **kwargs):
        return None

    monkeypatch.setattr(pyani_orm, "create_db", mock_return_none)
    pyani_script.run_main(args_createdb)


def test_download_single_genome(args_single_genome_download, mock_single_genome_dl):
    """Download a single genome.

    We mock the remote database access
    """
    pyani_script.run_main(args_single_genome_download)
