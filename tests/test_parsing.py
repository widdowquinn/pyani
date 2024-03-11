#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
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
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2029-2024 The University of Strathclude
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
"""Tests for pyani package intermediate file parsing.

These tests are intended to be run from the repository root using:

pytest -v
"""

from pyani.anim import parse_delta


def test_anim_delta(dir_anim_in):
    """Test parsing of NUCmer delta file.
    
    We expect (reference alignment length, query alignment 
    length, weighted average identity, and similarity error 
    count)
    """
    # TODO: allow for rounding error in the returned percentage
    #       identity (this may vary depending on machine/Python
    #       version)
    assert (4016947, 4017751, 0.9994572767031752, 2191) == parse_delta(dir_anim_in / "NC_002696_vs_NC_011916.delta")
