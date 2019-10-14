#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019 University of Strathclyde
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
"""Tests for pyani graphics.

These tests are intended to be run from the repository root using:

pytest -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.
"""

import os
import pandas as pd

from pyani import pyani_graphics, pyani_config, pyani_tools


# Work out where we are. We need to do this to find related data files
# for testing
curdir = os.path.dirname(os.path.abspath(__file__))

OUTDIR = os.path.join("tests", "test_graphics_output")


def define_inputs():
    return {
        "infilename": os.path.join(
            "tests", "target_ANIm_output", "ANIm_percentage_identity.tab"
        ),
        "labels": pyani_tools.get_labels(
            os.path.join("tests", "test_ani_data", "labels.tab"), logger=None
        ),
        "classes": pyani_tools.get_labels(
            os.path.join("tests", "test_ani_data", "classes.tab"), logger=None
        ),
    }


def draw_format_method(fmt, mth):
    """Render graphics format and method output."""
    inputs = define_inputs()
    outfilename = os.path.join(OUTDIR, "%s.%s" % (mth, fmt))
    stem = "ANIm_percentage_identity"
    df = pd.read_csv(inputs["infilename"], index_col=0, sep="\t")
    os.makedirs(OUTDIR, exist_ok=True)
    fn = {"mpl": pyani_graphics.heatmap_mpl, "seaborn": pyani_graphics.heatmap_seaborn}
    params = {"mpl": pyani_config.params_mpl, "seaborn": pyani_config.params_mpl}
    method_params = pyani_graphics.Params(
        params[mth](df)[stem], inputs["labels"], inputs["classes"]
    )
    fn[mth](df, outfilename, title=f"{mth}:{fmt} test", params=method_params)


def test_png_mpl():
    """Write .png graphics with mpl."""
    draw_format_method("png", "mpl")


def test_svg_mpl():
    """Write .svg graphics with mpl."""
    draw_format_method("svg", "mpl")


def test_pdf_mpl():
    """Write .pdf graphics with mpl."""
    draw_format_method("pdf", "mpl")


def test_png_seaborn():
    """Write .png graphics with seaborn."""
    draw_format_method("png", "seaborn")


def test_svg_seaborn():
    """Write .svg graphics with seaborn."""
    draw_format_method("svg", "seaborn")


def test_pdf_seaborn():
    """Write .pdf graphics with seaborn."""
    draw_format_method("pdf", "seaborn")
