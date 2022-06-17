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
"""Tests for pyani graphics.

These tests are intended to be run from the repository root using:

pytest -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.
"""

import pandas as pd

from pathlib import Path
from typing import Dict, NamedTuple

import pytest

from pyani import pyani_config, pyani_graphics
from pyani.pyani_tools import get_labels


class GraphicsTestInputs(NamedTuple):

    """Convenience struct for graphics test inputs."""

    filename: Path
    labels: Dict[str, str]
    classes: Dict[str, str]


@pytest.fixture
def graphics_inputs(dir_graphics_in):
    """Returns namedtuple of graphics inputs."""
    return GraphicsTestInputs(
        dir_graphics_in / "ANIm_percentage_identity.tab",
        get_labels(dir_graphics_in / "labels.tab"),
        get_labels(dir_graphics_in / "classes.tab"),
    )


def draw_format_method(fmt, mth, graphics_inputs, tmp_path):
    """Render graphics format and method output."""
    df = pd.read_csv(graphics_inputs.filename, index_col=0, sep="\t")
    fn = {"mpl": pyani_graphics.mpl.heatmap, "seaborn": pyani_graphics.sns.heatmap}
    sc = {"mpl": pyani_graphics.mpl.scatter, "seaborn": pyani_graphics.sns.scatter}
    params = {"mpl": pyani_config.params_mpl, "seaborn": pyani_config.params_mpl}
    method_params = pyani_graphics.Params(
        params[mth](df)["ANIm_percentage_identity"],
        graphics_inputs.labels,
        graphics_inputs.classes,
    )
    fn[mth](
        df, tmp_path / f"{mth}.{fmt}", title=f"{mth}:{fmt} test", params=method_params
    )
    sc[mth](
        df,
        df,
        tmp_path / f"{mth}.{fmt}",
        "matrix1",
        "matrix2",
        title=f"{mth}:{fmt} test",
        params=method_params,
    )


def test_png_mpl(graphics_inputs, tmp_path):
    """Write .png graphics with mpl."""
    draw_format_method("png", "mpl", graphics_inputs, tmp_path)


def test_svg_mpl(graphics_inputs, tmp_path):
    """Write .svg graphics with mpl."""
    draw_format_method("svg", "mpl", graphics_inputs, tmp_path)


def test_pdf_mpl(graphics_inputs, tmp_path):
    """Write .pdf graphics with mpl."""
    draw_format_method("pdf", "mpl", graphics_inputs, tmp_path)


def test_png_seaborn(graphics_inputs, tmp_path):
    """Write .png graphics with seaborn."""
    draw_format_method("png", "seaborn", graphics_inputs, tmp_path)


def test_svg_seaborn(graphics_inputs, tmp_path):
    """Write .svg graphics with seaborn."""
    draw_format_method("svg", "seaborn", graphics_inputs, tmp_path)


def test_pdf_seaborn(graphics_inputs, tmp_path):
    """Write .pdf graphics with seaborn."""
    draw_format_method("pdf", "seaborn", graphics_inputs, tmp_path)
