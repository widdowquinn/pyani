#!/usr/bin/env python

"""Tests for pyani graphics

These tests are intended to be run using the nose package
(see https://nose.readthedocs.org/en/latest/), from the repository root
directory.

If the test is run directly at the command-line, the output obtained by each
test is returned to STDOUT.
"""

import os
import pandas as pd
import shutil

from nose.tools import assert_equal, assert_less, nottest
from pyani import pyani_graphics, pyani_config, pyani_tools


# Work out where we are. We need to do this to find related data files
# for testing
curdir = os.path.dirname(os.path.abspath(__file__))

OUTDIR = os.path.join("tests", "test_graphics_output")


def define_inputs():
    return {'infilename': os.path.join("tests", "target_ANIm_output",
                                       "ANIm_percentage_identity.tab"),
            'labels': pyani_tools.get_labels(os.path.join("tests",
                                                          "test_ani_data",
                                                          "labels.tab"),
                                             logger=None),
            'classes': pyani_tools.get_labels(os.path.join("tests",
                                                           "test_ani_data",
                                                           "classes.tab"),
                                              logger=None)}
            

def draw_format_method(fmt, mth):
    """Render graphics format and method output."""
    inputs = define_inputs()
    outfilename = os.path.join(OUTDIR, "%s.%s" % (mth, fmt))
    stem = "ANIm_percentage_identity"
    df = pd.read_csv(inputs['infilename'], index_col=0, sep="\t")
    os.makedirs(OUTDIR, exist_ok=True)
    fn = {"mpl": pyani_graphics.heatmap_mpl,
          "seaborn": pyani_graphics.heatmap_seaborn}
    params = {"mpl": pyani_config.params_mpl,
              "seaborn": pyani_config.params_mpl}
    method_params = params[mth](df)[stem]
    fn[mth](df, outfilename, title="%s:%s test" % (mth, fmt),
            params=method_params,
            labels=inputs['labels'], classes=inputs['classes'])


def test_png_mpl():
    """Write .png graphics with mpl"""
    draw_format_method("png", "mpl")


def test_svg_mpl():
    """Write .svg graphics with mpl"""
    draw_format_method("svg", "mpl")


def test_pdf_mpl():
    """Write .pdf graphics with mpl"""
    draw_format_method("pdf", "mpl")


def test_png_seaborn():
    """Write .png graphics with seaborn"""
    draw_format_method("png", "seaborn")


def test_svg_seaborn():
    """Write .svg graphics with seaborn"""
    draw_format_method("svg", "seaborn")


def test_pdf_seaborn():
    """Write .pdf graphics with seaborn"""
    draw_format_method("pdf", "seaborn")
