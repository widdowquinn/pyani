#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) The University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
# Cathedral Street
# Glasgow
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2018 The James Hutton Institute
# (c) The University of Strathclyde 2019
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
"""Code to implement graphics output for ANI analyses."""

# Force matplotlib NOT to use an Xwindows backend on *nix, so that
# _tkinter.TclError is avoided when there is no $DISPLAY env: this can occur
# when running the package/script via ssh
# See http://stackoverflow.com/questions/2801882/\
#            generating-a-png-with-matplotlib-when-display-is-undefined
# This needs to be done before importing pyplot

import matplotlib  # pylint: disable=C0411

from . import mpl  # matplotlib wrappers
from . import sns  # seaborn wrappers

# Specify matplotlib backend. This *must* be done before pyplot import, but
# raises errors with flake8 etc. So we comment out the specific error
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402 # pylint: disable=wrong-import-position,wrong-import-order


# Convenience class to hold heatmap graphics parameters
class Params:  # pylint: disable=too-few-public-methods

    """Convenience class to hold heatmap rendering parameters."""

    def __init__(self, params, labels=None, classes=None):
        """Instantiate class.

        :param params:
        :param labels:
        :param classes:
        """
        self.cmap = plt.get_cmap(params[0])
        self.vmin = params[1]
        self.vmax = params[2]
        self.labels = labels
        self.classes = classes

    @property
    def vdiff(self):
        """Return difference between max and min values for presentation."""
        return max(0.01, self.vmax - self.vmin)


# helper for cleaning up matplotlib axes by removing ticks etc.
def clean_axis(axis):
    """Remove ticks, tick labels, and frame from axis.

    :param axis:
    """
    axis.get_xaxis().set_ticks([])
    axis.get_yaxis().set_ticks([])
    for spine in list(axis.spines.values()):
        spine.set_visible(False)
