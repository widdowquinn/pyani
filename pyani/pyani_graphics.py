#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) The University of Strathclude 2019
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
# (c) The University of Strathclude 2019
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

import warnings

from math import floor, log10

import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as distance
import seaborn as sns
import pandas as pd

from scipy.stats import gaussian_kde

from . import pyani_config

# Specify matplotlib backend. This *must* be done before pyplot import, but
# raises errors with flake8 etc. So we comment out the specific error
import matplotlib  # pylint: disable=C0411

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402 # pylint: disable=wrong-import-position,wrong-import-order
import matplotlib.gridspec as gridspec  # noqa: E402 # pylint: disable=wrong-import-position,wrong-import-order


# Register Matplotlib colourmaps
plt.register_cmap(cmap=pyani_config.CMAP_SPBND_BURD)
plt.register_cmap(cmap=pyani_config.CMAP_HADAMARD_BURD)
plt.register_cmap(cmap=pyani_config.CMAP_BURD)

# Matplotlib version dictates bug fixes
MPLVERSION = matplotlib.__version__


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
    """Remove ticks, tick labels, and frame from axis."""
    axis.get_xaxis().set_ticks([])
    axis.get_yaxis().set_ticks([])
    for spine in list(axis.spines.values()):
        spine.set_visible(False)


# Add classes colorbar to Seaborn plot
def get_seaborn_colorbar(dfr, classes):
    """Return a colorbar representing classes, for a Seaborn plot.

    The aim is to get a pd.Series for the passed dataframe columns,
    in the form:
    0    colour for class in col 0
    1    colour for class in col 1
    ...  colour for class in col ...
    n    colour for class in col n
    """
    levels = sorted(list(set(classes.values())))
    paldict = {
        lvl: pal
        for (lvl, pal) in zip(
            levels,
            sns.cubehelix_palette(
                len(levels), light=0.9, dark=0.1, reverse=True, start=1, rot=-2
            ),
        )
    }
    lvl_pal = {cls: paldict[lvl] for (cls, lvl) in list(classes.items())}
    # Have to use string conversion of the dataframe index, here
    col_cb = pd.Series([str(_) for _ in dfr.index]).map(lvl_pal)
    # The col_cb Series index now has to match the dfr.index, but
    # we don't create the Series with this (and if we try, it
    # fails) - so change it with this line
    col_cb.index = dfr.index
    return col_cb


# Add labels to the seaborn heatmap axes
def add_seaborn_labels(fig, params):
    """Add labels to Seaborn heatmap axes, in-place."""
    if params.labels:
        # If a label mapping is missing, use the key text as fall back
        [
            _.set_text(params.labels.get(_.get_text(), _.get_text()))
            for _ in fig.ax_heatmap.get_yticklabels()
        ]
        [
            _.set_text(params.labels.get(_.get_text(), _.get_text()))
            for _ in fig.ax_heatmap.get_xticklabels()
        ]
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xticklabels(), rotation=90)
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_yticklabels(), rotation=0)
    return fig


# Return a clustermap
def get_seaborn_clustermap(dfr, params, title=None, annot=True):
    """Return a Seaborn clustermap for the passed dataframe."""
    fig = sns.clustermap(
        dfr,
        cmap=params.cmap,
        vmin=params.vmin,
        vmax=params.vmax,
        col_colors=params.colorbar,
        row_colors=params.colorbar,
        figsize=(params.figsize, params.figsize),
        linewidths=params.linewidths,
        annot=annot,
    )

    # add labels for each of the input genomes
    add_seaborn_labels(fig, params)

    fig.cax.yaxis.set_label_position("left")
    if title:
        fig.cax.set_ylabel(title)

    # Return clustermap
    return fig


# Generate Seaborn heatmap output
def heatmap_seaborn(dfr, outfilename=None, title=None, params=None):
    """Return seaborn heatmap with cluster dendrograms.

    - dfr - pandas DataFrame with relevant data
    - outfilename - path to output file (indicates output format)
    """
    # Decide on figure layout size: a minimum size is required for
    # aesthetics, and a maximum to avoid core dumps on rendering.
    # If we hit the maximum size, we should modify font size.
    maxfigsize = 120
    calcfigsize = dfr.shape[0] * 1.1
    figsize = min(max(8, calcfigsize), maxfigsize)
    if figsize == maxfigsize:
        scale = maxfigsize / calcfigsize
        sns.set_context("notebook", font_scale=scale)

    # Add a colorbar?
    if params.classes is None:
        col_cb = None
    else:
        col_cb = get_seaborn_colorbar(dfr, params.classes)

    # Add attributes to parameter object, and draw heatmap
    params.colorbar = col_cb
    params.figsize = figsize
    params.linewidths = 0.25
    fig = get_seaborn_clustermap(dfr, params, title=title)

    # Save to file
    if outfilename:
        fig.savefig(outfilename)

    # Return clustermap
    return fig


def distribution_seaborn(dfr, outfilename, matname, title=None):
    """Return seaborn distribution plot for matrix.

    :param drf:  DataFrame with results matrix
    :param outfilename:  Path to output file for writing
    :param matname:  str, type of matrix being plotted
    :param title:  str, optional title
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(title)
    sns.distplot(
        dfr.values.flatten(), kde=False, rug=False, ax=axes[0], norm_hist=False
    )
    sns.distplot(
        dfr.values.flatten(), hist=False, rug=True, ax=axes[1], norm_hist=False
    )

    # Modify axes after data is plotted
    for ax in axes:
        if matname == "sim_errors":
            ax.set_xlim(0, ax.get_xlim()[1])
        elif matname in ["hadamard", "coverage"]:
            ax.set_xlim(0, 1.01)
        elif matname == "identity":
            ax.set_xlim(0.75, 1.01)

    # Tidy figure
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if outfilename:
        # For some reason seaborn gives us an AxesSubPlot with
        # sns.distplot, rather than a Figure, so we need this hack
        fig.savefig(outfilename)

    return fig


# Add dendrogram and axes to passed figure
def add_mpl_dendrogram(dfr, fig, params, heatmap_gs, orientation="col"):
    """Return a dendrogram and corresponding gridspec, attached to the fig.

    Modifies the fig in-place. Orientation is either 'row' or 'col' and
    determines location and orientation of the rendered dendrogram.

    We expect that the row/column index values should be ordered
    identically. If they are not, the dendrogram will not match labels
    """
    # Row or column axes?
    if orientation == "row":
        dists = distance.pdist(dfr)
        spec = heatmap_gs[1, 0]
        orient = "left"
        nrows, ncols = 1, 2
        height_ratios = [1]
    else:  # Column dendrogram
        dists = distance.pdist(dfr.T)
        spec = heatmap_gs[0, 1]
        orient = "top"
        nrows, ncols = 2, 1
        height_ratios = [1, 0.15]

    # Create row dendrogram axis
    gspec = gridspec.GridSpecFromSubplotSpec(
        nrows,
        ncols,
        subplot_spec=spec,
        wspace=0.0,
        hspace=0.1,
        height_ratios=height_ratios,
    )
    dend_axes = fig.add_subplot(gspec[0, 0])
    dend = sch.dendrogram(
        sch.linkage(distance.squareform(dists), method="complete"),
        color_threshold=np.inf,
        orientation=orient,
        labels=list(params.labels.values()),
        get_leaves=True,
    )
    clean_axis(dend_axes)

    return {"dendrogram": dend, "gridspec": gspec}


def distribution_mpl(dfr, outfilename, matname, title=None):
    """Return matplotlib distribution plot for matrix.

    :param drf:  DataFrame with results matrix
    :param outfilename:  Path to output file for writing
    :param matname:  str, type of matrix being plotted
    :param title:  str, optional title
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(title)
    data = dfr.values.flatten()
    xs = np.linspace(min(data), max(data), 200)
    # Plot histogram
    axes[0].hist(data, bins=50)
    # Plot density
    density = gaussian_kde(data)
    density._compute_covariance()
    axes[1].plot(xs, density(xs))

    # Modify axes after data is plotted
    for ax in axes:
        if matname == "sim_errors":
            ax.set_xlim(0, ax.get_xlim()[1])
        elif matname in ["hadamard", "coverage"]:
            ax.set_xlim(0, 1.01)
        elif matname == "identity":
            ax.set_xlim(ax.get_xlim()[0], 1.01)

    # Tidy figure
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if outfilename:
        # For some reason seaborn gives us an AxesSubPlot with
        # sns.distplot, rather than a Figure, so we need this hack
        fig.savefig(outfilename)

    return fig


# Create heatmap axes for Matplotlib output
def get_mpl_heatmap_axes(dfr, fig, heatmap_gs):
    """Return axis for Matplotlib heatmap."""
    # Create heatmap axis
    heatmap_axes = fig.add_subplot(heatmap_gs[1, 1])
    heatmap_axes.set_xticks(np.linspace(0, dfr.shape[0] - 1, dfr.shape[0]))
    heatmap_axes.set_yticks(np.linspace(0, dfr.shape[0] - 1, dfr.shape[0]))
    heatmap_axes.grid(False)
    heatmap_axes.xaxis.tick_bottom()
    heatmap_axes.yaxis.tick_right()
    return heatmap_axes


def add_mpl_colorbar(dfr, fig, dend, params, orientation="row"):
    """Add class colorbars to Matplotlib heatmap."""
    # Assign a numerical value to each class, for mpl
    classdict = {cls: idx for (idx, cls) in enumerate(params.classes.values())}

    # colourbar
    cblist = []
    for name in [str(_) for _ in dfr.index[dend["dendrogram"]["leaves"]]]:
        try:
            cblist.append(classdict[params.classes[name]])
        except KeyError:
            cblist.append(classdict[name])
    colbar = pd.Series(cblist)

    # Create colourbar axis - could capture if needed
    if orientation == "row":
        cbaxes = fig.add_subplot(dend["gridspec"][0, 1])
        cbaxes.imshow(
            [[cbar] for cbar in colbar.values],
            cmap=plt.get_cmap(pyani_config.MPL_CBAR),
            interpolation="nearest",
            aspect="auto",
            origin="lower",
        )
    else:
        cbaxes = fig.add_subplot(dend["gridspec"][1, 0])
        cbaxes.imshow(
            [colbar],
            cmap=plt.get_cmap(pyani_config.MPL_CBAR),
            interpolation="nearest",
            aspect="auto",
            origin="lower",
        )
    clean_axis(cbaxes)
    return colbar


# Add labels to the heatmap axes
def add_mpl_labels(heatmap_axes, rowlabels, collabels, params):
    """Add labels to Matplotlib heatmap axes, in-place."""
    if params.labels:
        # If a label mapping is missing, use the key text as fall back
        rowlabels = [params.labels.get(lab, lab) for lab in rowlabels]
        collabels = [params.labels.get(lab, lab) for lab in collabels]
    xlabs = heatmap_axes.set_xticklabels(collabels)
    ylabs = heatmap_axes.set_yticklabels(rowlabels)
    for label in xlabs:  # Rotate column labels
        label.set_rotation(90)
    for labset in (xlabs, ylabs):  # Smaller font
        for label in labset:
            label.set_fontsize(8)


# Add colour scale to heatmap
def add_mpl_colorscale(fig, heatmap_gs, ax_map, params, title=None):
    """Add colour scale to heatmap."""
    # Set tick intervals
    cbticks = [params.vmin + e * params.vdiff for e in (0, 0.25, 0.5, 0.75, 1)]
    if params.vmax > 10:
        exponent = int(floor(log10(params.vmax))) - 1
        cbticks = [int(round(e, -exponent)) for e in cbticks]

    scale_subplot = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=heatmap_gs[0, 0], wspace=0.0, hspace=0.0
    )
    scale_ax = fig.add_subplot(scale_subplot[0, 1])
    cbar = fig.colorbar(ax_map, scale_ax, ticks=cbticks)
    if title:
        cbar.set_label(title, fontsize=6)
    cbar.ax.yaxis.set_ticks_position("left")
    cbar.ax.yaxis.set_label_position("left")
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0)
    return cbar


# Generate Matplotlib heatmap output
def heatmap_mpl(dfr, outfilename=None, title=None, params=None):
    """Return matplotlib heatmap with cluster dendrograms.

    - dfr - pandas DataFrame with relevant data
    - outfilename - path to output file (indicates output format)
    - params - a list of parameters for plotting: [colormap, vmin, vmax]
    - labels - dictionary of alternative labels, keyed by default sequence
               labels
    - classes - dictionary of sequence classes, keyed by default sequence
                labels
    """
    # Sort rows by index - this ensures that labels match the dendrogram.
    # When recovering dataframes from the database, we get row
    # indexes/labels as integers, but out of order, resulting in a
    # mismatch of labels to leaves in the dendrogram. This line remedies
    # that.
    dfr = dfr.sort_index()

    # Layout figure grid and add title
    # Set figure size by the number of rows in the dataframe
    figsize = max(8, dfr.shape[0] * 0.175)
    fig = plt.figure(figsize=(figsize, figsize))
    # if title:
    #     fig.suptitle(title)
    heatmap_gs = gridspec.GridSpec(
        2, 2, wspace=0.0, hspace=0.0, width_ratios=[0.3, 1], height_ratios=[0.3, 1]
    )

    # Add column and row dendrograms/axes to figure
    coldend = add_mpl_dendrogram(dfr, fig, params, heatmap_gs, orientation="col")
    rowdend = add_mpl_dendrogram(dfr, fig, params, heatmap_gs, orientation="row")

    # Add heatmap axes to figure, with rows/columns as in the dendrograms
    heatmap_axes = get_mpl_heatmap_axes(dfr, fig, heatmap_gs)
    ax_map = heatmap_axes.imshow(
        dfr.iloc[rowdend["dendrogram"]["leaves"], coldend["dendrogram"]["leaves"]],
        interpolation="nearest",
        cmap=params.cmap,
        origin="lower",
        vmin=params.vmin,
        vmax=params.vmax,
        aspect="auto",
    )

    # Are there class colourbars to add?
    if params.classes is not None:
        add_mpl_colorbar(dfr, fig, coldend, params, orientation="col")
        add_mpl_colorbar(dfr, fig, rowdend, params, orientation="row")

    # Add heatmap labels
    add_mpl_labels(
        heatmap_axes, rowdend["dendrogram"]["ivl"], coldend["dendrogram"]["ivl"], params
    )

    # Add colour scale
    add_mpl_colorscale(fig, heatmap_gs, ax_map, params, title)

    # Return figure output, and write, if required
    plt.subplots_adjust(top=0.85)  # Leave room for title
    # fig.set_tight_layout(True)
    # We know that there is a UserWarning here about tight_layout and
    # using the Agg renderer on OSX, so catch and ignore it, for cleanliness.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        heatmap_gs.tight_layout(fig, h_pad=0.1, w_pad=0.5)
    if outfilename:
        fig.savefig(outfilename)
    return fig
