#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
"""Code to implement MatPLotLib graphics output for ANI analyses."""

import warnings

from math import floor, log10

import matplotlib  # pylint: disable=C0411
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as distance

from scipy.stats import gaussian_kde

from pyani import pyani_config

# Specify matplotlib backend. This *must* be done before pyplot import, but
# raises errors with flake8 etc. So we comment out the specific error
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402,E501 # pylint: disable=wrong-import-position,wrong-import-order,ungrouped-imports
import matplotlib.gridspec as gridspec  # noqa: E402,E501 # pylint: disable=wrong-import-position,wrong-import-order,ungrouped-imports


# Register Matplotlib colourmaps
plt.register_cmap(cmap=pyani_config.CMAP_SPBND_BURD)
plt.register_cmap(cmap=pyani_config.CMAP_HADAMARD_BURD)
plt.register_cmap(cmap=pyani_config.CMAP_BURD)

# Matplotlib version dictates bug fixes
MPLVERSION = matplotlib.__version__


# helper for cleaning up matplotlib axes by removing ticks etc.
def clean_axis(axis):
    """Remove ticks, tick labels, and frame from axis.

    :param axis:
    """
    axis.get_xaxis().set_ticks([])
    axis.get_yaxis().set_ticks([])
    for spine in list(axis.spines.values()):
        spine.set_visible(False)


# Add dendrogram and axes to passed figure
def add_dendrogram(dfr, fig, params, heatmap_gs, orientation="col"):
    """Return a dendrogram and corresponding gridspec, attached to the fig.

    :param dfr:  Pandas DataFrame describing input data
    :param fig:  matplotlib Fig that holds graphical output
    :param params:  pyani_graphics.Params object
    :param heatmap_gs:  matplotlib GridSpec for this dendrogram
    :param orientation:  str, "row" or "col"

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
    if len(list(params.labels.values())) == 0:
        labels = None
    else:
        labels = list(params.labels.values())
    dend = sch.dendrogram(
        sch.linkage(dists, method="complete"),
        color_threshold=np.inf,
        orientation=orient,
        labels=labels,
        get_leaves=True,
    )
    clean_axis(dend_axes)

    return {"dendrogram": dend, "gridspec": gspec}


def distribution(dfr, outfilename, matname, title=None):
    """Return matplotlib distribution plot for matrix.

    :param dfr:  DataFrame with results matrix
    :param outfilename:  Path to output file for writing
    :param matname:  str, type of matrix being plotted
    :param title:  str, optional title
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(title)
    data = dfr.values.flatten()
    xvals = np.linspace(min(data), max(data), 200)
    # Plot histogram
    axes[0].hist(data, bins=50)
    # Plot density
    density = gaussian_kde(data)
    density._compute_covariance()  # pylint: disable=protected-access
    axes[1].plot(xvals, density(xvals))

    # Modify axes after data is plotted
    for _ in axes:
        if matname == "sim_errors":
            _.set_xlim(0, _.get_xlim()[1])
        elif matname in ["hadamard", "coverage"]:
            _.set_xlim(0, 1.01)
        elif matname == "identity":
            _.set_xlim(_.get_xlim()[0], 1.01)

    # Tidy figure
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if outfilename:
        # For some reason seaborn gives us an AxesSubPlot with
        # sns.distplot, rather than a Figure, so we need this hack
        fig.savefig(outfilename)

    return fig


# Create heatmap axes for Matplotlib output
def get_heatmap_axes(dfr, fig, heatmap_gs):
    """Return axis for Matplotlib heatmap.

    :param dfr:
    :param fig:
    :param heatmap_gs:
    """
    # Create heatmap axis
    heatmap_axes = fig.add_subplot(heatmap_gs[1, 1])
    heatmap_axes.set_xticks(np.linspace(0, dfr.shape[0] - 1, dfr.shape[0]))
    heatmap_axes.set_yticks(np.linspace(0, dfr.shape[0] - 1, dfr.shape[0]))
    heatmap_axes.grid(False)
    heatmap_axes.xaxis.tick_bottom()
    heatmap_axes.yaxis.tick_right()
    return heatmap_axes


def add_colorbar(dfr, fig, dend, params, orientation="row"):
    """Add class colorbars to Matplotlib heatmap.

    :param dfr:
    :param fig:
    :param dent:
    :param params:
    :param orientation:
    """
    # Assign a numerical value to each class, for mpl
    classdict = {cls: idx for (idx, cls) in enumerate(params.classes.values())}

    # colourbar
    cblist = []
    for name in [str(_) for _ in dfr.index[dend["dendrogram"]["leaves"]]]:
        if name in params.classes:
            cblist.append(classdict[params.classes[name]])
        elif name in classdict:
            cblist.append(classdict[name])
        else:  # Catches genomes with no assigned class
            cblist.append(0)
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
def add_labels(heatmap_axes, rowlabels, collabels, params):
    """Add labels to Matplotlib heatmap axes, in-place.

    :param heatmap_axes:
    :param rowlabels:
    :param collabels:
    :param params:
    """
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
def add_colorscale(fig, heatmap_gs, ax_map, params, title=None):
    """Add colour scale to heatmap.

    :param fig:
    :param heatmap_gs:
    :param ax_map:
    :param params:
    :param title:
    """
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
def heatmap(dfr, outfilename=None, title=None, params=None):
    """Return matplotlib heatmap with cluster dendrograms.

    :param dfr:  pandas DataFrame with relevant data
    :param outfilename:  path to output file (indicates output format)
    :param params:  a list of parameters for plotting: [colormap, vmin, vmax]
    :param labels:  dictionary of alternative labels, keyed by default sequence labels
    :param classes:  dictionary of sequence classes, keyed by default sequence labels
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
    coldend = add_dendrogram(dfr, fig, params, heatmap_gs, orientation="col")
    rowdend = add_dendrogram(dfr, fig, params, heatmap_gs, orientation="row")

    # Add heatmap axes to figure, with rows/columns as in the dendrograms
    heatmap_axes = get_heatmap_axes(dfr, fig, heatmap_gs)
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
        add_colorbar(dfr, fig, coldend, params, orientation="col")
        add_colorbar(dfr, fig, rowdend, params, orientation="row")

    # Add heatmap labels
    add_labels(
        heatmap_axes, rowdend["dendrogram"]["ivl"], coldend["dendrogram"]["ivl"], params
    )

    # Add colour scale
    add_colorscale(fig, heatmap_gs, ax_map, params, title)

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


def scatter(
    dfr1,
    dfr2,
    outfilename=None,
    matname1="identity",
    matname2="coverage",
    title=None,
    params=None,
):
    """Return matplotlib scatterplot.

    :param dfr1:  pandas DataFrame with x-axis data
    :param dfr2:  pandas DataFrame with y-axis data
    :param outfilename:  path to output file (indicates output format)
    :param matname1:  name of x-axis data
    :param matname2:  name of y-axis data
    :param title:  title for the plot
    :param params:  a list of parameters for plotting: [colormap, vmin, vmax]
    """
    # Make an empty dataframe to collect the input data in
    combined = pd.DataFrame()

    # Add data
    combined[matname1] = dfr1.values.flatten()
    combined[matname2] = dfr2.values.flatten()

    # Add lable information, if available
    # if params.labels:
    #     hue = "labels"
    #  combined['labels'] =   #  add labels to dataframe; unsure of their configuration at this point
    # else:
    hue = None

    fig, ax = plt.subplots(figsize=(8, 8))
    fig.suptitle(title)
    ax.set_xlabel(f"{matname1.title()}")
    ax.set_ylabel(f"{matname2.title()}")

    plt.scatter(matname1, matname2, data=combined, c=hue, s=2)

    # Return figure output, and write, if required
    plt.subplots_adjust(top=0.85)  # Leave room for title
    fig.set_tight_layout(True)
    if outfilename:
        fig.savefig(outfilename)
    return fig


def bland_altman(
    dfr1, dfr2, outfilename, matname1, matname2, title=None, info=None, params=None
):
    """Return matplotlib Bland-Altman plot.

    :param dfr1:  pandas DataFrame with x-axis data
    :param dfr2:  pandas DataFrame with y-axis data
    :param outfilename:  path to output file (indicates output format)
    :param matname1:  name of x-axis data
    :param matname2:  name of y-axis data
    :param title:  title for the plot
    :param info:   information about the data in the plot
    :param params:  a list of parameters for plotting: [colormap, vmin, vmax]
    """
    # Make an empty dataframe to collect the input data in
    data = pd.DataFrame()

    # Add data
    data["avg"] = (dfr1 + dfr2).values.flatten() / 2
    data["AminusB"] = (dfr1 - dfr2).values.flatten()

    # Add lable information, if available
    # if params.labels:
    #     hue = "labels"
    #  combined['labels'] =   #  add labels to dataframe; unsure of their configuration at this point
    # else:
    hue = None

    fig, ax = plt.subplots(figsize=(8, 8))
    fig.suptitle(f"Bland-Altman plot for {matname1}")
    ax.set_xlabel(f"Average of run {matname1} scores")
    ax.set_ylabel(f"Difference between run {matname1} scores")

    plt.scatter("avg", "AminusB", data=data, c=hue, s=2)

    # Return figure output, and write, if required
    plt.subplots_adjust(top=0.85)  # Leave room for title
    fig.set_tight_layout(True)
    if outfilename:
        fig.savefig(outfilename)
    return fig
