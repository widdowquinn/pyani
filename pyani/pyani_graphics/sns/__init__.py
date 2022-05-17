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
"""Code to implement Seaborn graphics output for ANI analyses."""

import warnings

import matplotlib  # pylint: disable=C0411
import pandas as pd
import seaborn as sns
import logging
from typing import List, Dict
from pyani import pyani_config

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402,E501 # pylint: disable=wrong-import-position,wrong-import-order,ungrouped-imports


# Add classes colorbar to Seaborn plot
def get_colorbar(dfr, classes):
    """Return a colorbar representing classes, for a Seaborn plot.

    :param dfr:
    :param classes:

    The aim is to get a pd.Series for the passed dataframe columns,
    in the form:
    0    colour for class in col 0
    1    colour for class in col 1
    ...  colour for class in col ...
    n    colour for class in col n
    """
    levels = sorted(list(set(classes.values())))
    paldict = dict(
        zip(
            levels,
            sns.cubehelix_palette(
                len(levels), light=0.9, dark=0.1, reverse=True, start=1, rot=-2
            ),
        )
    )
    lvl_pal = {cls: paldict[lvl] for (cls, lvl) in list(classes.items())}
    # Have to use string conversion of the dataframe index, here
    col_cb = pd.Series([str(_) for _ in dfr.index]).map(lvl_pal)
    # The col_cb Series index now has to match the dfr.index, but
    # we don't create the Series with this (and if we try, it
    # fails) - so change it with this line
    col_cb.index = dfr.index
    return col_cb


# Add labels to the seaborn heatmap axes
def add_labels(fig, params):
    """Add labels to Seaborn heatmap axes, in-place.

    :param fig:
    :param params:
    """
    if params.labels:
        # If a label mapping is missing, use the key text as fall back
        for _ in fig.ax_heatmap.get_yticklabels():
            _.set_text(params.labels.get(_.get_text(), _.get_text()))
        for _ in fig.ax_heatmap.get_xticklabels():
            _.set_text(params.labels.get(_.get_text(), _.get_text()))
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xticklabels(), rotation=90)
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_yticklabels(), rotation=0)
    return fig


# Return a clustermap
def get_clustermap(dfr, params, title=None, annot=True):
    """Return a Seaborn clustermap for the passed dataframe.

    :param dfr:
    :param params:
    :param title:  str, plot title
    :param annot:  Boolean, add text for cell values?
    """

    # If we do not catch warnings here, then we often get the following warning:
    #   ClusterWarning: scipy.cluster: The symmetric non-negative hollow
    #   observation matrix looks suspiciously like an uncondensed distance matrix
    # The usual solution would be to convert the array with
    # scipy.spatial.distance.squareform(), but this requires that all values in
    # the main diagonal are zero, which is not the case for ANI.
    # As we know this is a (1-distance) matrix, we could just set the diagonal
    # to zero and fudge it, but this is not a good solution. Instead, we suppress
    # the warning in a context manager for this function call only, because we
    # know the warning is not relevant.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=(
                "scipy.cluster: The symmetric non-negative "
                "hollow observation matrix looks suspiciously like an "
                "uncondensed distance matrix"
            ),
        )
        fig = sns.clustermap(
            dfr,
            cmap=params.cmap,
            vmin=params.vmin,
            vmax=params.vmax,
            center=params.center,
            col_colors=params.colorbar,
            row_colors=params.colorbar,
            figsize=(params.figsize, params.figsize),
            linewidths=params.linewidths,
            annot=annot,
        )

    # add labels for each of the input genomes
    add_labels(fig, params)

    fig.cax.yaxis.set_label_position("left")
    if title:
        fig.ax_heatmap.set_title(title, pad=1000, fontdict={"fontsize": 75})
    # Return clustermap
    return fig


# Generate Seaborn heatmap output
def heatmap(dfr, outfilename=None, title=None, params=None):
    """Return seaborn heatmap with cluster dendrograms.

    :param dfr:  pandas DataFrame with relevant data
    :param outfilename:  path to output file (indicates output format)
    :param title:
    :param params:
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
        col_cb = get_colorbar(dfr, params.classes)

    # Add attributes to parameter object, and draw heatmap
    params.colorbar = col_cb
    params.figsize = figsize
    params.linewidths = 0.25
    fig = get_clustermap(dfr, params, title=title)

    # Save to file
    if outfilename:
        fig.savefig(outfilename)

    # Return clustermap
    return fig


def distribution(dfr, outfilename, matname, title=None):
    """Return seaborn distribution plot for matrix.

    :param drf:  DataFrame with results matrix
    :param outfilename:  Path to output file for writing
    :param matname:  str, type of matrix being plotted
    :param title:  str, optional title
    """
    fill = "#A6C8E0"
    rug = "#2678B2"
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(title)
    sns.histplot(
        dfr.values.flatten(),
        ax=axes[0],
        stat="count",
        element="step",
        color=fill,
        edgecolor=fill,
    )
    axes[0].set_ylim(ymin=0)
    sns.kdeplot(dfr.values.flatten(), ax=axes[1])
    sns.rugplot(dfr.values.flatten(), ax=axes[1], color=rug)

    # Modify axes after data is plotted
    for _ in axes:
        if matname.endswith("absdiffs"):
            _.set_xlim(0, _.get_xlim()[1])
        elif matname.endswith("diffs"):
            pass
        elif matname.split("_")[0] == "sim_errors":
            _.set_xlim(0, _.get_xlim()[1])
        elif matname.split("_")[0] in ["hadamard", "coverage"]:
            _.set_xlim(0, 1.01)
        elif matname.split("_")[0] == "identity":
            _.set_xlim(0.75, 1.01)

    # Tidy figure
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if outfilename:
        # For some reason seaborn gives us an AxesSubPlot with
        # sns.distplot, rather than a Figure, so we need this hack
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
    """Return seaborn scatterplot.

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

    # Create the plot
    fig = sns.lmplot(
        x=matname1,
        y=matname2,
        data=combined,
        hue=hue,
        fit_reg=False,
        scatter_kws={"s": 2},
    )
    fig.set(xlabel=matname1.title(), ylabel=matname2.title())
    plt.title(title)

    # Save to file
    if outfilename:
        fig.savefig(outfilename)

    # Return clustermap
    return fig


def bland_altman(
    dfr1,
    dfr2,
    outfilename,
    matname1,
    matname2,
    run_ids,
    title=None,
    info=None,
    params=None,
):
    """Return seaborn Bland-Altman plot.

    :param dfr1:  pandas DataFrame with x-axis data
    :param dfr2:  pandas DataFrame with y-axis data
    :param outfilename:  path to output file (indicates output format)
    :param matname1:  name of x-axis data
    :param matname2:  name of y-axis data
    :param run_ids:   tuple of run_ids (ref, query)
    :param title:  title for the plot
    :param info:   information about the data in the plot
    :param params:  a list of parameters for plotting: [colormap, vmin, vmax]
    """
    data = pd.DataFrame()
    ref_id, query_id = run_ids
    data["avg"] = (dfr1 + dfr2).values.flatten() / 2
    data["AminusB"] = (dfr1 - dfr2).values.flatten()

    fig = sns.lmplot(
        x="avg", y="AminusB", data=data, fit_reg=False, scatter_kws={"s": 2}, height=9
    )

    fig.ax.hlines(0, fig.ax.get_xbound()[0], fig.ax.get_xbound()[1], linewidths=1)
    fig.ax.margins(x=0)
    # fig.figtext(1, .5, info)
    fig.set(
        xlabel=f"Average of run {matname1} scores",
        ylabel=f"Difference between {matname1} scores (run {ref_id} - run {query_id})",
    )
    plt.title(title)

    fig.tight_layout()

    if outfilename:
        plt.savefig(outfilename)

    return fig
