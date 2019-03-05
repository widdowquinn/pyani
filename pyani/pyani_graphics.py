# Copyright 2013-2019, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to implement graphics output for ANI analyses."""

# Force matplotlib NOT to use an Xwindows backend on *nix, so that
# _tkinter.TclError is avoided when there is no $DISPLAY env: this can occur
# when running the package/script via ssh
# See http://stackoverflow.com/questions/2801882/\
#            generating-a-png-with-matplotlib-when-display-is-undefined
# This needs to be done before importing pyplot

from math import floor, log10
import warnings

import matplotlib

# Specify matplotlib backend
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np

import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as distance

import seaborn as sns
import pandas as pd

from . import pyani_config


# Register Matplotlib colourmaps
plt.register_cmap(cmap=pyani_config.CMAP_SPBND_BURD)
plt.register_cmap(cmap=pyani_config.CMAP_HADAMARD_BURD)
plt.register_cmap(cmap=pyani_config.CMAP_BURD)


# Convenience class to hold heatmap graphics parameters
class Params(object):  # pylint: disable=too-few-public-methods
    """Convenience class to hold heatmap rendering parameters."""

    def __init__(self, params, labels=None, classes=None):
        self.cmap = plt.get_cmap(params[0])
        self.vmin = params[1]
        self.vmax = params[2]
        self.labels = labels
        self.classes = classes

    @property
    def vdiff(self):
        """Returns difference between max and min values for presentation"""
        return max(0.01, self.vmax - self.vmin)


# helper for cleaning up matplotlib axes by removing ticks etc.
def clean_axis(axis):
    """Remove ticks, tick labels, and frame from axis"""
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
    col_cb = pd.Series(dfr.index).map(lvl_pal)
    # The col_cb Series index now has to match the dfr.index, but
    # we don't create the Series with this (and if we try, it
    # fails) - so change it with this line
    col_cb.index = dfr.index
    return col_cb


# Get safe Seaborn labels
def get_safe_seaborn_labels(dfr, labels):
    """Returns labels guaranteed to correspond to the dataframe."""
    if labels is not None:
        return [labels.get(i, i) for i in dfr.index]
    return [i for i in dfr.index]


# Return a clustermap
def get_seaborn_clustermap(dfr, params, title=None, annot=True):
    """Returns a Seaborn clustermap."""
    fig = sns.clustermap(
        dfr,
        cmap=params.cmap,
        vmin=params.vmin,
        vmax=params.vmax,
        col_colors=params.colorbar,
        row_colors=params.colorbar,
        figsize=(params.figsize, params.figsize),
        linewidths=params.linewidths,
        xticklabels=params.labels,
        yticklabels=params.labels,
        annot=annot,
    )
    fig.cax.yaxis.set_label_position("left")
    if title:
        fig.cax.set_ylabel(title)

    # Rotate ticklabels
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xticklabels(), rotation=90)
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_yticklabels(), rotation=0)

    # Return clustermap
    return fig


# Generate Seaborn heatmap output
def heatmap_seaborn(dfr, outfilename=None, title=None, params=None):
    """Returns seaborn heatmap with cluster dendrograms.

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

    # Labels are defined before we build the clustering
    # If a label mapping is missing, use the key text as fall back
    params.labels = get_safe_seaborn_labels(dfr, params.labels)

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


# Add dendrogram and axes to passed figure
def add_mpl_dendrogram(dfr, fig, heatmap_gs, orientation="col"):
    """Return a dendrogram and corresponding gridspec, attached to the fig

    Modifies the fig in-place. Orientation is either 'row' or 'col' and
    determines location and orientation of the rendered dendrogram.
    """
    # Row or column axes?
    if orientation == "row":
        dists = distance.squareform(distance.pdist(dfr))
        spec = heatmap_gs[1, 0]
        orient = "left"
        nrows, ncols = 1, 2
        height_ratios = [1]
    else:  # Column dendrogram
        dists = distance.squareform(distance.pdist(dfr.T))
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
    )
    clean_axis(dend_axes)
    return {"dendrogram": dend, "gridspec": gspec}


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
    for name in dfr.index[dend["dendrogram"]["leaves"]]:
        if name not in params.classes:
            params.classes[name] = name

    # Assign a numerical value to each class, for mpl
    classdict = {cls: idx for (idx, cls) in enumerate(params.classes.values())}

    # colourbar
    cblist = []
    for name in dfr.index[dend["dendrogram"]["leaves"]]:
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
    """Returns matplotlib heatmap with cluster dendrograms.

    - dfr - pandas DataFrame with relevant data
    - outfilename - path to output file (indicates output format)
    - params - a list of parameters for plotting: [colormap, vmin, vmax]
    - labels - dictionary of alternative labels, keyed by default sequence
               labels
    - classes - dictionary of sequence classes, keyed by default sequence
                labels
    """
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
    coldend = add_mpl_dendrogram(dfr, fig, heatmap_gs, orientation="col")
    rowdend = add_mpl_dendrogram(dfr, fig, heatmap_gs, orientation="row")

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
        heatmap_axes,
        dfr.index[rowdend["dendrogram"]["leaves"]],
        dfr.index[coldend["dendrogram"]["leaves"]],
        params,
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
