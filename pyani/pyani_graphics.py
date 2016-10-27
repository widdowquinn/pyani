# Copyright 2013-2016, The James Hutton Insitute
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

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
    
import numpy as np
import warnings

import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as distance

import seaborn as sns
import pandas as pd

from . import pyani_config

from math import floor, log10

# Define custom matplotlib colourmaps
# 1a) Map for species boundaries (95%: 0.95), blue for values at
# 0.9 or below, red for values at 1.0; white at 0.95.
# Also, anything below 0.7 is 70% grey
cdict_spbnd_BuRd = {'red': ((0.0, 0.0, 0.7),
                            (0.7, 0.7, 0.0),
                            (0.9, 0.0, 0.0),
                            (0.95, 1.0, 1.0),
                            (1.0, 1.0, 1.0)),
                    'green': ((0.0, 0.0, 0.7),
                              (0.7, 0.7, 0.0),
                              (0.9, 0.0, 0.0),
                              (0.95, 1.0, 1.0),
                              (1.0, 0.0, 0.0)),
                    'blue': ((0.0, 0.0, 0.7),
                             (0.7, 0.7, 1.0),
                             (0.95, 1.0, 1.0),
                             (1.0, 0.0, 0.0))}
cmap_spbnd_BuRd = LinearSegmentedColormap("spbnd_BuRd", cdict_spbnd_BuRd)
plt.register_cmap(cmap=cmap_spbnd_BuRd)

# 1b) Map for species boundaries (95%: 0.95), blue for values at
# 0.9 or below, red for values at 1.0; white at 0.9.
# Also, anything below 0.8 is 70% grey
cdict_hadamard_BuRd = {'red': ((0.0, 0.0, 0.7),
                            (0.8, 0.7, 0.0),
                            (0.9, 0.0, 0.0),
                            (0.9, 1.0, 1.0),
                            (1.0, 1.0, 1.0)),
                    'green': ((0.0, 0.0, 0.7),
                              (0.8, 0.7, 0.0),
                              (0.9, 0.0, 0.0),
                              (0.9, 1.0, 1.0),
                              (1.0, 0.0, 0.0)),
                    'blue': ((0.0, 0.0, 0.7),
                             (0.8, 0.7, 1.0),
                             (0.9, 1.0, 1.0),
                             (1.0, 0.0, 0.0))}
cmap_hadamard_BuRd = LinearSegmentedColormap("hadamard_BuRd",
                                             cdict_hadamard_BuRd)
plt.register_cmap(cmap=cmap_hadamard_BuRd)

# 2) Blue for values at 0.0, red for values at 1.0; white at 0.5
cdict_BuRd = {'red': ((0.0, 0.0, 0.0),
                      (0.5, 1.0, 1.0),
                      (1.0, 1.0, 1.0)),
              'green': ((0.0, 0.0, 0.0),
                        (0.5, 1.0, 1.0),
                        (1.0, 0.0, 0.0)),
              'blue': ((0.0, 1.0, 1.0),
                       (0.5, 1.0, 1.0),
                       (1.0, 0.0, 0.0))}
cmap_BuRd = LinearSegmentedColormap("BuRd", cdict_BuRd)
plt.register_cmap(cmap=cmap_BuRd)


# helper for cleaning up matplotlib axes by removing ticks etc.
def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in list(ax.spines.values()):
        sp.set_visible(False)


# Generate Seaborn heatmap output
def heatmap_seaborn(df, outfilename=None, title=None, cmap=None,
                    vmin=None, vmax=None, labels=None, classes=None):
    """Returns seaborn heatmap with cluster dendrograms.

    - df - pandas DataFrame with relevant data
    - outfilename - path to output file (indicates output format)
    - cmap - colourmap option
    - vmin - float, minimum value on the heatmap scale
    - vmax - float, maximum value on the heatmap scale
    - labels - dictionary of alternative labels, keyed by default sequence
               labels
    - classes - dictionary of sequence classes, keyed by default sequence
                labels
    """
    # Obtain colour map
    cmap = plt.get_cmap(cmap)

    # Decide on figure layout size: a minimum size is required for
    # aesthetics, and a maximum to avoid core dumps on rendering.
    # If we hit the maximum size, we should modify font size.
    maxfigsize = 120
    calcfigsize = df.shape[0] * 1.1
    figsize = min(max(8, calcfigsize), maxfigsize)
    if figsize == maxfigsize:
        scale = maxfigsize/calcfigsize
        sns.set_context("notebook", font_scale=scale)

    # Add class colour bar. The aim is to get a pd.Series for the columns
    # of the form:
    # 0    colour for class in col 0
    # 1    colour for class in col 1
    # ...  colour for class in col ...
    # n    colour for class in col n
    # This is in col_cb when we're finished
    if classes is not None:
        levels = sorted(list(set(classes.values())))
        pal = sns.cubehelix_palette(len(levels),
                                    light=.9, dark=.1, reverse=True,
                                    start=1, rot=-2)
        paldict = {lvl: pal for (lvl, pal) in zip(levels, pal)}
        lvl_pal = {cls: paldict[lvl] for (cls, lvl) in list(classes.items())}
        col_cb = pd.Series(df.index).map(lvl_pal)
        # The col_cb Series index now has to match the df.index, but
        # we don't create the Series with this (and if we try, it
        # fails) - so change it, here.
        col_cb.index = df.index
    else:
        col_cb = None

    # Labels are defined before we build the clustering
    # If a label mapping is missing, use the key text as fall back
    if labels is not None:
        newlabels = [labels.get(i, i) for i in df.index]
    else:
        newlabels = [i for i in df.index]

    # Plot heatmap
    fig = sns.clustermap(df,
                         cmap=cmap, vmin=vmin, vmax=vmax,
                         col_colors=col_cb, row_colors=col_cb,
                         figsize=(figsize, figsize),
                         linewidths=0.25,
                         xticklabels=newlabels,
                         yticklabels=newlabels,
                         annot=True)
    fig.cax.yaxis.set_label_position('left')
    fig.cax.set_ylabel(title)

    # Rotate ticklabels
    fig.ax_heatmap.set_xticklabels(fig.ax_heatmap.get_xticklabels(),
                                   rotation=90)
    fig.ax_heatmap.set_yticklabels(fig.ax_heatmap.get_yticklabels(),
                                   rotation=0)

    # Save to file
    if outfilename:
        fig.savefig(outfilename)
    return fig


# Generate Matplotlib heatmap output
def heatmap_mpl(df, outfilename=None, title=None, cmap=None,
                vmin=None, vmax=None, labels=None, classes=None):
    """Returns matplotlib heatmap with cluster dendrograms.

    - df - pandas DataFrame with relevant data
    - outfilename - path to output file (indicates output format)
    - cmap - colourmap option
    - vmin - float, minimum value on the heatmap scale
    - vmax - float, maximum value on the heatmap scale
    - labels - dictionary of alternative labels, keyed by default sequence
               labels
    - classes - dictionary of sequence classes, keyed by default sequence
                labels
    """
    # Get indication of dataframe size and, if necessary, max and
    # min values for colormap
    dfsize = df.shape[0]
    if vmin is None:
        vmin = df.values.min()
    if vmax is None:
        vmax = df.values.max()
    # a vdiff of zero causes display problems in matplotlib, so we set a minval
    vdiff = max(vmax - vmin, 0.01)
    cbticks = [vmin + e * vdiff for e in (0, 0.25, 0.5, 0.75, 1)]
    if vmax > 10:
        exponent = int(floor(log10(vmax))) - 1
        cbticks = [int(round(e, -exponent)) for e in cbticks]

    # Obtain appropriate colour map
    cmap = plt.get_cmap(cmap)

    # Layout figure grid and add title
    # Set figure size by the number of rows in the dataframe
    figsize = max(8, dfsize * 0.175)
    fig = plt.figure(figsize=(figsize, figsize))
    # if title:
    #     fig.suptitle(title)
    heatmapGS = gridspec.GridSpec(2, 2, wspace=0.0, hspace=0.0,
                                  width_ratios=[0.3, 1],
                                  height_ratios=[0.3, 1])

    # Calculate column pairwise distances and dendrogram
    coldists = distance.squareform(distance.pdist(df.T))
    colclusters = sch.linkage(coldists, method='complete')

    # Create column dendrogram axis
    colGS = gridspec.GridSpecFromSubplotSpec(2, 1,
                                             subplot_spec=heatmapGS[0, 1],
                                             wspace=0.0, hspace=0.1,
                                             height_ratios=[1, 0.15])
    coldend_axes = fig.add_subplot(colGS[0, 0])
    coldend = sch.dendrogram(colclusters, color_threshold=np.inf)
    clean_axis(coldend_axes)

    # Calculate row pairwise distances and dendrogram
    rowdists = distance.squareform(distance.pdist(df))
    rowclusters = sch.linkage(rowdists, method='complete')

    # Create row dendrogram axis
    rowGS = gridspec.GridSpecFromSubplotSpec(1, 2,
                                             subplot_spec=heatmapGS[1, 0],
                                             wspace=0.1, hspace=0.0,
                                             width_ratios=[1, 0.15])
    rowdend_axes = fig.add_subplot(rowGS[0, 0])
    rowdend = sch.dendrogram(rowclusters, color_threshold=np.inf,
                             orientation="left")
    clean_axis(rowdend_axes)

    # Create heatmap axis
    heatmap_axes = fig.add_subplot(heatmapGS[1, 1])
    ax_map = heatmap_axes.imshow(df.ix[rowdend['leaves'],
                                       coldend['leaves']],
                                 interpolation='nearest',
                                 cmap=cmap, origin='lower',
                                 vmin=vmin, vmax=vmax,
                                 aspect='auto')
    heatmap_axes.set_xticks(np.linspace(0, dfsize-1, dfsize))
    heatmap_axes.set_yticks(np.linspace(0, dfsize-1, dfsize))
    heatmap_axes.grid('off')
    heatmap_axes.xaxis.tick_bottom()
    heatmap_axes.yaxis.tick_right()

    # Are there class colourbars to add?
    if classes is not None:
        # If not all dendrogram leaves are present, extend classes
        for name in df.index[coldend['leaves']]:
            if name not in classes:
                classes[name] = name

        # Assign a numerical value to each class, for mpl
        classdict = {cls: idx for (idx, cls) in enumerate(classes.values())}
        # Column colourbar
        col_cblist = []
        for name in df.index[coldend['leaves']]:
            try:
                col_cblist.append(classdict[classes[name]])
            except KeyError:
                col_cblist.append(classdict[name])
        col_cb = pd.Series(col_cblist)

        # Row colourbar
        row_cblist = []
        for name in df.index[rowdend['leaves']]:
            try:
                row_cblist.append(classdict[classes[name]])
            except KeyError:
                row_cblist.append(classdict[name])
        row_cb = pd.Series(row_cblist)

        # Create column colourbar axis
        col_cbaxes = fig.add_subplot(colGS[1, 0])
        col_axi = col_cbaxes.imshow([col_cb],
                                    cmap=plt.get_cmap(pyani_config.MPL_CBAR),
                                    interpolation='nearest', aspect='auto',
                                    origin='lower')
        clean_axis(col_cbaxes)
        # Create row colourbar axis
        row_cbaxes = fig.add_subplot(rowGS[0, 1])
        row_axi = row_cbaxes.imshow([[x] for x in row_cb.values],
                                    cmap=plt.get_cmap(pyani_config.MPL_CBAR),
                                    interpolation='nearest', aspect='auto',
                                    origin='lower')
        clean_axis(row_cbaxes)

    # Add heatmap labels
    rowticklabels = df.index[rowdend['leaves']]
    colticklabels = df.index[coldend['leaves']]
    if labels:
        # If a label mapping is missing, use the key text as fall back
        rowticklabels = [labels.get(lab, lab) for lab in rowticklabels]
        colticklabels = [labels.get(lab, lab) for lab in colticklabels]
    xlabs = heatmap_axes.set_xticklabels(colticklabels)
    ylabs = heatmap_axes.set_yticklabels(rowticklabels)
    for label in xlabs:  # Rotate column labels
        label.set_rotation(90)
    for labset in (xlabs, ylabs):  # Smaller font
        for label in labset:
            label.set_fontsize(8)

    # Add colour scale
    scale_subplot =\
        gridspec.GridSpecFromSubplotSpec(1, 3,
                                         subplot_spec=heatmapGS[0, 0],
                                         wspace=0.0, hspace=0.0)
    scale_ax = fig.add_subplot(scale_subplot[0, 1])
    cb = fig.colorbar(ax_map, scale_ax, ticks=cbticks)
    if title:
        cb.set_label(title, fontsize=6)
    cb.ax.yaxis.set_ticks_position('left')
    cb.ax.yaxis.set_label_position('left')
    cb.ax.tick_params(labelsize=6)
    cb.outline.set_linewidth(0)

    # Return figure output, and write, if required
    plt.subplots_adjust(top=0.85)  # Leave room for title
    # fig.set_tight_layout(True)
    # We know that there is a UserWarning here about tight_layout and
    # using the Agg renderer on OSX, so catch and ignore it, for cleanliness.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        heatmapGS.tight_layout(fig, h_pad=0.1, w_pad=0.5)
    if outfilename:
        fig.savefig(outfilename)
    return fig
