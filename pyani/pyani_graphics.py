# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence. 
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to implement graphics output for ANI analyses."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as distance

from matplotlib.colors import LinearSegmentedColormap
from math import floor, log10

# Define custom matplotlib colourmaps
# 1) Map for species boundaries (95%: 0.95), blue for values at
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
    for sp in ax.spines.values():
        sp.set_visible(False)


def heatmap_mpl(df, outfilename, title=None, cmap=None,
                vmin=None, vmax=None):
    """Write matplotlib heatmap with cluster dendrograms.

    - df - pandas DataFrame with relevant data
    - outfilename - path to output file (indicates output format)
    - cmap - colourmap option
    """
    # Get indication of dataframe size and, if necessary, max and 
    # min values for colormap
    size = df.shape[0]  
    if vmin == None:
        vmin = df.values.min()
    if vmax == None:
        vmax = df.values.max()
    vdiff = vmax - vmin
    cbticks = [vmin + e * vdiff for e in (0, 0.25, 0.5, 0.75, 1)]
    if vmax > 10:
        exponent = int(floor(log10(vmax)))
        cbticks = [int(round(e, -exponent)) for e in cbticks]

    # Obtain appropriate colour map
    cmap = plt.get_cmap(cmap)

    # Layout figure grid and add title
    fig = plt.figure()
    if title:
        fig.suptitle(title)
    heatmapGS = gridspec.GridSpec(2, 2, wspace=0.0, hspace=0.0,
                                  width_ratios=[0.3, 1],
                                  height_ratios=[0.3, 1])

    # Calculate column pairwise distances and dendrogram
    coldists = distance.squareform(distance.pdist(df.T))
    colclusters = sch.linkage(coldists, method='complete')

    # Create column dendrogram axis
    coldend_axes = fig.add_subplot(heatmapGS[0,1])
    coldend = sch.dendrogram(colclusters, color_threshold=np.inf)
    clean_axis(coldend_axes)

    # Calculate row pairwise distances and dendrogram
    rowdists = distance.squareform(distance.pdist(df))
    rowclusters = sch.linkage(rowdists, method='complete')

    # Create column dendrogram axis
    rowdend_axes = fig.add_subplot(heatmapGS[1,0])
    rowdend = sch.dendrogram(rowclusters, color_threshold=np.inf,
                             orientation="right")
    clean_axis(rowdend_axes)
        
    # Create heatmap axis
    heatmap_axes = fig.add_subplot(heatmapGS[1,1])
    ax_map = heatmap_axes.imshow(df.ix[rowdend['leaves'],
                                       coldend['leaves']],
                                 interpolation='nearest',
                                 cmap=cmap, origin='lower',
                                 vmin=vmin, vmax=vmax,
                                 aspect='auto')
    heatmap_axes.set_xticks(np.linspace(0, size-1, size))
    heatmap_axes.set_yticks(np.linspace(0, size-1, size))
    heatmap_axes.set_xticklabels(df.index[coldend['leaves']])
    heatmap_axes.set_yticklabels(df.index[rowdend['leaves']])
    heatmap_axes.grid('off')
    heatmap_axes.xaxis.tick_bottom()
    heatmap_axes.yaxis.tick_right()

    # Add colour scale
    scale_subplot =\
                    gridspec.GridSpecFromSubplotSpec(1, 2, 
                        subplot_spec=heatmapGS[0,0], wspace=0.0, hspace=0.0)
    scale_ax = fig.add_subplot(scale_subplot[0,1])
    cb = fig.colorbar(ax_map, scale_ax, ticks=cbticks, )
    cb.set_label('Scale')
    cb.ax.yaxis.set_ticks_position('left')
    cb.ax.yaxis.set_label_position('left')
    cb.outline.set_linewidth(0)

    # Write resulting figure
    fig.set_tight_layout(True)
    fig.savefig(outfilename)
