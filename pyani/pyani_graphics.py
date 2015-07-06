# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to implement graphics output for ANI analyses."""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as distance
# import rpy2.robjects as robjects
import pyani_config
import pandas as pd
import warnings

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
    vdiff = vmax - vmin
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
    #if title:
    #    fig.suptitle(title)
    heatmapGS = gridspec.GridSpec(2, 2, wspace=0.0, hspace=0.0,
                                  width_ratios=[0.3, 1],
                                  height_ratios=[0.3, 1])

    # Calculate column pairwise distances and dendrogram
    coldists = distance.squareform(distance.pdist(df.T))
    colclusters = sch.linkage(coldists, method='complete')

    # Create column dendrogram axis
    colGS = gridspec.GridSpecFromSubplotSpec(2, 1,
                                             subplot_spec = heatmapGS[0, 1],
                                             wspace = 0.0, hspace = 0.1,
                                             height_ratios = [1, 0.15])
    coldend_axes = fig.add_subplot(colGS[0, 0])
    coldend = sch.dendrogram(colclusters, color_threshold=np.inf)
    clean_axis(coldend_axes)

    # Calculate row pairwise distances and dendrogram
    rowdists = distance.squareform(distance.pdist(df))
    rowclusters = sch.linkage(rowdists, method='complete')

    # Create row dendrogram axis
    rowGS = gridspec.GridSpecFromSubplotSpec(1, 2,
                                             subplot_spec = heatmapGS[1, 0],
                                             wspace = 0.1, hspace = 0.0,
                                             width_ratios = [1, 0.15])
    rowdend_axes = fig.add_subplot(rowGS[0, 0])
    rowdend = sch.dendrogram(rowclusters, color_threshold=np.inf,
                             orientation="right")
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
        classdict = {cls:idx for (idx, cls) in enumerate(classes.values())}
        # Column colourbar
        col_cb = pd.Series([classdict[classes[name]] for name in
                            df.index[coldend['leaves']]])
        # Row colourbar
        row_cb = pd.Series([classdict[classes[name]] for name in
                            df.index[rowdend['leaves']]])

        # Create column colourbar axis
        col_cbaxes = fig.add_subplot(colGS[1, 0])
        col_axi = col_cbaxes.imshow([col_cb],
                                    cmap=plt.get_cmap(pyani_config.MPL_CBAR),
                                    interpolation = 'nearest', aspect = 'auto',
                                    origin='lower')
        clean_axis(col_cbaxes)
        # Create row colourbar axis
        row_cbaxes = fig.add_subplot(rowGS[0, 1])
        row_axi = row_cbaxes.imshow([[x] for x in row_cb.values],
                                    cmap=plt.get_cmap(pyani_config.MPL_CBAR),
                                    interpolation = 'nearest', aspect = 'auto',
                                    origin='lower')
        clean_axis(row_cbaxes)

    # Add heatmap labels
    rowticklabels = df.index[rowdend['leaves']]
    colticklabels = df.index[coldend['leaves']]
    if labels:
        rowticklabels = [labels[lab] for lab in rowticklabels]
        colticklabels = [labels[lab] for lab in colticklabels]
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
    #fig.set_tight_layout(True)
    # We know that there is a UserWarning here about tight_layout and
    # using the Agg renderer on OSX, so catch and ignore it, for cleanliness.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        heatmapGS.tight_layout(fig, h_pad=0.1, w_pad=0.5)
    if outfilename:
        fig.savefig(outfilename)
    return fig


# Draw heatmap with R
# def heatmap_r(infilename, outfilename, title=None, cmap="bluered",
#               vmin=None, vmax=None, gformat=None, labels=None, classes=None):
#     """Uses R to draw heatmap, and returns R code used for rendering.
#
#     - infilename - path to tab-separated table with data
#     - outfilename - path to output file
#     - cmap - colourmap option
#     - vmin - float, minimum value on the heatmap scale
#     - vmax - float, maximum value on the heatmap scale
#     - gformat - string indicating graphics output format
#     - labels - dictionary of alternative labels, keyed by default sequence
#                labels
#     - classes - dictionary of sequence classes, keyed by default sequence
#                 labels
#     """
#     vdiff = vmax - vmin
#     vstep = 0.001 * vdiff
#
#     # Read the data in to get row/column information for labels and classes,
#     # and so we can guesstimate output image size
#     df = pd.DataFrame.from_csv(infilename, sep="\t", header=0)
#
#     # Prepare R code
#     rstr = ["library(gplots)", "library(RColorBrewer)"]  # R import
#     rstr.append("ani = read.table('%s', header=T, sep='\\t', row.names=1)" %
#                 infilename)
#     rstr.append("%s('%s')" % (gformat, outfilename))
#     rstr.append("par(cex.main=0.75)")
#     if classes:  # Define colour list
#         lbls, cls = [], []
#         for lbl, cla in classes.items():
#             lbls.append(lbl)
#             cls.append(cla)
#         lablist = ','.join(['"%s"' % l for l in lbls])
#         clslist = ','.join(['"%s"' % c for c in cls])
#         rstr.append("labels = data.frame(row.names=c(%s), class=c(%s))" %
#                     (lablist, clslist))
#         rstr.append('lablist = sort(unique(labels[,"class"]))')
#         rstr.append("colourlist = rainbow(length(lablist))")
#         rstr.append("label_colours = data.frame(row.names=sort(lablist), " +
#                     "colours=colourlist)")
#         rstr.append("labels$colours = label_colours[labels$class,]")
#         rstr.append("rowlabels = labels[rownames(ani),]")
#         rstr.append("collabels = labels[colnames(ani),]")
#     cmd = "heatmap.2(as.matrix(ani), col=%s, " % cmap +\
#           "breaks=seq(%.2f, %.2f, %f), " % (vmin, vmax, vstep) +\
#           "trace='none', " +\
#           "margins=c(15, 15), cexCol=0.05 + 1/log10(2 * ncol(ani)), " +\
#           "cexRow=0.05 + 1/log10(2 * nrow(ani)), " +\
#           "main='%s'" % title
#     if labels:  # Change labels from data to passed dict
#         labrow = ','.join(['"%s"' % labels[rowname] for rowname in df.index])
#         labcol = ','.join(['"%s"' % labels[colname] for colname in df.columns])
#         cmd += ", labRow=c(%s), labCol=c(%s)" % (labrow, labcol)
#     if classes:
#         cmd += ', RowSideColors=as.vector(rowlabels[,"colours"])'
#         cmd += ', ColSideColors=as.vector(collabels[,"colours"])'
#     cmd += ")"
#     rstr.append(cmd)
#     if classes:  # Add legend for colour labels
#         rstr.append('lablist = sort(unique(labels[,"class"]))')
#         rstr.append("par(lend=1)")
#         rstr.append('legend("topright", legend=sort(lablist), ' +
#                     "col=colourlist, lty=1, lwd=10)")
#     rstr.append("dev.off()")
#
#     # Execute R code
#     rstr = '\n'.join(rstr)
#     robjects.r(rstr)
#     return rstr + '\n'
