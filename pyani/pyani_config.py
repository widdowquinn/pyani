# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Configuration settings for the pyani package.
"""

from matplotlib.colors import LinearSegmentedColormap

# Defaults assume that common binaries are on the $PATH
NUCMER_DEFAULT = "nucmer"
FILTER_DEFAULT = "delta-filter"
BLASTN_DEFAULT = "blastn"
MAKEBLASTDB_DEFAULT = "makeblastdb"
BLASTALL_DEFAULT = "blastall"
FORMATDB_DEFAULT = "formatdb"
QSUB_DEFAULT = "qsub"

# Stems for output files
ANIM_FILESTEMS = ("ANIm_alignment_lengths", "ANIm_percentage_identity",
                  "ANIm_alignment_coverage", "ANIm_similarity_errors",
                  "ANIm_hadamard")
ANIB_FILESTEMS = ("ANIb_alignment_lengths", "ANIb_percentage_identity",
                  "ANIb_alignment_coverage", "ANIb_similarity_errors",
                  "ANIb_hadamard")
TETRA_FILESTEMS = ("TETRA_correlations",)
ANIBLASTALL_FILESTEMS = ("ANIblastall_alignment_lengths",
                         "ANIblastall_percentage_identity",
                         "ANIblastall_alignment_coverage",
                         "ANIblastall_similarity_errors",
                         "ANIblastall_hadamard")

# Output subdirectory names for each method
ALIGNDIR = {'ANIm': 'nucmer_output',
            'ANIb': 'blastn_output',
            'ANIblastall': 'blastall_output'}

# Any valid matplotlib colour map can be used here
# See, e.g. http://matplotlib.org/xkcd/examples/color/colormaps_reference.html
MPL_CBAR = 'Spectral'

# Parameters for analyses
FRAGSIZE = 1020  # Default ANIb fragment size

# SGE/OGE scheduler parameters
SGE_WAIT = 0.01  # Base unit of time (s) to wait between polling SGE

# Custom Matplotlib colourmaps
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
CMAP_SPBND_BURD = LinearSegmentedColormap("spbnd_BuRd",
                                          cdict_spbnd_BuRd)

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
CMAP_HADAMARD_BURD = LinearSegmentedColormap("hadamard_BuRd",
                                             cdict_hadamard_BuRd)

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
CMAP_BURD = LinearSegmentedColormap("BuRd", cdict_BuRd)


# Graphics parameters for each output file. Note that this should be
# in sync with the output file stems above
def params_mpl(df):
    """Returns dict of matplotlib parameters, dependent on dataframe."""
    return {'ANIb_alignment_lengths': ('afmhot', df.values.min(),
                                       df.values.max()),
            'ANIb_percentage_identity': ('spbnd_BuRd', 0, 1),
            'ANIb_alignment_coverage': ('BuRd', 0, 1),
            'ANIb_hadamard': ('hadamard_BuRd', 0, 1),
            'ANIb_similarity_errors': ('afmhot', df.values.min(),
                                       df.values.max()),
            'ANIm_alignment_lengths': ('afmhot', df.values.min(),
                                       df.values.max()),
            'ANIm_percentage_identity': ('spbnd_BuRd', 0, 1),
            'ANIm_alignment_coverage': ('BuRd', 0, 1),
            'ANIm_hadamard': ('hadamard_BuRd', 0, 1),
            'ANIm_similarity_errors': ('afmhot', df.values.min(),
                                       df.values.max()),
            'TETRA_correlations': ('spbnd_BuRd', 0, 1),
            'ANIblastall_alignment_lengths': ('afmhot', df.values.min(),
                                              df.values.max()),
            'ANIblastall_percentage_identity': ('spbnd_BuRd', 0, 1),
            'ANIblastall_alignment_coverage': ('BuRd', 0, 1),
            'ANIblastall_hadamard': ('hadamard_BuRd', 0, 1),
            'ANIblastall_similarity_errors': ('afmhot', df.values.min(),
                                              df.values.max())}
