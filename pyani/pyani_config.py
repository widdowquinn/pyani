# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Configuration settings for the pyani package.
"""

# Defaults assume that common binaries are on the $PATH
NUCMER_DEFAULT = "nucmer"
BLASTN_DEFAULT = "blastn"
MAKEBLASTDB_DEFAULT = "makeblastdb"
BLASTALL_DEFAULT = "blastall"
FORMATDB_DEFAULT = "formatdb"
QSUB_DEFAULT = "qsub"

# Stems for output files
ANIM_FILESTEMS = ("ANIm_alignment_lengths", "ANIm_percentage_identity",
                  "ANIm_alignment_coverage", "ANIm_similarity_errors")
ANIB_FILESTEMS = ("ANIb_alignment_lengths", "ANIb_percentage_identity",
                  "ANIb_alignment_coverage", "ANIb_similarity_errors")
TETRA_FILESTEMS = ("TETRA_correlations",)
ANIBLASTALL_FILESTEMS = ("ANIblastall_alignment_lengths",
                         "ANIblastall_percentage_identity",
                         "ANIblastall_alignment_coverage",
                         "ANIblastall_similarity_errors")

# Output subdirectory names for each method
ALIGNDIR = {'ANIm': 'nucmer_output',
            'ANIb': 'blastn_output',
            'ANIblastall': 'blastall_output'}

# Colour gradients for use in R and Matplotlib:
# R_AFMHOT is a custom colour palette
R_AFMHOT = 'colorRampPalette(c("black","red","yellow","white"))'
# Any valid matplotlib colour map can be used here
# See, e.g. http://matplotlib.org/xkcd/examples/color/colormaps_reference.html
MPL_CBAR = 'Spectral'

# Parameters for analyses
FRAGSIZE = 1020  # Default ANIb fragment size

# SGE/OGE scheduler parameters
SGE_WAIT = 0.01  # Base unit of time (s) to wait between polling SGE


# Graphics parameters for each output file. Note that this should be
# in sync with the output file stems above
def params_mpl(df):
    """Returns dict of matplotlib parameters, dependent on dataframe."""
    return {'ANIb_alignment_lengths': ('afmhot', df.values.min(),
                                       df.values.max()),
            'ANIb_percentage_identity': ('spbnd_BuRd', 0, 1),
            'ANIb_alignment_coverage': ('BuRd', 0, 1),
            'ANIb_similarity_errors': ('afmhot', df.values.min(),
                                       df.values.max()),
            'ANIm_alignment_lengths': ('afmhot', df.values.min(),
                                       df.values.max()),
            'ANIm_percentage_identity': ('spbnd_BuRd', 0, 1),
            'ANIm_alignment_coverage': ('BuRd', 0, 1),
            'ANIm_similarity_errors': ('afmhot', df.values.min(),
                                       df.values.max()),
            'TETRA_correlations': ('spbnd_BuRd', 0, 1),
            'ANIblastall_alignment_lengths': ('afmhot', df.values.min(),
                                              df.values.max()),
            'ANIblastall_percentage_identity': ('spbnd_BuRd', 0, 1),
            'ANIblastall_alignment_coverage': ('BuRd', 0, 1),
            'ANIblastall_similarity_errors': ('afmhot', df.values.min(),
                                              df.values.max())}


def params_r(df):
    """Returns dictionary of R rendering parameters, dependent on dataframe."""
    return {'ANIb_alignment_lengths': (R_AFMHOT, df.values.min(),
                                       df.values.max()),
            'ANIb_percentage_identity': ('bluered', 0.9, 1),
            'ANIb_alignment_coverage': ('bluered', 0, 1),
            'ANIb_similarity_errors': (R_AFMHOT, df.values.min(),
                                       df.values.max()),
            'ANIm_alignment_lengths': (R_AFMHOT, df.values.min(),
                                       df.values.max()),
            'ANIm_percentage_identity': ('bluered', 0.9, 1),
            'ANIm_alignment_coverage': ('bluered', 0, 1),
            'ANIm_similarity_errors': (R_AFMHOT, df.values.min(),
                                       df.values.max()),
            'TETRA_correlations': ('bluered', 0.9, 1),
            'ANIblastall_alignment_lengths': (R_AFMHOT, df.values.min(),
                                              df.values.max()),
            'ANIblastall_percentage_identity': ('bluered', 0.9, 1),
            'ANIblastall_alignment_coverage': ('bluered', 0, 1),
            'ANIblastall_similarity_errors': (R_AFMHOT, df.values.min(),
                                              df.values.max())}
