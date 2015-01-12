# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence. 
# Please see the LICENSE file that should have been included as part of
# this package.

"""Configuration settings for the pyani package.
"""

# defaults assume that common binaries are on the $PATH
NUCMER_DEFAULT="nucmer"
BLASTN_DEFAULT="blastn"
MAKEBLASTDB_DEFAULT="makeblastdb"

# stems for output files
ANIM_FILESTEMS = ("ANIm_alignment_lengths", "ANIm_percentage_identity",
                  "ANIm_alignment_coverage", "ANIm_similarity_errors")
ANIB_FILESTEMS = ("ANIb_alignment_lengths", "ANIb_percentage_identity",
                  "ANIb_alignment_coverage", "ANIb_similarity_errors")

