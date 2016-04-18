# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to implement the TETRA average nucleotide identity method.

Provides functions for calculation of TETRA as described in:

Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for the
prokaryotic species definition. Proc Natl Acad Sci USA 106: 19126-19131.
doi:10.1073/pnas.0906412106.

and

Teeling et al. (2004) Application of tetranucleotide frequencies for the
assignment of genomic fragments. Env. Microbiol. 6(9): 938-947.
doi:10.1111/j.1462-2920.2004.00624.x
"""

import pandas as pd

import collections
import os
import math

from Bio import SeqIO


# Calculate tetranucleotide Z-score for a set of input sequences
def calculate_tetra_zscores(infilenames):
    """Returns dictionary of TETRA Z-scores for each input file.

    - infilenames - collection of paths to sequence files
    """
    org_tetraz = {}
    for filename in infilenames:
        org = os.path.splitext(os.path.split(filename)[-1])[0]
        org_tetraz[org] = calculate_tetra_zscore(filename)
    return org_tetraz


# Calculate tetranucleotide Z-score for a single sequence file
def calculate_tetra_zscore(filename):
    """Returns TETRA Z-score for the sequence in the passed file.

    - filename - path to sequence file

    Calculates mono-, di-, tri- and tetranucleotide frequencies
    for each sequence, on each strand, and follows Teeling et al. (2004)
    in calculating a corresponding Z-score for each observed
    tetranucleotide frequency, dependent on the mono-, di- and tri-
    nucleotide frequencies for that input sequence.
    """
    # For the Teeling et al. method, the Z-scores require us to count
    # mono, di, tri and tetranucleotide sequences
    monocnt, dicnt, tricnt, tetracnt = (collections.defaultdict(int),
                                        collections.defaultdict(int),
                                        collections.defaultdict(int),
                                        collections.defaultdict(int))
    for rec in SeqIO.parse(filename, 'fasta'):
        for s in [str(rec.seq).upper(),
                  str(rec.seq.reverse_complement()).upper()]:
            # The Teeling et al. algorithm requires us to consider
            # both strand orientations, so monocounts are easy
            monocnt['G'] += s.count('G')
            monocnt['C'] += s.count('C')
            monocnt['T'] += s.count('T')
            monocnt['A'] += s.count('A')
            # For di, tri and tetranucleotide counts, loop over the
            # sequence and its reverse complement, until near the end:
            for i in range(len(s[:-4])):
                di, tri, tetra = s[i:i+2], s[i:i+3], s[i:i+4]
                dicnt[str(di)] += 1
                tricnt[str(tri)] += 1
                tetracnt[str(tetra)] += 1
            # Then clean up the straggling bit at the end:
            tricnt[str(s[-4:-1])] += 1
            tricnt[str(s[-3:])] += 1
            dicnt[str(s[-4:-2])] += 1
            dicnt[str(s[-3:-1])] += 1
            dicnt[str(s[-2:])] += 1
    # Following Teeling (2004), calculate expected frequencies for each
    # tetranucleotide; we ignore ambiguity symbols
    tetra_exp = {}
    for t in [tet for tet in tetracnt if tetra_clean(tet)]:
        tetra_exp[t] = 1.*tricnt[t[:3]]*tricnt[t[1:]]/dicnt[t[1:3]]
    # TODO: Looks like we can collapse the two loops below into one
    # Following Teeling (2004) we approximate the std dev of each
    # tetranucleotide
    tetra_sd = {}
    for t, exp in list(tetra_exp.items()):
        den = dicnt[t[1:3]]
        tetra_sd[t] = math.sqrt(exp * (den - tricnt[t[:3]]) *
                                (den - tricnt[t[1:]]) / (den * den))
    # Following Teeling (2004) we calculate the Z-score for each
    # tetranucleotide
    tetra_z = {}
    for t, exp in list(tetra_exp.items()):
        try:
            tetra_z[t] = (tetracnt[t] - exp)/tetra_sd[t]
        except ZeroDivisionError:
            # We hit a zero in the estimation of variance
            zeroes = [k for k, v in list(tetra_sd.items()) if v == 0]
            tetra_z[t] = 1 / (dicnt[t[1:3]] * dicnt[t[1:3]])
    return tetra_z


# Returns true if the passed string contains only A, C, G or T
def tetra_clean(s):
    """ Checks that a passed string contains only unambiguous IUPAC nucleotide
        symbols. We are assuming that a low frequency of IUPAC ambiguity
        symbols doesn't affect our calculation.
    """
    if not len(set(s) - set('ACGT')):
        return True
    return False


# Calculate Pearson's correlation coefficient from the Z-scores for each
# tetranucleotide. If we're forcing rpy2, might as well use that, though...
def calculate_correlations(tetra_z):
    """Returns dataframe of Pearson correlation coefficients.

    - tetra_z - dictionary of Z-scores, keyed by sequence ID

    Calculates Pearson correlation coefficient from Z scores for each
    tetranucleotide. This is done longhand here, which is fast enough,
    but for robustness we might want to do something else... (TODO).

    Note that we report a correlation by this method, rather than a
    percentage identity.
    """
    orgs = sorted(tetra_z.keys())
    correlations = pd.DataFrame(index=orgs, columns=orgs,
                                dtype=float).fillna(1.0)
    for idx, o1 in enumerate(orgs[:-1]):
        for o2 in orgs[idx+1:]:
            assert sorted(tetra_z[o1].keys()) == sorted(tetra_z[o2].keys())
            tets = sorted(tetra_z[o1].keys())
            z1 = [tetra_z[o1][t] for t in tets]
            z2 = [tetra_z[o2][t] for t in tets]
            z1_mean = sum(z1) / len(z1)
            z2_mean = sum(z2) / len(z2)
            z1diffs = [z - z1_mean for z in z1]
            z2diffs = [z - z2_mean for z in z2]
            diffprods = sum([z1diffs[i] * z2diffs[i] for i in
                             range(len(z1diffs))])
            z1diffs2 = sum([z * z for z in z1diffs])
            z2diffs2 = sum([z * z for z in z2diffs])
            correlations[o1][o2] = diffprods/math.sqrt(z1diffs2 * z2diffs2)
            correlations[o2][o1] = correlations[o1][o2]
    return correlations
