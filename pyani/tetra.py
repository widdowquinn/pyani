# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2024
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2019-2024 University of Strathclyde
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

import collections
import math

from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd  # type: ignore

from Bio import SeqIO  # type: ignore


# Calculate tetranucleotide Z-score for a set of input sequences
def calculate_tetra_zscores(infilenames: Iterable) -> Dict[str, Dict[str, float]]:
    """Return dictionary of TETRA Z-scores for each input file.

    :param infilenames:  iterable of paths to input sequence files
    """
    org_tetraz = {}  # type: Dict[str, Dict[str, float]]
    for filename in infilenames:
        org_tetraz[filename.stem] = calculate_tetra_zscore(filename)
    return org_tetraz


# Calculate tetranucleotide Z-score for a single sequence file
def calculate_tetra_zscore(filename: Path) -> Dict[str, float]:
    """Return TETRA Z-score for the sequence in the passed file.

    :param filename:  path to sequence file

    Calculates mono-, di-, tri- and tetranucleotide frequencies
    for each sequence, on each strand, and follows Teeling et al. (2004)
    in calculating a corresponding Z-score for each observed
    tetranucleotide frequency, dependent on the mono-, di- and tri-
    nucleotide frequencies for that input sequence.
    """
    # For the Teeling et al. method, the Z-scores require us to count
    # mono, di, tri and tetranucleotide sequences - these are stored
    # (in order) in the counts tuple
    counts = (
        collections.defaultdict(int),
        collections.defaultdict(int),
        collections.defaultdict(int),
        collections.defaultdict(int),
    )  # type: Tuple
    for rec in SeqIO.parse(filename, "fasta"):
        for seq in [str(rec.seq).upper(), str(rec.seq.reverse_complement()).upper()]:
            # The Teeling et al. algorithm requires us to consider
            # both strand orientations, so monocounts are easy
            for base in ("G", "C", "T", "A"):
                counts[0][base] += seq.count(base)
            # For di, tri and tetranucleotide counts, loop over the
            # sequence and its reverse complement, until near the end:
            for i in range(len(seq[:-4])):
                din, tri, tetra = seq[i : i + 2], seq[i : i + 3], seq[i : i + 4]
                counts[1][str(din)] += 1
                counts[2][str(tri)] += 1
                counts[3][str(tetra)] += 1
            # Then clean up the straggling bit at the end:
            counts[2][str(seq[-4:-1])] += 1
            counts[2][str(seq[-3:])] += 1
            counts[1][str(seq[-4:-2])] += 1
            counts[1][str(seq[-3:-1])] += 1
            counts[1][str(seq[-2:])] += 1
    # Following Teeling (2004), calculate expected frequencies for each
    # tetranucleotide; we ignore ambiguity symbols
    tetra_exp = {}
    for tet in [tetn for tetn in counts[3] if tetra_clean(tetn)]:
        tetra_exp[tet] = (
            1.0 * counts[2][tet[:3]] * counts[2][tet[1:]] / counts[1][tet[1:3]]
        )
    # Following Teeling (2004) we approximate the std dev and Z-score for each
    # tetranucleotide
    tetra_sd = {}
    tetra_z = {}
    for tet, exp in list(tetra_exp.items()):
        den = counts[1][tet[1:3]]
        tetra_sd[tet] = math.sqrt(
            exp * (den - counts[2][tet[:3]]) * (den - counts[2][tet[1:]]) / (den * den)
        )
        try:
            tetra_z[tet] = (counts[3][tet] - exp) / tetra_sd[tet]
        except ZeroDivisionError:
            # To record if we hit a zero in the estimation of variance
            # zeroes = [k for k, v in list(tetra_sd.items()) if v == 0]
            tetra_z[tet] = 1 / (counts[1][tet[1:3]] * counts[1][tet[1:3]])
    return tetra_z


# Returns true if the passed string contains only A, C, G or T
def tetra_clean(instr: str) -> bool:
    """Return True if string contains only unambiguous IUPAC nucleotide symbols.

    :param instr:  str, nucleotide sequence

    We are assuming that a low frequency of IUPAC ambiguity symbols doesn't
    affect our calculation.
    """
    if set(instr) - set("ACGT"):
        return False
    return True


# Calculate Pearson's correlation coefficient from the Z-scores for each
# tetranucleotide. If we're forcing rpy2, might as well use that, though...
def calculate_correlations(tetra_z: Dict[str, Dict[str, float]]) -> pd.DataFrame:
    """Return dataframe of Pearson correlation coefficients.

    :param tetra_z:  dict, Z-scores, keyed by sequence ID

    Calculates Pearson correlation coefficient from Z scores for each
    tetranucleotide. This is done longhand here, which is fast enough,
    but for robustness we might want to do something else... (TODO).

    Note that we report a correlation by this method, rather than a
    percentage identity.
    """
    orgs = sorted(tetra_z.keys())
    correlations = pd.DataFrame(index=orgs, columns=orgs, dtype=float).fillna(1.0)
    for idx, org1 in enumerate(orgs[:-1]):
        for org2 in orgs[idx + 1 :]:
            if not sorted(tetra_z[org1].keys()) == sorted(tetra_z[org2].keys()):
                raise AssertionError()
            tets = sorted(tetra_z[org1].keys())
            zscores = [
                [tetra_z[org1][t] for t in tets],
                [tetra_z[org2][t] for t in tets],
            ]
            zmeans = [sum(zscore) / len(zscore) for zscore in zscores]
            zdiffs = [
                [z - zmeans[0] for z in zscores[0]],
                [z - zmeans[1] for z in zscores[1]],
            ]
            diffprods = sum(
                [zdiffs[0][i] * zdiffs[1][i] for i in range(len(zdiffs[0]))]
            )
            zdiffs2 = [sum([z * z for z in zdiffs[0]]), sum([z * z for z in zdiffs[1]])]
            correlations.loc[org1, org2] = diffprods / math.sqrt(
                zdiffs2[0] * zdiffs2[1]
            )
            correlations.loc[org2, org1] = correlations.loc[org1, org2]
    return correlations
