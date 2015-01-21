# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to implement the ANIm average nucleotide identity method.

Calculates ANI by the ANIm method, as described in Richter et al (2009)
Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

All input FASTA format files are compared against each other, pairwise,
using NUCmer (binary location must be provided). NUCmer output will be stored
in a specified output directory.

The NUCmer .delta file output is parsed to obtain an alignment length
and similarity error count for every unique region alignment. These are
processed to give matrices of aligned sequence lengths, similarity error
counts, average nucleotide identity (ANI) percentages, and minimum aligned
percentage (of whole genome) for each pairwise comparison.
"""

import pandas as pd

import os

import pyani_config
import pyani_files


# Generate list of NUCmer pairwise comparison command lines from
# passed sequence filenames
def generate_nucmer_commands(filenames, outdir='.',
                             nucmer_exe=pyani_config.NUCMER_DEFAULT,
                             maxmatch=False):
    """Return a list of NUCmer command-lines for ANIm

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files generating NUCmer command lines for each
    pairwise comparison.
    """
    cmdlines = []
    for idx, f1 in enumerate(filenames[:-1]):
        cmdlines.extend([construct_nucmer_cmdline(f1, f2, outdir,
                                                  nucmer_exe) for
                         f2 in filenames[idx+1:]])
    return cmdlines


# Generate single NUCmer pairwise comparison command line from pair of
# input filenames
def construct_nucmer_cmdline(fname1, fname2, outdir='.',
                             nucmer_exe=pyani_config.NUCMER_DEFAULT,
                             maxmatch=False):
    """Returns a single NUCmer pairwise comparison command.

    - fname1 - query FASTA filepath
    - fname2 - subject FASTA filepath
    - outdir - path to output directory
    - maxmatch - Boolean flag indicating whether to use NUCmer's -maxmatch
    option. If not, the -mum option is used instead
    """
    outprefix = os.path.join(outdir, "%s_vs_%s" %
                             (os.path.splitext(os.path.split(fname1)[-1])[0],
                              os.path.splitext(os.path.split(fname2)[-1])[0]))
    if maxmatch:
        mode = "-maxmatch"
    else:
        mode = "-mum"
    return "{0} {1} -p {2} {3} {4}".format(nucmer_exe, mode, outprefix,
                                           fname1, fname2)


# Parse NUCmer delta file to get total alignment length and total sim_errors
def parse_delta(filename):
    """Returns (alignment length, similarity errors) tuple from passed .delta.

    - filename - path to the input .delta file

    Extracts the aligned length and number of similarity errors for each
    aligned uniquely-matched region, and returns the cumulative total for
    each as a tuple.
    """
    aln_length, sim_errors = 0, 0
    for line in [l.strip().split() for l in open(filename, 'rU').readlines()]:
        if line[0] == 'NUCMER' or line[0].startswith('>'):  # Skip headers
            continue
        # We only process lines with seven columns:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0]))
            sim_errors += int(line[4])
    return aln_length, sim_errors


# Parse all the .delta files in the passed directory
def process_deltadir(delta_dir, org_lengths):
    """Returns a tuple of ANIm results for .deltas in passed directory.

    - delta_dir - path to the directory containing .delta files
    - org_lengths - dictionary of total sequence lengths, keyed by sequence

    Returns the following pandas dataframes in a tuple; query sequences are
    rows, subject sequences are columns:

    - alignment_lengths - symmetrical: total length of alignment
    - percentage_identity - symmetrical: percentage identity of alignment
    - alignment_coverage - non-symmetrical: coverage of query and subject
    - similarity_errors - symmetrical: count of similarity errors

    May throw a ZeroDivisionError if one or more NUCmer runs failed, or a
    very distant sequence was included in the analysis.
    """
    # Process directory to identify input files
    deltafiles = pyani_files.get_input_files(delta_dir, '.delta')
    labels = org_lengths.keys()
    # Hold data in pandas dataframe
    alignment_lengths = pd.DataFrame(index=labels, columns=labels,
                                     dtype=float)
    similarity_errors = pd.DataFrame(index=labels, columns=labels,
                                     dtype=float).fillna(0)
    percentage_identity = pd.DataFrame(index=labels, columns=labels,
                                       dtype=float).fillna(1.0)
    alignment_coverage = pd.DataFrame(index=labels, columns=labels,
                                      dtype=float).fillna(1.0)
    # Fill diagonal NA values for alignment_length with org_lengths
    for org, length in org_lengths.items():
        alignment_lengths[org][org] = length
    # Process .delta files assuming that the filename format holds:
    # org1_vs_org2.delta
    for deltafile in deltafiles:
        qname, sname = \
            os.path.splitext(os.path.split(deltafile)[-1])[0].split('_vs_')
        tot_length, tot_sim_error = parse_delta(deltafile)
        query_cover = float(tot_length) / org_lengths[qname]
        sbjct_cover = float(tot_length) / org_lengths[sname]
        # Calculate percentage ID of aligned length. This may fail if
        # total length is zero.
        # The ZeroDivisionError that would arise should be handled
        # Common causes are that a NUCmer run failed, or that a very
        # distant sequence was included in the analysis.
        perc_id = 1 - float(tot_sim_error) / tot_length
        # Populate dataframes: when assigning data, pandas dataframes
        # take column, index order, i.e. df['column']['row'] - this only
        # matters for asymmetrical data
        alignment_lengths.loc[qname, sname] = tot_length
        alignment_lengths.loc[sname, qname] = tot_length
        similarity_errors.loc[qname, sname] = tot_sim_error
        similarity_errors.loc[sname, qname] = tot_sim_error
        percentage_identity.loc[qname, sname] = perc_id
        percentage_identity.loc[sname, qname] = perc_id
        alignment_coverage.loc[sname, qname] = query_cover
        alignment_coverage.loc[qname, sname] = sbjct_cover
    return(alignment_lengths, percentage_identity, alignment_coverage,
           similarity_errors)
