# Copyright 2013-2016, The James Hutton Insitute
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

import os

from . import pyani_config
from . import pyani_files
from . import pyani_jobs
from .pyani_tools import ANIResults


# Generate list of Job objects, one per NUCmer run
def generate_nucmer_jobs(filenames, outdir='.',
                         nucmer_exe=pyani_config.NUCMER_DEFAULT,
                         filter_exe=pyani_config.FILTER_DEFAULT,
                         maxmatch=False,
                         jobprefix="ANINUCmer"):
    """Return a list of Jobs describing NUCmer command-lines for ANIm

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files, generating Jobs describing NUCmer command lines
    for each pairwise comparison.
    """
    cmdlines = generate_nucmer_commands(filenames, outdir, nucmer_exe,
                                        filter_exe, maxmatch)
    joblist = []
    for idx, cmd in enumerate(cmdlines):
        joblist.append(pyani_jobs.Job("%s_%06d" % (jobprefix, idx), cmd))
    return joblist


# Generate list of NUCmer pairwise comparison command lines from
# passed sequence filenames
def generate_nucmer_commands(filenames, outdir='.',
                             nucmer_exe=pyani_config.NUCMER_DEFAULT,
                             filter_exe=pyani_config.FILTER_DEFAULT,
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
    for idx, fname1 in enumerate(filenames[:-1]):
        cmdlines.extend([construct_nucmer_cmdline(fname1, fname2, outdir,
                                                  nucmer_exe, filter_exe,
                                                  maxmatch) for
                         fname2 in filenames[idx+1:]])
    return cmdlines


# Generate single NUCmer pairwise comparison command line from pair of
# input filenames
def construct_nucmer_cmdline(fname1, fname2, outdir='.',
                             nucmer_exe=pyani_config.NUCMER_DEFAULT,
                             filter_exe=pyani_config.FILTER_DEFAULT,
                             maxmatch=False):
    """Returns a single NUCmer pairwise comparison command.

    NOTE: This command-line writes output data to a subdirectory of the passed
    outdir, called "nucmer_output".

    - fname1 - query FASTA filepath
    - fname2 - subject FASTA filepath
    - outdir - path to output directory
    - maxmatch - Boolean flag indicating whether to use NUCmer's -maxmatch
    option. If not, the -mum option is used instead
    """
    outsubdir = os.path.join(outdir, pyani_config.ALIGNDIR['ANIm'])
    outprefix = os.path.join(outsubdir, "%s_vs_%s" %
                             (os.path.splitext(os.path.split(fname1)[-1])[0],
                              os.path.splitext(os.path.split(fname2)[-1])[0]))
    if maxmatch:
        mode = "--maxmatch"
    else:
        mode = "--mum"
    nucmercmd = "{0} {1} -p {2} {3} {4}".format(nucmer_exe, mode, outprefix,
                                                fname1, fname2)
    filtercmd = "{0} -1 {1} > {2}".format(filter_exe,
                                          outprefix + '.delta',
                                          outprefix + '.filter')
    return "{0}; {1}".format(nucmercmd, filtercmd)


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
def process_deltadir(delta_dir, org_lengths, logger=None):
    """Returns a tuple of ANIm results for .deltas in passed directory.

    - delta_dir - path to the directory containing .delta files
    - org_lengths - dictionary of total sequence lengths, keyed by sequence

    Returns the following pandas dataframes in an ANIResults object;
    query sequences are rows, subject sequences are columns:

    - alignment_lengths - symmetrical: total length of alignment
    - percentage_identity - symmetrical: percentage identity of alignment
    - alignment_coverage - non-symmetrical: coverage of query and subject
    - similarity_errors - symmetrical: count of similarity errors

    May throw a ZeroDivisionError if one or more NUCmer runs failed, or a
    very distant sequence was included in the analysis.
    """
    # Process directory to identify input files - as of v0.2.4 we use the
    # .filter files that result from delta-filter (1:1 alignments)
    deltafiles = pyani_files.get_input_files(delta_dir, '.filter')

    # Hold data in ANIResults object
    results = ANIResults(list(org_lengths.keys()), "ANIm")

    # Fill diagonal NA values for alignment_length with org_lengths
    for org, length in list(org_lengths.items()):
        results.alignment_lengths[org][org] = length

    # Process .delta files assuming that the filename format holds:
    # org1_vs_org2.delta
    for deltafile in deltafiles:
        qname, sname = \
            os.path.splitext(os.path.split(deltafile)[-1])[0].split('_vs_')

        # We may have .delta files from other analyses in the same directory
        # If this occurs, we raise a warning, and skip the .delta file
        if qname not in list(org_lengths.keys()):
            if logger:
                logger.warning("Query name %s not in input " % qname +
                               "sequence list, skipping %s" % deltafile)
            continue
        if sname not in list(org_lengths.keys()):
            if logger:
                logger.warning("Subject name %s not in input " % sname +
                               "sequence list, skipping %s" % deltafile)
            continue
        tot_length, tot_sim_error = parse_delta(deltafile)
        if tot_length == 0 and logger is not None:
            if logger:
                logger.warning("Total alignment length reported in " +
                               "%s is zero!" % deltafile)
        query_cover = float(tot_length) / org_lengths[qname]
        sbjct_cover = float(tot_length) / org_lengths[sname]

        # Calculate percentage ID of aligned length. This may fail if
        # total length is zero.
        # The ZeroDivisionError that would arise should be handled
        # Common causes are that a NUCmer run failed, or that a very
        # distant sequence was included in the analysis.
        try:
            perc_id = 1 - float(tot_sim_error) / tot_length
        except ZeroDivisionError:
            perc_id = 0  # set arbitrary value of zero identity
            results.zero_error = True

        # Populate dataframes: when assigning data from symmetrical MUMmer
        # output, both upper and lower triangles will be populated
        results.add_tot_length(qname, sname, tot_length)
        results.add_sim_errors(qname, sname, tot_sim_error)
        results.add_pid(qname, sname, perc_id)
        results.add_coverage(qname, sname, query_cover, sbjct_cover)
    return results
