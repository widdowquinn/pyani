# -*- coding: utf-8 -*-
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

(c) The James Hutton Institute 2016-2019
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2016-2019 The James Hutton Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
import re
import subprocess

from . import pyani_config
from . import pyani_files
from . import pyani_jobs

from .pyani_tools import ANIResults


# Get a list of FASTA files from the input directory
def get_fasta_files(dirname=None):
    """Return list of FASTA files in the passed directory.

    - dirname - path to input directory
    """
    if dirname is None:
        dirname = "."
    infiles = pyani_files.get_input_files(
        dirname, ".fasta", ".fas", ".fa", ".fna", ".fsa_nt"
    )
    return infiles


# Get NUCmer version
def get_version(nucmer_exe=pyani_config.NUCMER_DEFAULT):
    """Return NUCmer package version as a string.

    We expect NUCmer to return a string on STDERR as

    nucmer
    NUCmer (NUCleotide MUMmer) version 3.1
    """
    cmdline = [nucmer_exe, "-V"]
    result = subprocess.run(
        cmdline, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    return re.search(r"(?<=version\s)[0-9\.]*", str(result.stderr, "utf-8")).group()


# Generate list of Job objects, one per NUCmer run
def generate_nucmer_jobs(
    filenames,
    outdir=".",
    nucmer_exe=pyani_config.NUCMER_DEFAULT,
    filter_exe=pyani_config.FILTER_DEFAULT,
    maxmatch=False,
    jobprefix="ANINUCmer",
):
    """Return list of Jobs describing NUCmer command-lines for ANIm.

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files, generating Jobs describing NUCmer command lines
    for each pairwise comparison.
    """
    ncmds, fcmds = generate_nucmer_commands(
        filenames, outdir, nucmer_exe, filter_exe, maxmatch
    )
    joblist = []
    for idx, ncmd in enumerate(ncmds):
        njob = pyani_jobs.Job("%s_%06d-n" % (jobprefix, idx), ncmd)
        fjob = pyani_jobs.Job("%s_%06d-f" % (jobprefix, idx), fcmds[idx])
        fjob.add_dependency(njob)
        joblist.append(fjob)
    return joblist


# Generate list of NUCmer pairwise comparison command lines from
# passed sequence filenames
def generate_nucmer_commands(
    filenames,
    outdir=".",
    nucmer_exe=pyani_config.NUCMER_DEFAULT,
    filter_exe=pyani_config.FILTER_DEFAULT,
    maxmatch=False,
):
    """Return list of NUCmer command-lines for ANIm.

    The first element returned is a list of NUCmer commands, and the
    second a corresponding list of delta_filter_wrapper.py commands.
    The NUCmer commands should each be run before the corresponding
    delta-filter command.

    TODO: This return value needs to be reworked as a collection.

    - filenames - a list of paths to input FASTA files
    - outdir - path to output directory
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files generating NUCmer command lines for each
    pairwise comparison.
    """
    nucmer_cmdlines, delta_filter_cmdlines = [], []
    for idx, fname1 in enumerate(filenames[:-1]):
        for fname2 in filenames[idx + 1 :]:
            ncmd, dcmd = construct_nucmer_cmdline(
                fname1, fname2, outdir, nucmer_exe, filter_exe, maxmatch
            )
            nucmer_cmdlines.append(ncmd)
            delta_filter_cmdlines.append(dcmd)
    return (nucmer_cmdlines, delta_filter_cmdlines)


# Generate single NUCmer pairwise comparison command line from pair of
# input filenames
def construct_nucmer_cmdline(
    fname1,
    fname2,
    outdir=".",
    nucmer_exe=pyani_config.NUCMER_DEFAULT,
    filter_exe=pyani_config.FILTER_DEFAULT,
    maxmatch=False,
):
    """Return a tuple of corresponding NUCmer and delta-filter commands

    The split into a tuple was made necessary by changes to SGE/OGE.
    The delta-filter command must now be run as a dependency of the NUCmer
    command, and be wrapped in a Python script to capture STDOUT.

    NOTE: This command-line writes output data to a subdirectory of the passed
    outdir, called "nucmer_output".

    - fname1 - query FASTA filepath
    - fname2 - subject FASTA filepath
    - outdir - path to output directory
    - maxmatch - Boolean flag indicating whether to use NUCmer's -maxmatch
    option. If not, the -mum option is used instead
    """
    outsubdir = os.path.join(outdir, pyani_config.ALIGNDIR["ANIm"])
    outprefix = os.path.join(
        outsubdir,
        "%s_vs_%s"
        % (
            os.path.splitext(os.path.split(fname1)[-1])[0],
            os.path.splitext(os.path.split(fname2)[-1])[0],
        ),
    )
    if maxmatch:
        mode = "--maxmatch"
    else:
        mode = "--mum"
    nucmercmd = "{0} {1} -p {2} {3} {4}".format(
        nucmer_exe, mode, outprefix, fname1, fname2
    )
    filtercmd = "delta_filter_wrapper.py " + "{0} -1 {1} {2}".format(
        filter_exe, outprefix + ".delta", outprefix + ".filter"
    )
    return (nucmercmd, filtercmd)


# Parse NUCmer delta file to get total alignment length and total sim_errors
def parse_delta(filename):
    """Return (alignment length, similarity errors) tuple from passed .delta.

    - filename - path to the input .delta file

    Extracts the aligned length and number of similarity errors for each
    aligned uniquely-matched region, and returns the cumulative total for
    each as a tuple.

    Delta file format has seven numbers in the lines of interest:
    - start on query
    - end on query
    - start on target
    - end on target
    - error count (non-identical, plus indels)
    - similarity errors (non-positive match scores)
    - stop codons (always zero for nucmer)
    """
    aln_length, sim_errors = 0, 0
    for line in [l.strip().split() for l in open(filename, "r").readlines()]:
        if line[0] == "NUCMER" or line[0].startswith(">"):  # Skip headers
            continue
        # We only process lines with seven columns:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0]) + 1)
            indels = int(line[4]) - int(line[5])
            aln_length -= indels
            sim_errors += int(line[5])
    return aln_length, sim_errors


# Parse all the .delta files in the passed directory
def process_deltadir(delta_dir, org_lengths, logger=None):
    """Return tuple of ANIm results for .deltas in passed directory.

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
    deltafiles = pyani_files.get_input_files(delta_dir, ".filter")

    # Hold data in ANIResults object
    results = ANIResults(list(org_lengths.keys()), "ANIm")

    # Fill diagonal NA values for alignment_length with org_lengths
    for org, length in list(org_lengths.items()):
        results.alignment_lengths[org][org] = length

    # Process .delta files assuming that the filename format holds:
    # org1_vs_org2.delta
    for deltafile in deltafiles:
        qname, sname = os.path.splitext(os.path.split(deltafile)[-1])[0].split("_vs_")

        # We may have .delta files from other analyses in the same directory
        # If this occurs, we raise a warning, and skip the .delta file
        if qname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Query name %s not in input " % qname
                    + "sequence list, skipping %s" % deltafile
                )
            continue
        if sname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Subject name %s not in input " % sname
                    + "sequence list, skipping %s" % deltafile
                )
            continue
        tot_length, tot_sim_error = parse_delta(deltafile)
        if tot_length == 0 and logger is not None:
            if logger:
                logger.warning(
                    "Total alignment length reported in " + "%s is zero!" % deltafile
                )
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
