# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019
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
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019 University of Strathclyde
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
import platform
import re
import shutil
import subprocess

from logging import Logger
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from . import pyani_config
from . import pyani_files
from . import pyani_jobs

from .pyani_tools import ANIResults


# Get a list of FASTA files from the input directory
def get_fasta_files(dirname: Path = Path(".")) -> Iterable:
    """Return iterable of FASTA files in the passed directory.

    :param dirname:  str, path to input directory
    """
    infiles = pyani_files.get_input_files(
        dirname, ".fasta", ".fas", ".fa", ".fna", ".fsa_nt"
    )
    return infiles


# Get NUCmer version
def get_version(nucmer_exe: Path = pyani_config.NUCMER_DEFAULT) -> str:
    """Return NUCmer package version as a string.

    :param nucmer_exe:  path to NUCmer executable

    We expect NUCmer to return a string on STDERR as

    .. code-block:: bash

        $ nucmer
        NUCmer (NUCleotide MUMmer) version 3.1

    we concatenate this with the OS name.

    The following circumstances are explicitly reported as strings

    - no executable at passed path
    - non-executable file at passed path
    - no version info returned
    """
    nucmer_path = Path(shutil.which(nucmer_exe))  # type:ignore

    if not os.path.isfile(nucmer_path):  # no executable
        return f"No nucmer executable at {nucmer_path}"

    if not os.access(nucmer_path, os.X_OK):  # file exists but not executable
        return f"nucmer exists at {nucmer_path} but not executable"

    cmdline = [nucmer_exe, "-V"]  # type: List
    result = subprocess.run(
        cmdline, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    match = re.search(r"(?<=version\s)[0-9\.]*", str(result.stderr, "utf-8"))
    version = match.group()  # type: ignore

    if 0 == len(version.strip()):
        return f"nucmer exists at {nucmer_path} but could not retrieve version"

    return f"{platform.system()}_{version} ({nucmer_path})"


# Generate list of Job objects, one per NUCmer run
def generate_nucmer_jobs(
    filenames: List[Path],
    outdir: Path = Path("."),
    nucmer_exe: Path = pyani_config.NUCMER_DEFAULT,
    filter_exe: Path = pyani_config.FILTER_DEFAULT,
    maxmatch: bool = False,
    jobprefix: str = "ANINUCmer",
):
    """Return list of Jobs describing NUCmer command-lines for ANIm.

    :param filenames:  Iterable, Paths to input FASTA files
    :param outdir:  str, path to output directory
    :param nucmer_exe:  str, location of the nucmer binary
    :param filter_exe:
    :param maxmatch:  Boolean flag indicating to use NUCmer's -maxmatch option
    :param jobprefix:

    Loop over all FASTA files, generating Jobs describing NUCmer command lines
    for each pairwise comparison.
    """
    ncmds, fcmds = generate_nucmer_commands(
        filenames, outdir, nucmer_exe, filter_exe, maxmatch
    )
    joblist = []
    for idx, ncmd in enumerate(ncmds):
        njob = pyani_jobs.Job(f"{jobprefix}_{idx:06d}-n", ncmd)
        fjob = pyani_jobs.Job(f"{jobprefix}_{idx:06d}-f", fcmds[idx])
        fjob.add_dependency(njob)
        joblist.append(fjob)
    return joblist


# Generate list of NUCmer pairwise comparison command lines from
# passed sequence filenames
def generate_nucmer_commands(
    filenames: List[Path],
    outdir: Path = Path("."),
    nucmer_exe: Path = pyani_config.NUCMER_DEFAULT,
    filter_exe: Path = pyani_config.FILTER_DEFAULT,
    maxmatch: bool = False,
) -> Tuple[List, List]:
    """Return list of NUCmer command-lines for ANIm.

    :param filenames:  a list of paths to input FASTA files
    :param outdir:  path to output directory
    :param nucmer_exe:  location of the nucmer binary
    :param maxmatch:  Boolean flag indicating to use NUCmer's -maxmatch option

    The first element returned is a list of NUCmer commands, and the
    second a corresponding list of delta_filter_wrapper.py commands.
    The NUCmer commands should each be run before the corresponding
    delta-filter command.

    TODO: This return value needs to be reworked as a collection.

    Loop over all FASTA files generating NUCmer command lines for each
    pairwise comparison.
    """
    nucmer_cmdlines, delta_filter_cmdlines = [], []
    filenames = sorted(filenames)  # enforce ordering of filenames
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
    fname1: Path,
    fname2: Path,
    outdir: Path = Path("."),
    nucmer_exe: Path = pyani_config.NUCMER_DEFAULT,
    filter_exe: Path = pyani_config.FILTER_DEFAULT,
    maxmatch: bool = False,
) -> Tuple[str, str]:
    """Return a tuple of corresponding NUCmer and delta-filter commands.

    :param fname1:  path to query FASTA file
    :param fname2:  path to subject FASTA file
    :param outdir:  path to output directory
    :param nucmer_exe:
    :param filter_exe:
    :param maxmatch:  Boolean flag indicating whether to use NUCmer's -maxmatch
    option. If not, the -mum option is used instead

    The split into a tuple was made necessary by changes to SGE/OGE.
    The delta-filter command must now be run as a dependency of the NUCmer
    command, and be wrapped in a Python script to capture STDOUT.

    NOTE: This command-line writes output data to a subdirectory of the passed
    outdir, called "nucmer_output".
    """
    # Cast path strings to pathlib.Path for safety
    fname1, fname2 = Path(fname1), Path(fname2)

    # Compile commands
    outsubdir = outdir / pyani_config.ALIGNDIR["ANIm"]
    outprefix = outsubdir / f"{fname1.stem}_vs_{fname2.stem}"
    if maxmatch:
        mode = "--maxmatch"
    else:
        mode = "--mum"
    nucmercmd = "{0} {1} -p {2} {3} {4}".format(
        nucmer_exe, mode, outprefix, fname1, fname2
    )
    # There's a subtle pathlib.Path issue, here. We must use string concatenation to add suffixes
    # to the outprefix files, as using path.with_suffix() instead can replace part of the filestem
    # in those cases where there is a period in the stem (this occurs frequently as it is part
    # of the NCBI notation for genome assembly versions)
    filtercmd = (
        f"delta_filter_wrapper.py {filter_exe} -1 {str(outprefix) + '.delta'} "
        f"{str(outprefix) + '.filter'}"
    )
    return (nucmercmd, filtercmd)


# Parse NUCmer delta file to get total alignment length and total sim_errors
def parse_delta(filename: Path) -> Tuple[int, int]:
    """Return (alignment length, similarity errors) tuple from passed .delta.

    :param filename:  Path, path to the input .delta file

    Extracts the aligned length and number of similarity errors for each
    aligned uniquely-matched region, and returns the cumulative total for
    each as a tuple.

    Similarity errors are defined in the .delta file spec (see below) as
    non-positive match scores. For NUCmer output, this is identical to the
    number of errors (non-identities and indels).

    Delta file format has seven numbers in the lines of interest:
    see http://mummer.sourceforge.net/manual/ for specification

    - start on query
    - end on query
    - start on target
    - end on target
    - error count (non-identical, plus indels)
    - similarity errors (non-positive match scores)
        [NOTE: with PROmer this is equal to error count]
    - stop codons (always zero for nucmer)

    To calculate alignment length, we take the length of the aligned region of
    the reference (no gaps), and process the delta information. This takes the
    form of one value per line, following the header sequence. Positive values
    indicate an insertion in the reference; negative values a deletion in the
    reference (i.e. an insertion in the query). The total length of the alignment
    is then:

    reference_length + insertions - deletions

    For example:

    A = ABCDACBDCAC$
    B = BCCDACDCAC$
    Delta = (1, -3, 4, 0)
    A = ABC.DACBDCAC$
    B = .BCCDAC.DCAC$

    A is the reference and has length 11. There are two insertions (positive delta),
    and one deletion (negative delta). Alignment length is then 11 + 1 = 12.
    """
    in_aln, aln_length, sim_errors = False, 0, 0
    for line in [_.strip().split() for _ in filename.open("r").readlines()]:
        if line[0] == "NUCMER" or line[0].startswith(">"):  # Skip headers
            continue
        # Lines with seven columns are alignment region headers:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0])) + 1  # reference length
            sim_errors += int(line[4])  # count of non-identities and indels
            in_aln = True
        # Lines with a single column (following a header) report numbers of symbols
        # until next insertion (+ve) or deletion (-ve) in the reference; one line per
        # insertion/deletion; the alignment always ends with 0
        if in_aln and line[0].startswith("0"):
            in_aln = False
        elif in_aln:
            # Add one to the alignment length for each reference insertion; subtract
            # one for each deletion
            val = int(line[0])
            if val < 1:  # deletion in reference
                aln_length += 1
            elif val == 0:  # ends the alignment entry
                in_aln = False
    return aln_length, sim_errors


# Parse all the .delta files in the passed directory
def process_deltadir(
    delta_dir: Path, org_lengths: Dict, logger: Optional[Logger] = None
) -> ANIResults:
    """Return tuple of ANIm results for .deltas in passed directory.

    :param delta_dir:  Path, path to the directory containing .delta files
    :param org_lengths:  dictionary of total sequence lengths, keyed by sequence

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
        qname, sname = deltafile.stem.split("_vs_")

        # We may have .delta files from other analyses in the same directory
        # If this occurs, we raise a warning, and skip the .delta file
        if qname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Query name %s not in input sequence list, skipping %s",
                    qname,
                    deltafile,
                )
            continue
        if sname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Subject name %s not in input sequence list, skipping %s",
                    sname,
                    deltafile,
                )
            continue
        tot_length, tot_sim_error = parse_delta(deltafile)
        if tot_length == 0 and logger is not None:
            if logger:
                logger.warning(
                    "Total alignment length reported in %s is zero!", deltafile
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
