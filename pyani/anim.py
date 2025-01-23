# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
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
# Copyright (c) 2016-2019 The James Hutton Institute
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

import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import intervaltree
from collections import defaultdict

from logging import Logger
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from . import pyani_config
from . import pyani_files
from . import pyani_jobs
from . import PyaniException

from .pyani_tools import ANIResults


class PyaniANImException(PyaniException):
    """ANIm-specific exception for pyani."""


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
    - non-executable file at passed path (this includes cases where the user doesn't
      have execute permissions on the file)
    - no version info returned
    """
    try:
        # Returns a TypeError if `nucmer_exe` is None
        try:
            nucmer_path = shutil.which(nucmer_exe)  # type:ignore
        except TypeError:
            return f"expected path to nucmer executable; received {nucmer_exe}"
        # Returns a TypeError if `nucmer_path` is not on the PATH
        nucmer_path = Path(nucmer_path)
    except TypeError:
        return f"{nucmer_exe} is not found in $PATH"

    if not nucmer_path.is_file():  # no executable
        return f"No nucmer executable at {nucmer_path}"

    # This should catch cases when the file can't be executed by the user
    if not os.access(nucmer_path, os.X_OK):  # file exists but not executable
        return f"nucmer exists at {nucmer_path} but not executable"

    cmdline = [nucmer_exe, "-V"]  # type: List
    result = subprocess.run(
        cmdline,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )

    # version information appears in different places for
    # different nucmer releases
    if result.stderr:  # expected to work for <= MUMmer3
        match = re.search(
            r"(?<=version\s)[0-9\.]*", str(result.stderr + result.stdout, "utf-8")
        )
    elif result.stdout:  # expected to work for MUMmer4
        match = re.search(r"[0-9a-z\.]*", str(result.stdout, "utf-8"))

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
    logger = logging.getLogger(__name__)

    ncmds, fcmds = generate_nucmer_commands(
        filenames, outdir, nucmer_exe, filter_exe, maxmatch
    )
    joblist = []
    for idx, ncmd in enumerate(ncmds):
        logger.info("Working with nucmer command: %s" % ncmd)
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

            # Comparisions A_vs_B
            ncmd, dcmd = construct_nucmer_cmdline(
                fname1, fname2, outdir, nucmer_exe, filter_exe, maxmatch
            )
            nucmer_cmdlines.append(ncmd)
            delta_filter_cmdlines.append(dcmd)

            # Comparions B_vs_A
            ncmd_rvs, dcmd_rvs = construct_nucmer_cmdline(
                fname2, fname1, outdir, nucmer_exe, filter_exe, maxmatch
            )
            nucmer_cmdlines.append(ncmd_rvs)
            delta_filter_cmdlines.append(dcmd_rvs)

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
    fname1, fname2 = [Path(fname1), Path(fname2)]

    # Compile commands
    # Nested output folders to avoid N^2 scaling in files-per-folder
    # Create folders incrementally (want an error if outdir does not exist)
    outsubdir = outdir / pyani_config.ALIGNDIR["ANIm"]
    outsubdir.mkdir(exist_ok=True)
    outsubdir = outdir / pyani_config.ALIGNDIR["ANIm"] / fname1.stem
    outsubdir.mkdir(exist_ok=True)
    outprefix = outsubdir / f"{fname1.stem}_vs_{fname2.stem}"
    if maxmatch:
        mode = "--maxmatch"
    else:
        mode = "--mum"
    nucmercmd = f"{nucmer_exe} {mode} -p {outprefix} {fname1} {fname2}"
    # There's a subtle pathlib.Path issue, here. We must use string concatenation to add suffixes
    # to the outprefix files, as using path.with_suffix() instead can replace part of the filestem
    # in those cases where there is a period in the stem (this occurs frequently as it is part
    # of the NCBI notation for genome assembly versions)
    filtercmd = (
        f"delta_filter_wrapper.py {filter_exe} -1 {str(outprefix) + '.delta'} "
        f"{str(outprefix) + '.filter'}"
    )
    return (nucmercmd, filtercmd)


def parse_delta(filename: Path) -> Tuple[int, int, float, int]:
    """Return (reference alignment length, query alignment length, average identity, similarity erors)

    :param filename: Path to the input .delta file

    Calculates the aligned lengths for reference and query and average nucleotide
    identity, and returns the cumulative total for each as a tuple.

    The delta file format contains seven numbers in the lines of interest:
    see http://mummer.sourceforge.net/manual/ for specification

    - start on query
    - end on query
    - start on target
    - end on target
    - error count (non-identical, plus indels)
    - similarity errors (non-positive match scores)
        [NOTE: with PROmer this is equal to error count]
    - stop codons (always zero for nucmer)

    We report ANIm identity by finding an average across all alignments using
    the following formula:

    sum of weighted identical bases / sum of aligned bases from each fragment

    For example:

    reference.fasta query.fasta
    NUCMER
    >ref_seq_A ref_seq_B 40 40
    1 10 1 11 5 5 0
    -1
    0
    15 20 25 30 0 0 0

    The delta file tells us there are two alignments. The first alignment runs from base 1
    to base 10 in the reference sequence, and from base 1 to 11 in the query sequence
    with a similarity error of 5. The second alignment runs from base 15 to 20 in
    the reference, and base 25 to 30 in the query with 0 similarity errors. To calculate
    the %ID, we can:

    - Find the number of all aligned bases from each sequence:
    aligned reference bases region 1 = 10 - 1 + 1 = 10
    aligned query bases region 1 = 11 - 1 + 1 = 11
    aligned reference bases region 2 = 20 - 15 + 1 = 6
    aligned query bases region 2 = 30 - 25 + 1 = 6

    - Find weighted identical bases
    alignment 1 identity weighted = (10 + 11) - (2 * 5) = 11
    alignment 2 identity weighted = (6 + 6) - (2 * 0) = 12

    - Calculate %ID
    (11 + 12) / (10 + 11 + 6 + 6) = 0.696969696969697

    To calculate alignment lengths, we extract the regions of each alignment
    (either for query or reference) provided in the .delta file and merge the overlapping
    regions with IntervalTree. Then, we calculate the total sum of all aligned regions.
    """

    current_ref, current_qry, raln_length, qaln_length, sim_error, avrg_ID = (
        None,
        None,
        0,
        0,
        0,
        0.0,
    )

    regions_ref = defaultdict(list)  # Hold a dictionary for query regions
    regions_qry = defaultdict(list)  # Hold a dictionary for query regions

    aligned_bases = []  # Hold a list for aligned bases for each sequence
    weighted_identical_bases = []  # Hold a list for weighted identical bases

    for line in [_.strip().split() for _ in filename.open("r").readlines()]:
        if line[0] == "NUCMER":  # Skip headers
            continue
        # Lines starting with ">" indicate which sequences are aligned
        if line[0].startswith(">"):
            current_ref = line[0].strip(">")
            current_qry = line[1]
        # Lines with seven columns are alignment region headers:
        if len(line) == 7:
            # Obtaining aligned regions needed to check for overlaps
            regions_ref[current_ref].append(
                tuple(sorted(list([int(line[0]), int(line[1])])))
            )  # aligned regions reference
            regions_qry[current_qry].append(
                tuple(sorted(list([int(line[2]), int(line[3])])))
            )  # aligned regions qry

            # Calculate aligned bases for each sequence
            ref_aln_lengths = abs(int(line[1]) - int(line[0])) + 1
            qry_aln_lengths = abs(int(line[3]) - int(line[2])) + 1
            aligned_bases.append(ref_aln_lengths)
            aligned_bases.append(qry_aln_lengths)

            # Calculate weighted identical bases
            sim_error += int(line[4])
            weighted_identical_bases.append(
                (ref_aln_lengths + qry_aln_lengths) - (2 * int(line[4]))
            )

    # Calculate average %ID
    avrg_ID = sum(weighted_identical_bases) / sum(aligned_bases)

    # Calculate total aligned bases (no overlaps)
    for seq_id in regions_qry:
        qry_tree = intervaltree.IntervalTree.from_tuples(regions_qry[seq_id])
        qry_tree.merge_overlaps(strict=False)
        for interval in qry_tree:
            qaln_length += interval.end - interval.begin + 1

    for seq_id in regions_ref:
        ref_tree = intervaltree.IntervalTree.from_tuples(regions_ref[seq_id])
        ref_tree.merge_overlaps(strict=False)
        for interval in ref_tree:
            raln_length += interval.end - interval.begin + 1

    return (raln_length, qaln_length, avrg_ID, sim_error)


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
    logger = logging.getLogger(__name__)

    # Process directory to identify input files - as of v0.2.4 we use the
    # .filter files that result from delta-filter (1:1 alignments)
    deltafiles = sorted(delta_dir.glob("*/*.filter"))

    logger.info("%s has %d files to load", delta_dir, len(deltafiles))
    if not deltafiles:
        logger.error("%s empty? No filter files found", delta_dir)
        raise PyaniANImException(f"{delta_dir} contains no filter files.")

    # Hold data in ANIResults object
    results = ANIResults(list(org_lengths.keys()), "ANIm")

    # Fill diagonal NA values for alignment_length with org_lengths
    for org, length in list(org_lengths.items()):
        results.alignment_lengths.loc[org, org] = length

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
        (
            query_tot_length,
            subject_tot_length,
            weighted_identity,
            tot_sim_error,
        ) = parse_delta(deltafile)
        if subject_tot_length == 0 and logger is not None:
            if logger:
                logger.warning(
                    "Total alignment length reported in %s is zero!", deltafile
                )
            sys.exit("Zero length alignment!")
        query_cover = float(query_tot_length) / org_lengths[qname]
        sbjct_cover = float(subject_tot_length) / org_lengths[sname]
        perc_id = weighted_identity

        # Populate dataframes: when assigning data from symmetrical MUMmer
        # output, both upper and lower triangles will be populated
        results.add_tot_length(qname, sname, query_tot_length, subject_tot_length)
        results.add_sim_errors(qname, sname, tot_sim_error)
        results.add_pid(qname, sname, perc_id)
        results.add_coverage(qname, sname, query_cover, sbjct_cover)
    return results
