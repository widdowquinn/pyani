# -*- coding: utf-8 -*-
# (c) The University of Strathclyde 2021–Present
# Author: Bailey Harrington
#
# Contact: bailey.harrington@strath.ac.uk
#
# Bailey Harrington,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G4 0RE,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2021–Present University of Strathclyde
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
"""Code to implement the fastANI average nucleotide identity method."""

import logging
import platform
import os
import re
import shutil
import subprocess

from logging import Logger
from pathlib import Path
from typing import NamedTuple, Callable, Dict, List, Optional, Tuple

# import pandas as pd

from Bio import SeqIO

from . import pyani_config
from . import pyani_files
from . import pyani_jobs
from . import PyaniException
from .pyani_tools import ANIResults, BLASTcmds, BLASTexes, BLASTfunctions


class PyaniFastANIException(PyaniException):

    """Exception raised when there is a problem with fastANI"""


class ComparisonResult(NamedTuple):
    reference: Path
    query: Path
    ani: float
    matches: int
    fragments: int


def get_version(fastani_exe: Path = pyani_config.FASTANI_DEFAULT) -> str:
    """Return FastANI package version as a string.

    :param fastani_exe: path to FastANI executable

    We expect fastANI to return a string on STDOUT as

    .. code-block:: bash

        $ ./fastANI -v
        version 1.32

    we concatenate this with the OS name.

    The following circumstances are explicitly reported as strings:

    - no executable at passed path
    - non-executable file at passed path (this includes cases where the user doesn't have execute permissions on the file)
    - no version info returned
    """
    try:
        fastani_path = Path(shutil.which(fastani_exe))  # type:ignore
    # Returns a TypeError if `fastani_exe` is None
    except TypeError:
        return f"expected file location; received {fastani_exe}"

    # If a string that is not an executable is passed to
    # shutil.which(), the return value will be None, so
    # this check is still needed
    if fastani_path is None:
        return f"{fastani_exe} is not found in $PATH"

    if not fastani_path.is_file():  # no executable
        return f"No fastANI executable at {fastani_path}"

    # This should catch cases when the file can't be executed by the user
    if not os.access(fastani_path, os.X_OK):  # file exists but not executable
        return f"fastANI exists at {fastani_path} but not executable"

    cmdline = [fastani_exe, "-v"]  # type: List
    result = subprocess.run(
        cmdline,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )  # type CompletedProcess

    match = re.search(
        r"(?<=version\s)[0-9\.]*", str(result.stderr + result.stdout, "utf-8")
    )
    version = match.group()  # type: ignore

    if 0 == len(version.strip()):
        return f"fastANI exists at {fastani_path} but could not retrieve version"

    return f"{platform.system()}_{version} ({fastani_path})"


# Generate list of Job objects, one per fastANI run -> this is maybe not necessary
def generate_fastani_jobs(
    filenames: List[Path],
    outdir: Path = Path("."),
    fastani_exe: Path = pyani_config.FASTANI_DEFAULT,
    fragLen: int = 3000,
    kmerSize: int = 16,
    minFraction: float = 0.2,
    jobprefix: str = "fastANI",
):  # should this have a hint for return type?
    """Return list of Jobs describing fastANI command lines.

    :param filenames:   a list of paths to input FASTA files
    :param outdir:      path to output directory
    :param fastani_exe: location of the fastANI binary
    :param fragLen:     fragment length to use
    :param kmerSize:        kmer size to use
    :param minFraction: minimum portion of the genomes that must match to trust ANI
    :param jobprefix:

    Loop over all FASTA files, generating Jobs describing fastANI command lines
    for each pairwise comparison.
    """
    fastcmds = generate_fastani_commands(
        filenames, outdir, fastani_exe, fragLen, kmerSize, minFraction
    )
    joblist = []
    for idx, fastcmd in enumerate(fastcmds):
        fastjob = pyani_jobs.Job(f"{jobprefix}_{idx:06d}", fastcmd)
        joblist.append(fastjob)
    return joblist


def generate_fastani_commands(
    filenames: List[Path],
    outdir: Path = Path("."),
    fastani_exe: Path = pyani_config.FASTANI_DEFAULT,
    fragLen: int = 3000,
    kmerSize: int = 16,
    minFraction: float = 0.2,
) -> List[str]:
    """Return list of fastANI command lines.

    :param filenames:   a list of paths to input FASTA files
    :param outdir:      path to output directory
    :param fastani_exe: location of the fastANI binary
    :param fragLen:     fragment length to use
    :param kmerSize:        kmer size to use
    :param minFraction: minimum portion of the genomes that must match to trust ANI

    Loop over all FASTA files generating fastANI command lines for each pairwise comparison.
    """
    fastani_cmdlines = []
    filenames = sorted(filenames)  # enforce ordering of filenames
    for idx, query in enumerate(filenames):
        for index, ref in enumerate(filenames):
            if idx == index:
                pass  # do we want to compare things to themselves?
            fastcmd = construct_fastani_cmdline(
                query, ref, outdir, fastani_exe, fragLen, kmerSize, minFraction
            )
            fastani_cmdlines.append(fastcmd)
    return fastani_cmdlines


# ¶ If int type is specified here, can a value of None be passed?
# ¶ I have added default values because I am not sure.
# ¶ This currently does not use the list option of fastANI
def construct_fastani_cmdline(
    query: Path,
    ref: Path,
    outdir: Path = Path("."),
    fastani_exe: Path = pyani_config.FASTANI_DEFAULT,
    fragLen: int = 3000,
    kmerSize: int = 16,
    minFraction: float = 0.2,
) -> str:
    """Will return a fastcmd item

    :param query:           path to query file
    :param ref:             path to reference file
    :param outdir:          path to output directory
    :param fastani_exe:     path to fastANI executable
    :param fragLen:         fragment length to use
    :param kmerSize:        kmer size to use
    :param minFraction:     minimum portion of the genomes that must match to trust ANI
    """
    logger = logging.getLogger(__name__)

    # Cast path strings to pathlib.Path for safety
    query, ref = Path(query), Path(ref)

    # Compile commands
    outsubdir = outdir / pyani_config.ALIGNDIR["fastANI"]
    outsubdir.mkdir(exist_ok=True)
    outfile = outsubdir / f"{query.stem}_vs_{ref.stem}.fastani"
    fastcmd = f"{fastani_exe} -q {query} -r {ref} -o {outfile} --fragLen {fragLen} -k {kmerSize} --minFraction {minFraction}"

    logger.debug("Compiled command: %s", fastcmd)

    return fastcmd


def parse_fastani_file(filename: Path) -> ComparisonResult:
    """
    Return (ref genome, query genome, ANI estimate, orthologous matches,
    sequence fragments) tuple.

    :param filename: Path, path to the input file

    Extracts the ANI estimate, the number of orthologous matches, and the
    number of sequence fragments considered from the fastANI output file.

    We assume that all fastANI comparisons are pairwise: one query and
    one reference file. The fastANI file should contain a single line.

    fsatANI *can* produce multi-line output, if a list of query/reference
    files is given to it.
    """
    # ¶ Example code from a different project
    # def add_snp(holder, type, key, *value):
    #    holder[key] = type(*value)
    # Create some sort of holder:
    # ¶ The following is for an output file with multiple lines
    # results = []
    # for line in [_.strip().split() for _ in open(filename, "r").readlines()]:
    #     if len(line) == 5:
    #         # Convert types from string to numeric
    #         line[2] = float(line[2]) / 100  # ANI value
    #         line[3] = int(line[3])  # number of matching fragments
    #         line[4] = int(line[4])  # total number of fragments
    #         results.append(ComparisonResult(*line))
    #     else:
    #         raise ValueError(f"Line contains too many/too few items: {line}")
    #         continue
    # return results
    line = open(filename, "r").readline().strip().split()
    if not line:  # No file content; either run failed or no detectable similarity
        raise PyaniFastANIException(f"Input file {filename} is empty")
    return ComparisonResult(
        line[0], line[1], 0.01 * float(line[2]), int(line[3]), int(line[4])
    )


def process_files(outdir: Path, org_lengths: Dict) -> ANIResults:
    """Return tuple of fastANI results for files in passed directory.

    :param outdir:  Path, path to the directory containing output files
    :param org_lengths: dictionary of total sequence lengths, keyed by sequence

    Returns the following pandas dataframes in an ANIResults object;
    query sequences are rows, reference sequences are columns:

    - alignment_lengths - asymmetrical: total length of alignment
    - percentage_identity - asymmetrical: percentage identity of alignment
    - alignment_coverage - asymmetrical: coverage of query and reference
    - similarity_errors - asymmetrical: count of similarity errors

    May throw a ZeroDivisionError if one or more fastANI runs failed, or a
    very distant sequence was included in the analysis.
    """
    logger = logging.getLogger(__name__)

    # Process directory to identify input files
    outfiles = pyani_files.get_input_files(outdir, ".out")

    # Hold data in ANIResults object
    results = ANIResults(list(org_lengths.keys()), "fastANI")

    # Fill diagonal NA values for alignment_length with org_lengths
    # ¶ Is this necessary? Or can ANIm not do org X org comparisons? Will the ANIResults object always have NAs on the diagonal?
    for org, length in list(org_lengths.items()):
        results.alignment_lengths[org][org] = length

    # Process .out files assuming that the filename format holds:
    # org1_vs_org2.out
    for outfile in outfiles:
        qname, rname = outfile.stem.split("_vs_")

        # We may have .out files from other analyses int eh same directory
        # If this occurs, we raise a warning and skip the .out file
        if qname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Query name %s not in input sequence list, skipping %s",
                    qname,
                    outfile,
                )
            continue
        if rname not in list(org_lengths.keys()):
            if logger:
                logger.warning(
                    "Reference name %s not in input sequence list, skipping %s",
                    rname,
                    outfile,
                )
            continue
        resultList = parse_fastani_file(
            outfile
        )  # Returns a list of ComparisonResult objects
        for result in resultList:
            if result.matches == 0 and logger is not None:
                if logger:
                    logger.warning(
                        "Total alignment length reported in %s is zero!", outfile
                    )
            tot_length = result.matches
            sim_errors = result.fragments - tot_length
            query_cover = float(tot_length) / org_lengths[qname]

            # Calculate percentage ID of aligned length. This may fail if
            # total length is zero.
            # The ZeroDivisionError that would arise should be handled
            # Common causes are that a NUCmer run failed, or that a very
            # distant sequence was included in the analysis.
            try:
                perc_id = 1 - float(sim_errors) / tot_length
            except ZeroDivisionError:
                perc_id = 0  # set arbitrary value of zero identity
                results.zero_error = True
            results.add_tot_length(qname, rname, tot_length)
            results.add_sim_errors(qname, rname, sim_errors)
            results.add_pid(qname, rname, perc_id)
            results.add_coverage(qname, rname, query_cover, None)
    return results


# """
# class ComparisonResult(NamedTuple):
#    reference: Path
#    query: Path
#    ani: float
#    matches: int
#    fragments: int
# fastANI is a fast alignment-free implementation for computing whole-genome
# Average Nucleotide Identity (ANI) between genomes
# -----------------
# Example usage:
# $ fastANI -q genome1.fa -r genome2.fa -o output.txt
# $ fastANI -q genome1.fa --rl genome_list.txt -o output.txt
#
# Available options
# -----------------
# -h, --help
#    Print this help page
#
# -r <value>, --ref <value>
#    reference genome (fasta/fastq)[.gz]
#
# --refList <value>, --rl <value>
#    a file containing list of reference genome files, one genome per line
#
# -q <value>, --query <value>
#    query genome (fasta/fastq)[.gz]
#
# --ql <value>, --queryList <value>
#    a file containing list of query genome files, one genome per line
#
# -k <value>, --kmer <value>
#    kmer size <= 16 [default : 16]
#
# -t <value>, --threads <value>
#    thread count for parallel execution [default : 1]
#
# --fragLen <value>
#    fragment length [default : 3,000]
#
# --minFraction <value>
#    minimum fraction of genome that must be shared for trusting ANI. If
#    reference and query genome size differ, smaller one among the two is
#    considered. [default : 0.2]
#
# --visualize
#    output mappings for visualization, can be enabled for single genome to
#    single genome comparison only [disabled by default]
#
# --matrix
#    also output ANI values as lower triangular matrix (format inspired from
#    phylip). If enabled, you should expect an output file with .matrix
#    extension [disabled by default]
#
# -o <value>, --output <value> [required]
#    output file name
#
# -v, --version
#    Show version
#
# """
