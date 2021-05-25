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
# Copyright (c) 2021–Present university of Strathclyde
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

import platform
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
from .pyani_tools import ANIResults, BLASTcmds, BLASTexes, BLASTfunctions


class ComparisonResult(NamedTuple):
    reference: Path
    query: Path
    ani: float
    orthologs: int
    fragments: int


def get_version(fastani_exe: Path = pyani_config.FASTANI_DEFAULT) -> str:
    """Return FastANI package version as a string.

    :param fastani_exe: path to FastANI executable

    We expect fastANI to return a string on STDOUT as

    .. code-block:: bash

        $ ./fastANI -v
        version 1.32

    we concatenate this with the OS name.
    """
    cmdline = [fastani_exe, "-v"]  # type: List
    result = subprocess.run(
        cmdline, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
    )
    # type CompletedProcess
    match = re.search(r"(?<=version\s)[0-9\.]*", str(result.stderr, "utf-8"))
    version = match.group()  # type: ignore
    return f"{platform.system()}_{version}"


# Generate list of Job objects, one per fastANI run -> this is maybe not necessary
def generate_fastani_jobs(
    reference: Path,
    query: Path,
    outdir: Path = Path("."),
    fastani_exe: Path = pyani_config.FASTANI_DEFAULT,
    # kmer
):
    pass


def generate_commands():
    pass


def construct_fastani_cmdline():
    """Will return a fastcmd item"""
    # return fastcmd
    pass


def parse_output_file(filename: Path) -> ComparisonResult:
    """
    Return (ref genome, query genome, ANI estimate, orthologous matches,
    sequence fragments) tuple.

    :param filename: Path, path to the input file

    Extracts the ANI estimate, the number of orthologous matches, and the
    number of sequence fragments considered from the fastANI output file.
    """
    # ¶ Example code from a different project
    # def add_snp(holder, type, key, *value):
    #    holder[key] = type(*value)
    # Create some sort of holder:
    box = []
    for line in [_.strip().split() for _ in open(filename, "r").readlines()]:
        if len(line) == 5:
            box.append(ComparisonResult(*line))
        else:
            raise ValueError("Line contains too many/too few items: %s" % line)
            continue
    return box


def process_files():
    pass


# """
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
