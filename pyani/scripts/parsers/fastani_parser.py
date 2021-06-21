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

"""Provides parser for fastani subcommand."""

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, _SubParsersAction
from pathlib import Path
from typing import List, Optional

from pyani import pyani_config
from pyani.scripts import subcommands


def build(
    subps: _SubParsersAction,
    parents: Optional[List[ArgumentParser]] = None,
) -> None:
    """Return a command-line parser for the fastani subcommand.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    The terminology may be confusing, but in practice the main parser collects
    command-line arguments that are then available to this parser, which inherits
    options from the parsers in `parents` in addition to those defined below.
    """
    parser = subps.add_parser(
        "fastani", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required arguments: input and output directories
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        action="store",
        dest="indir",
        default=None,
        type=Path,
        help="input genome directory",
    )
    parser.add_argument(  # I think this should be handled with pyani's overall output
        "-o",
        "--out",
        required=True,
        action="store",
        dest="outdir",
        default=None,
        type=Path,
        help="output analysis results directory",
    )
    # Optional arguments
    parser.add_argument(  # keep ;)
        "--dbpath",
        action="store",
        dest="dbpath",
        default=Path(".pyani/pyanidb"),
        type=Path,
        help="path to pyani database",
    )
    parser.add_argument(
        "--fastani_exe",
        dest="fastani_exe",
        action="store",
        default=pyani_config.FASTANI_DEFAULT,
        type=Path,
        help="path to fastani executable",
    )
    parser.add_argument(
        "-k",
        "--kmer",
        dest="kmerSize",
        action="store",
        default=16,
        type=int,
        help="kmer size <= 16",
    )
    parser.add_argument(
        "--fragLen",
        dest="fragLen",
        action="store",
        default=3000,
        type=int,
        help="fragment length",
    )
    parser.add_argument(  # should be under control of scheduling parser?
        "-t",
        "--threads",
        dest="threads",
        action="store",
        default=1,
        type=int,
        help="thread count for parallel execution",
    )
    parser.add_argument(
        "--minFraction",
        dest="minFraction",
        action="store",
        default=0.2,
        type=float,
        help="Minimum fraction of genome that must be shared for trusting ANI. If reference and query genome size differ, smaller one among the two is considered.",
    )
    parser.add_argument(  # should be handled by the report subcommand; BUT this would be useful for concordance testing - does the extra pyani processing distort fastANI's output?
        "--matrix",
        dest="matrix",
        action="store_true",
        default=False,
        help="Also output ANI values as lower triangular matrix (format inspired from phylip). If enabled, you should expect an output file with .matrix extension.",
    )
    """parser.add_argument(
        "--visualize",
        dest="visualize",
        action="store_true",
        default=False,
        help="Output mappings for visualization, can be enabled for single genome to single genome comparison only",
    )"""
    parser.set_defaults(func=subcommands.subcmd_fastani)


# """
# -----------------
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
