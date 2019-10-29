# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2013-2019
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
# Copyright (c) 2013-2019 The James Hutton Institute
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
"""Code to handle files for average nucleotide identity calculations."""

import os

from Bio import SeqIO

from pyani import PyaniException


# General exception for scripts
class PyaniFilesException(PyaniException):

    """Exception raised by pyani when file interaction goes bad."""


# Get a list of FASTA files from the input directory
def get_fasta_files(dirname=None):
    """Return a list of FASTA files in the passed directory.

    :param dirname: - path to input directory
    """
    if dirname is None:
        dirname = "."
    infiles = get_input_files(dirname, ".fasta", ".fas", ".fa", ".fna", ".fsa_nt")
    return infiles


# Return a list of paths to FASTA files in a directory
def get_fasta_paths(dirname, extlist=None):
    """Return a list of paths to files matching a list of FASTA extensions.

    :param dirname:
    :param extlist:

    Returns the full path to each file.
    """
    extlist = extlist or [".fna", ".fa", ".fasta", ".fas"]
    return [
        os.path.join(dirname, fname)
        for fname in os.listdir(dirname)
        if os.path.isfile(os.path.join(dirname, fname))
        and os.path.splitext(fname)[-1] in extlist
    ]


# Get a list of FASTA files and corresponding hashes from the input directory
def get_fasta_and_hash_paths(dirname="."):
    """Return a list of (FASTA file, hash file) tuples in passed directory.

    :param dirname:  str, path to input directory

    Raises an IOError if the corresponding hash for a FASTA file does not exist
    """
    infiles = get_fasta_paths(dirname)
    outfiles = []
    for infile in infiles:
        hashfile = os.path.splitext(infile)[0] + ".md5"
        if not os.path.isfile(hashfile):
            raise IOError("Hashfile %s does not exist" % hashfile)
        outfiles.append((infile, hashfile))
    return outfiles


# Get list of FASTA files in a directory
def get_input_files(dirname, *ext):
    """Return files in passed directory, filtered by extension.

    :param dirname:  str, path to input directory
    :param *ext:  list of arguments describing permitted file extensions
    """
    filelist = [f for f in os.listdir(dirname) if os.path.splitext(f)[-1] in ext]
    return [os.path.join(dirname, f) for f in filelist]


# Get lengths of input sequences
def get_sequence_lengths(fastafilenames):
    """Return dictionary of sequence lengths, keyed by organism.

    :param fastafilenames:  str, path to input FASTA file

    Biopython's SeqIO module is used to parse all sequences in the FASTA
    file corresponding to each organism, and the total base count in each
    is obtained.

    NOTE: ambiguity symbols are not discounted.
    """
    tot_lengths = {}
    for fname in fastafilenames:
        tot_lengths[os.path.splitext(os.path.split(fname)[-1])[0]] = sum(
            [len(s) for s in SeqIO.parse(fname, "fasta")]
        )
    return tot_lengths


# Get hash string from hash file
def read_hash_string(filename):
    """Return the hash and file strings from the passed hash file.

    :param filename:  str, path to file containing hash string
    """
    try:
        with open(filename, "r") as ifh:
            data = ifh.read().strip().split()
    except Exception:
        raise PyaniFilesException(f"Could not load hash file {filename}.")

    # We expect the first string in the file to be the hash (the second is the
    # filename)
    return (data[0], data[1])


# Get description string from FASTA file
def read_fasta_description(filename):
    """Return the first description string from a FASTA file.

    :param filename:  str, path to FASTA file
    """
    for data in SeqIO.parse(filename, "fasta"):
        if data.description:
            return data.description
    raise PyaniFilesException(f"No sequences in {filename} contain a description.")


# Load class or label file as dictionary
def load_classes_labels(path):
    r"""Return a dictionary of genome classes or labels keyed by hash.

    :param path:  str, path to classes or labels file

    The expected format of the classes and labels files is:
    <HASH>\t<FILESTEM>\t<CLASS>|<LABEL>,
    where <HASH> is the MD5 hash of the genome data (this is not checked);
    <FILESTEM> is the path to the genome file (this is intended to be a
    record for humans to audit, it's not needed for the database interaction;
    and <CLASS>|<LABEL> is the class or label associated with that genome.
    """
    datadict = {}
    with open(path, "r") as ifh:
        for line in ifh.readlines():
            genomehash, _, data = line.strip().split("\t")
            datadict[genomehash] = data
    return datadict


# Collect existing output files when in recovery mode
def collect_existing_output(dirpath, program, args):
    """Return a list of existing output files at dirpath.

    :param dirpath:  str, path to existing output directory
    :param args:  Namespace, command-line arguments for the run
    """
    # Obtain collection of expected output files already present in directory
    if program == "nucmer":
        if args.nofilter:
            suffix = ".delta"
        else:
            suffix = ".filter"
    elif program == "blastn":
        suffix = ".blast_tab"
    existingfiles = [
        fname for fname in os.listdir(dirpath) if os.path.splitext(fname)[-1] == suffix
    ]
    return existingfiles
