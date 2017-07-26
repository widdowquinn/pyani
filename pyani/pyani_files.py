# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to help handle files for average nucleotide identity calculations."""

import os

from Bio import SeqIO

from . import pyani_tools


# Get a list of FASTA files from the input directory
def get_fasta_files(dirname=None):
    """Returns a list of FASTA files in the passed directory
    - dirname - path to input directory
    """
    if dirname is None:
        dirname = '.'
    infiles = get_input_files(dirname, '.fasta', '.fas', '.fa', '.fna',
                              '.fsa_nt')
    return infiles


# Return a list of paths to FASTA files in a directory
def get_fasta_paths(dirname, extlist=['.fna', '.fa', '.fasta', '.fas']):
    """Returns a list of paths to files matching a list of FASTA extensions.

    Returns the full path to each file.
    """
    return [os.path.join(dirname, fname) for fname in os.listdir(dirname) if
            os.path.isfile(os.path.join(dirname, fname)) and
            os.path.splitext(fname)[-1] in extlist]


# Get a list of FASTA files and corresponding hashes from the input directory
def get_fasta_and_hash_paths(dirname='.'):
    """Returns a list of (FASTA file, hash file) tuples in passed directory

    - dirname             - path to input directory

    Raises an IOError if the corresponding hash for a FASTA file does not exist
    """
    infiles = get_fasta_paths(dirname)
    outfiles = []
    for infile in infiles:
        hashfile = os.path.splitext(infile)[0] + '.md5'
        if not os.path.isfile(hashfile):
            raise IOError("Hashfile %s does not exist" % hashfile)
        outfiles.append((infile, hashfile))
    return outfiles


# Get list of FASTA files in a directory
def get_input_files(dirname, *ext):
    """Returns files in passed directory, filtered by extension.

    - dirname - path to input directory
    - *ext - list of arguments describing permitted file extensions
    """
    filelist = [f for f in os.listdir(dirname) if
                os.path.splitext(f)[-1] in ext]
    return [os.path.join(dirname, f) for f in filelist]


# Get lengths of input sequences
def get_sequence_lengths(fastafilenames):
    """Returns dictionary of sequence lengths, keyed by organism.

    Biopython's SeqIO module is used to parse all sequences in the FASTA
    file corresponding to each organism, and the total base count in each
    is obtained.

    NOTE: ambiguity symbols are not discounted.
    """
    tot_lengths = {}
    for fn in fastafilenames:
        tot_lengths[os.path.splitext(os.path.split(fn)[-1])[0]] = \
            sum([len(s) for s in SeqIO.parse(fn, 'fasta')])
    return tot_lengths


# Get hash string from hash file
def read_hash_string(filename):
    """Returns the hash and file strings from the passed hash file."""
    with open(filename, 'r') as ifh:
        data = ifh.read().strip().split()

    # We expect the first string in the file to be the hash (the second is the
    # filename)
    return (data[0], data[1])


# Get description string from FASTA file
def read_fasta_description(filename):
    """Returns the first description string from a FASTA file."""
    for data in SeqIO.parse(filename, 'fasta'):
        if data.description:
            return data.description
