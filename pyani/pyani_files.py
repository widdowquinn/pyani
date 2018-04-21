# -*- coding: utf-8 -*-
"""Code to help handle files for average nucleotide identity calculations.

(c) The James Hutton Institute 2013-2017
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

Copyright (c) 2013-2017 The James Hutton Institute

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

from Bio import SeqIO


# Get a list of FASTA files from the input directory
def get_fasta_files(dirname=None):
    """Return a list of FASTA files in the passed directory.

    - dirname - path to input directory
    """
    if dirname is None:
        dirname = '.'
    infiles = get_input_files(dirname, '.fasta', '.fas', '.fa', '.fna',
                              '.fsa_nt')
    return infiles


# Return a list of paths to FASTA files in a directory
def get_fasta_paths(dirname, extlist=['.fna', '.fa', '.fasta', '.fas']):
    """Return a list of paths to files matching a list of FASTA extensions.

    Returns the full path to each file.
    """
    return [os.path.join(dirname, fname) for fname in os.listdir(dirname) if
            os.path.isfile(os.path.join(dirname, fname)) and
            os.path.splitext(fname)[-1] in extlist]


# Get a list of FASTA files and corresponding hashes from the input directory
def get_fasta_and_hash_paths(dirname='.'):
    """Return a list of (FASTA file, hash file) tuples in passed directory.

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
    """Return files in passed directory, filtered by extension.

    - dirname - path to input directory
    - *ext - list of arguments describing permitted file extensions
    """
    filelist = [f for f in os.listdir(dirname) if
                os.path.splitext(f)[-1] in ext]
    return [os.path.join(dirname, f) for f in filelist]


# Get lengths of input sequences
def get_sequence_lengths(fastafilenames):
    """Return dictionary of sequence lengths, keyed by organism.

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
    """Return the hash and file strings from the passed hash file."""
    with open(filename, 'r') as ifh:
        data = ifh.read().strip().split()

    # We expect the first string in the file to be the hash (the second is the
    # filename)
    return (data[0], data[1])


# Get description string from FASTA file
def read_fasta_description(filename):
    """Return the first description string from a FASTA file."""
    for data in SeqIO.parse(filename, 'fasta'):
        if data.description:
            return data.description
