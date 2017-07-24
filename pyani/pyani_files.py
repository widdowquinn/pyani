# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to help handle files for average nucleotide identity calculations."""

import os

from . import pyani_tools


# Return a list of paths to FASTA files in a directory
def get_fasta_paths(dirname, extlist=['.fna', '.fa', '.fasta', '.fas']):
    """Returns a list of paths to files matching a list of FASTA extensions.

    Returns the full path to each file.
    """
    return [os.path.join(dirname, fname) for fname in os.listdir(dirname) if
            os.path.isfile(os.path.join(dirname, fname)) and
            os.path.splitext(fname)[-1] in extlist]


# Get a list of FASTA files and corresponding hashes from the input directory
def get_fasta_and_hash_paths(dirname='.', create_hash=True):
    """Returns a list of (FASTA file, hash file) tuples in passed directory

    - dirname             - path to input directory
    - create_hash         - Boolean; if True and the FASTA file has no
                            corresponding hash ('.md5') file, a hash file
                            is created.

    If the FASTA file does not have a corresponding hash, then the
    corresponding hash is created, when create_hash=True. Otherwise the file
    path is not returned.
    """
    infiles = pyani_tools.get_fasta_paths(dirname)
    for file in infiles:
        hashfile = os.path.splitext(infile) + '.md5'
        if not os.path.isfile(hashfile):
            raise IOError("Hashfile %s does not exist" % hashfile)
    return infiles



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
