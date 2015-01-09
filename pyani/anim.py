# Copyright 2013-2015, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence. 
# Please see the LICENSE file that should have been included as part of
# this package.

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

import pyani_config
import os


# Generate list of NUCmer pairwise comparison command lines from
# passed sequence filenames
def generate_nucmer_commands(filenames,
                             nucmer_exe=pyani_config.NUCMER_DEFAULT,
                             maxmatch=False):
    """Return a list of NUCmer command-lines for ANIm

    - filenames - a list of paths to input FASTA files
    - nucmer_exe - location of the nucmer binary
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option

    Loop over all FASTA files generating NUCmer command lines for each
    pairwise comparison.
    """
    cmdlines = []
    for idx, f1 in enumerate(filenames[:-1]):
        cmdlines.extend([construct_nucmer_cmdline(f1, f2, nucmer_exe) for 
                         f2 in filenames[idx+1:]])
    return cmdlines


# Generate single NUCmer pairwise comparison command line from pair of
# input filenames
def construct_nucmer_cmdline(filename1, filename2, 
                             nucmer_exe=pyani_config.NUCMER_DEFAULT,
                             maxmatch=False):
    """Returns a single NUCmer pairwise comparison command.

    - filename1 - query FASTA filepath
    - filename2 - subject FASTA filepath
    - maxmatch - Boolean flag indicating to use NUCmer's -maxmatch option
    """
    outprefix = "%s_vs_%s" %\
                (os.path.splitext(os.path.split(filename1)[-1])[0],
                 os.path.splitext(os.path.split(filename2)[-1])[0])
    if maxmatch:
        mode = "-maxmatch"
    else:
        mode = "-mum"
    return "{0} {1} -p {2} {3} {4}".format(nucmer_exe, mode, outprefix,
                                           filename1, filename2)
    
# Parse NUCmer delta file to get total alignment length and total sim_errors
def parse_delta(filename):
    """Parses a NUCmer output .delta file

    - filename - path to the input .delta file

    Extracts the aligned length and number of similarity errors for each
    aligned uniquely-matched region, and returns the cumulative total for
    each as a tuple.
    """
    aln_length, sim_errors = 0, 0
    for line in [l.strip().split() for l in open(filename, 'rU').readlines()]:
        if line[0] == 'NUCMER' or line[0].startswith('>'):  # Skip headers
            continue
        # We only process lines with seven columns:
        if len(line) == 7:
            aln_length += abs(int(line[1]) - int(line[0]))
            sim_errors += int(line[4])
    return aln_length, sim_errors


    
