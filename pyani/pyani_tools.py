#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Code to support pyani.

(c) The James Hutton Institute 2013-2019
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

Copyright (c) 2013-2019 The James Hutton Institute

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

import shutil

from collections import namedtuple
from functools import lru_cache

import pandas as pd

from Bio import SeqIO

from . import pyani_config, download

# Convenience struct for matrix data returned by ORM
MatrixData = namedtuple("MatrixData", "name data graphic_args")

# Convenience struct for third-party dependency presence
Dependencies = namedtuple("Dependencies", "blast legacy_blast mummer")


# CLASSES
# =======

# Class to hold ANI dataframe results
class ANIResults:
    """Holds ANI dataframe results."""

    def __init__(self, labels, mode):
        """Initialise with four empty, labelled dataframes."""
        self.alignment_lengths = pd.DataFrame(index=labels, columns=labels, dtype=float)
        self.similarity_errors = pd.DataFrame(
            index=labels, columns=labels, dtype=float
        ).fillna(0)
        self.percentage_identity = pd.DataFrame(
            index=labels, columns=labels, dtype=float
        ).fillna(1.0)
        self.alignment_coverage = pd.DataFrame(
            index=labels, columns=labels, dtype=float
        ).fillna(1.0)
        self.zero_error = False
        self.mode = mode

    def add_tot_length(self, qname, sname, value, sym=True):
        """Add a total length value to self.alignment_lengths."""
        self.alignment_lengths.loc[qname, sname] = value
        if sym:
            self.alignment_lengths.loc[sname, qname] = value

    def add_sim_errors(self, qname, sname, value, sym=True):
        """Add a similarity error value to self.similarity_errors."""
        self.similarity_errors.loc[qname, sname] = value
        if sym:
            self.similarity_errors.loc[sname, qname] = value

    def add_pid(self, qname, sname, value, sym=True):
        """Add a percentage identity value to self.percentage_identity."""
        self.percentage_identity.loc[qname, sname] = value
        if sym:
            self.percentage_identity.loc[sname, qname] = value

    def add_coverage(self, qname, sname, qcover, scover=None):
        """Add percentage coverage values to self.alignment_coverage."""
        self.alignment_coverage.loc[qname, sname] = qcover
        if scover:
            self.alignment_coverage.loc[sname, qname] = scover

    @property
    def hadamard(self):
        """Return Hadamard matrix (identity * coverage)."""
        return self.percentage_identity * self.alignment_coverage

    @property
    def data(self):
        """Return list of (dataframe, filestem) tuples."""
        stemdict = {
            "ANIm": pyani_config.ANIM_FILESTEMS,
            "ANIb": pyani_config.ANIB_FILESTEMS,
            "ANIblastall": pyani_config.ANIBLASTALL_FILESTEMS,
        }
        return zip(
            (
                self.alignment_lengths,
                self.percentage_identity,
                self.alignment_coverage,
                self.similarity_errors,
                self.hadamard,
            ),
            stemdict[self.mode],
        )
        # return [(self.alignment_lengths, "ANIm_alignment_lengths"),
        #        (self.percentage_identity, "ANIm_percentage_identity"),
        #        (self.alignment_coverage, "ANIm_alignment_coverage"),
        #        (self.similarity_errors, "ANIm_similarity_errors"),
        #        (self.hadamard, "ANIm_hadamard")]


# Class to hold BLAST functions
class BLASTfunctions:
    """Class to hold BLAST functions."""

    def __init__(self, db_func, blastn_func):
        """Instantiate class."""
        self.db_func = db_func
        self.blastn_func = blastn_func


# Class to hold BLAST executables
class BLASTexes:
    """Class to hold BLAST functions."""

    def __init__(self, format_exe, blast_exe):
        """Instantiate class."""
        self.format_exe = format_exe
        self.blast_exe = blast_exe


# Class to hold/build BLAST commands
class BLASTcmds:
    """Class for construction of BLASTN and database formatting commands."""

    def __init__(self, funcs, exes, prefix, outdir):
        """Instantiate class."""
        self.funcs = funcs
        self.exes = exes
        self.prefix = prefix
        self.outdir = outdir

    def build_db_cmd(self, fname):
        """Return database format/build command."""
        return self.funcs.db_func(fname, self.outdir, self.exes.format_exe)[0]

    def get_db_name(self, fname):
        """Return database filename."""
        return self.funcs.db_func(fname, self.outdir, self.exes.format_exe)[1]

    def build_blast_cmd(self, fname, dbname):
        """Return BLASTN command."""
        return self.funcs.blastn_func(fname, dbname, self.outdir, self.exes.blast_exe)


# UTILITY FUNCTIONS
# =================

# Make a dictionary of assembly download info
def make_asm_dict(taxon_ids, retries):
    """Return dict of assembly UIDs, keyed by each passed taxon ID."""
    asm_dict = dict()

    for tid in taxon_ids:
        asm_uids = download.get_asm_uids(tid, retries)
        asm_dict[tid] = asm_uids.asm_ids

    return asm_dict


# Read sequence annotations in from file
def get_labels(filename, logger=None):
    r"""Return dictionary of alternative sequence labels, or None.

    - filename - path to file containing tab-separated table of labels
    Input files should be formatted as <key>\t<label>, one pair per line.
    """
    labeldict = {}
    if filename is not None:
        if logger:
            logger.info("Reading labels from %s", filename)
        with open(filename, "r") as ifh:
            count = 0
            for line in ifh.readlines():
                count += 1
                try:
                    key, label = line.strip().split("\t")
                except ValueError:
                    if logger:
                        logger.warning("Problem with class file: %s", filename)
                        logger.warning("%d: %s", (count, line.strip()))
                        logger.warning("(skipping line)")
                    continue
                else:
                    labeldict[key] = label
    return labeldict


# Return the total length of sequences in a passed FASTA file
def get_genome_length(filename):
    """Return total length of all sequences in a FASTA file."""
    with open(filename, "r") as ifh:
        return sum([len(record) for record in SeqIO.parse(ifh, "fasta")])


# Helper function to label results matrices from Run objects
def label_results_matrix(matrix, labels):
    """Return results matrix dataframe with labels

    :param matrix:  results dataframe deriving from Run object
    :param labels:  dictionary of genome labels
        labels must be keyed by index/col values from matrix

    Applies the labels from the dictionary to the dataframe in
    matrix, and returns the result.
    """
    matrix.columns = [f"{labels.get(_, _)}:{_}" for _ in matrix.columns]
    matrix.index = [f"{labels.get(_, _)}:{_}" for _ in matrix.index]
    return matrix


# Helper function that establishes whether dependencies are present
# This caches the most recent result
def has_dependencies():
    """Returns a namedtuple indicating whether third-party dependencies
    are available.

    An LRU cache stores the last set of test results, for convenience
    """
    return Dependencies(
        shutil.which("blastn"), shutil.which("blastall"), shutil.which("nucmer")
    )
