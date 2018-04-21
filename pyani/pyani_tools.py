# -*- coding: utf-8 -*-
"""Code to support pyani.

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

import sys
import traceback

import pandas as pd
from . import pyani_config, pyani_db, download

from Bio import SeqIO

# EXCEPTIONS
# ==========

# General exception for scripts


class PyaniException(Exception):
    """General exception for pyani."""

    def __init__(self, msg="Error in pyani module"):
        """Instantiate class."""
        Exception.__init__(self, msg)


# Report last exception as string
def last_exception():
    """Return last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# CLASSES
# =======

# Class to hold ANI dataframe results
class ANIResults(object):
    """Holds ANI dataframe results."""

    def __init__(self, labels, mode):
        """Initialise with four empty, labelled dataframes."""
        self.alignment_lengths = pd.DataFrame(index=labels, columns=labels,
                                              dtype=float)
        self.similarity_errors = pd.DataFrame(index=labels, columns=labels,
                                              dtype=float).fillna(0)
        self.percentage_identity = pd.DataFrame(index=labels, columns=labels,
                                                dtype=float).fillna(1.0)
        self.alignment_coverage = pd.DataFrame(index=labels, columns=labels,
                                               dtype=float).fillna(1.0)
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
        stemdict = {"ANIm": pyani_config.ANIM_FILESTEMS,
                    "ANIb": pyani_config.ANIB_FILESTEMS,
                    "ANIblastall": pyani_config.ANIBLASTALL_FILESTEMS}
        return zip((self.alignment_lengths, self.percentage_identity,
                    self.alignment_coverage, self.similarity_errors,
                    self.hadamard), stemdict[self.mode])
        # return [(self.alignment_lengths, "ANIm_alignment_lengths"),
        #        (self.percentage_identity, "ANIm_percentage_identity"),
        #        (self.alignment_coverage, "ANIm_alignment_coverage"),
        #        (self.similarity_errors, "ANIm_similarity_errors"),
        #        (self.hadamard, "ANIm_hadamard")]

        
# Class to hold BLAST functions
class BLASTfunctions(object):
    """Class to hold BLAST functions."""

    def __init__(self, db_func, blastn_func):
        """Instantiate class."""
        self.db_func = db_func
        self.blastn_func = blastn_func


# Class to hold BLAST executables
class BLASTexes(object):
    """Class to hold BLAST functions."""

    def __init__(self, format_exe, blast_exe):
        """Instantiate class."""
        self.format_exe = format_exe
        self.blast_exe = blast_exe


# Class to hold/build BLAST commands
class BLASTcmds(object):
    """Class for construction of BLASTN and database formatting commands."""

    def __init__(self, funcs, exes, prefix, outdir):
        """Instantiate class."""
        self.funcs = funcs
        self.exes = exes
        self.prefix = prefix
        self.outdir = outdir

    def build_db_cmd(self, fname):
        """Return database format/build command."""
        return self.funcs.db_func(fname, self.outdir,
                                  self.exes.format_exe)[0]

    def get_db_name(self, fname):
        """Return database filename."""
        return self.funcs.db_func(fname, self.outdir,
                                  self.exes.format_exe)[1]

    def build_blast_cmd(self, fname, dbname):
        """Return BLASTN command."""
        return self.funcs.blastn_func(fname, dbname, self.outdir,
                                      self.exes.blast_exe)


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
        with open(filename, 'rU') as ifh:
            count = 0
            for line in ifh.readlines():
                count += 1
                try:
                    key, label = line.strip().split('\t')
                except ValueError:
                    if logger:
                        logger.warning("Problem with class file: %s",
                                       filename)
                        logger.warning("%d: %s", (count, line.strip()))
                        logger.warning("(skipping line)")
                    continue
                else:
                    labeldict[key] = label
    return labeldict


# Return the total length of sequences in a passed FASTA file
def get_genome_length(filename):
    """Return total length of all sequences in a FASTA file."""
    with open(filename, 'r') as ifh:
        return sum([len(record) for record in SeqIO.parse(ifh, 'fasta')])


# Add the contents of a labels file to the pyani database for a given run
def add_dblabels(dbpath, run_id, labelspath):
    """Add the contents of a labels file to the pyani database.

    - dbpath       path to the pyani database
    - run_id       the ID of this run (for the database)
    - labelspath   path to the file with labels inforamtion

    The expected format of the labels file is: <HASH>\t<FILESTEM>\t<LABEL>,
    where <HASH> is the MD5 hash of the genome data (this is not checked);
    <FILESTEM> is the path to the genome file (this is intended to be a
    record for humans to audit, it's not needed for the database interaction;
    and <LABEL> is the label associated with that genome.

    If the genome referred to in the label file is not present in the database,
    that line is skipped.

    Returns a list of IDs for each label
    """
    label_ids = []
    with open(labelspath, 'r') as lfh:
        for line in lfh.readlines():
            hash, stem, label = line.strip().split('\t')
            try:
                genome_id = pyani_db.get_genome(dbpath, hash)[0][0]
            except IndexError:
                continue
            label_ids.append(pyani_db.add_genome_label(dbpath, genome_id,
                                                       run_id, label))
    return label_ids


# Add the contents of a classes file to the pyani database for a given run
def add_dbclasses(dbpath, run_id, classespath):
    """Add the contents of a classes file to the pyani database.

    - dbpath       path to the pyani database
    - run_id       the ID of this run (for the database)
    - classespath  path to the file with classes inforamtion

    The expected format of the classes file is: <HASH>\t<FILESTEM>\t<CLASS>,
    where <HASH> is the MD5 hash of the genome data (this is not checked);
    <FILESTEM> is the path to the genome file (this is intended to be a
    record for humans to audit, it's not needed for the database interaction;
    and <CLASS> is the class associated with that genome.

    If the genome referred to in the label file is not present in the database,
    that line is skipped.

    Returns a list of IDs for each class
    """
    class_ids = []
    with open(classespath, 'r') as lfh:
        for line in lfh.readlines():
            hash, stem, gclass = line.strip().split('\t')
            try:
                genome_id = pyani_db.get_genome(dbpath, hash)[0][0]
            except IndexError:
                continue
            class_ids.append(pyani_db.add_genome_class(dbpath, genome_id,
                                                       run_id, gclass))
    return class_ids
