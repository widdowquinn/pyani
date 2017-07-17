# Copyright 2016, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to support pyani."""

import logging
import os
import re
import shutil
import sys
import time
import traceback

import pandas as pd
from . import pyani_config, download


# EXCEPTIONS
#============

# General exception for scripts
class PyaniException(Exception):
    """General exception for pyani.py script"""
    def __init__(self, msg="Error in pyani.py script"):
        Exception.__init__(self, msg)

        
# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# CLASSES
#=========

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
        #return [(self.alignment_lengths, "ANIm_alignment_lengths"),
        #        (self.percentage_identity, "ANIm_percentage_identity"),
        #        (self.alignment_coverage, "ANIm_alignment_coverage"),
        #        (self.similarity_errors, "ANIm_similarity_errors"),
        #        (self.hadamard, "ANIm_hadamard")]


# Class to hold BLAST functions
class BLASTfunctions(object):
    """Class to hold BLAST functions."""
    def __init__(self, db_func, blastn_func):
        self.db_func = db_func
        self.blastn_func = blastn_func


# Class to hold BLAST executables
class BLASTexes(object):
    """Class to hold BLAST functions."""
    def __init__(self, format_exe, blast_exe):
        self.format_exe = format_exe
        self.blast_exe = blast_exe


# Class to hold/build BLAST commands
class BLASTcmds(object):
    """Class to hold BLAST command data for construction of BLASTN and
    database formatting commands.
    """
    def __init__(self, funcs, exes, prefix, outdir):
        self.funcs = funcs
        self.exes = exes
        self.prefix = prefix
        self.outdir = outdir

    def build_db_cmd(self, fname):
        """Return database format/build command"""
        return self.funcs.db_func(fname, self.outdir,
                                  self.exes.format_exe)[0]

    def get_db_name(self, fname):
        """Return database filename"""
        return self.funcs.db_func(fname, self.outdir,
                                  self.exes.format_exe)[1]

    def build_blast_cmd(self, fname, dbname):
        """Return BLASTN command"""
        return self.funcs.blastn_func(fname, dbname, self.outdir,
                                      self.exes.blast_exe)


# UTILITY FUNCTIONS
#===================   

# Create a directory (handling force/noclobber options)
def make_outdir(outdir, force, noclobber, logger):
    """Create output directory (allows for force and noclobber).

    The intended outcomes are:
    outdir doesn't exist: create outdir
    outdir exists: raise exception
    outdir exists, --force only: remove the directory tree
    outdir exists, --force --noclobber: continue with existing directory tree
                                        but do not overwrite files

    So long as the outdir is created with this function, we need only check
    for args.noclobber elsewhere to see how to proceed when a file exists.
    """
    if os.path.isdir(outdir):
        logger.warning("Output directory %s exists", outdir)
        if not force:
            raise PyaniException("Will not modify existing directory %s" %
                                 outdir)
        elif force and not noclobber:
            # Delete old directory and start again
            logger.warning("Overwrite forced. Removing %s and everything " +
                           "below it (--force)", outdir)
            shutil.rmtree(outdir)
        else:
            logger.warning("Keeping existing directory, skipping existing " +
                           "files (--force --noclobber).")
    os.makedirs(outdir, exist_ok=True)


# Make a dictionary of assembly download info
def make_asm_dict(taxon_ids, retries):
    """Return a dict of assembly UIDs, keyed by each passed taxon ID."""
    asm_dict = dict()

    for tid in taxon_ids:
        asm_uids = download.get_asm_uids(tid, retries)
        asm_dict[tid] = asm_uids.asm_ids

    return asm_dict


# Download the genome and MD5 hash from NCBI
def download_genome_and_hash(filestem, suffix, ftpstem, outdir, timeout,
                             logger):
    """Download genome and accompanying MD5 hash from NCBI.

    This function tries the (assumed to be passed) RefSeq FTP URL first and,
    if that fails, then attempts to download the corresponding GenBank data.

    We attempt to gracefully skip genomes with download errors.
    """
    # First attempt: RefSeq download
    dlstatus = download.retrieve_genome_and_hash(filestem, suffix,
                                                 ftpstem, outdir, timeout)
    if dlstatus.error is not None:  # Something went awry
        logger.warning("RefSeq download failed: skipping!\n%s", dlstatus.error)
        # Second attempt: GenBank download
        logger.warning("Trying GenBank alternative assembly")
        gbfilestem = re.sub('^GCF_', 'GCA_', filestem)
        logger.info("Retrieving URLs for %s", gbfilestem)
        gbdlstatus = download.retrieve_genome_and_hash(gbfilestem, suffix,
                                                       ftpstem, outdir,
                                                       timeout)
        if gbdlstatus.error:  # Something went awry again
            logger.error("GenBank download failed: skipping!\n%s",
                         gbdlstatus.error)
            dlstatus = gbdlstatus
            dlstatus.skipped = True

    return dlstatus


# Write class and label files
def write_classes_labels(classes, labels, outdir, classfname, labelfname,
                         noclobber, logger):
    """Write classes and labels files for the downloads."""
    # Write classes
    classfname = os.path.join(outdir, classfname)
    logger.info("Writing classes file to %s", classfname)
    if os.path.exists(classfname) and noclobber:
        logger.warning("Class file %s exists, not overwriting", classfname)
    else:
        with open(classfname, "w") as cfh:
            cfh.write('\n'.join(classes) + '\n')

    # Write labels
    labelfname = os.path.join(outdir, labelfname)
    logger.info("Writing labels file to %s", labelfname) 
    if os.path.exists(labelfname) and noclobber:
        logger.warning("Label file %s exists, not overwriting", labelfname)
    else:
        with open(labelfname, "w") as lfh:
            lfh.write('\n'.join(labels) + '\n')     

            
# Read sequence annotations in from label file
def get_labels(filename, logger=None):
    """Returns a dictionary of alternative sequence labels, or None

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
