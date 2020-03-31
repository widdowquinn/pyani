#!/usr/bin/env python
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2013-2019
# (c) The University of Strathclude 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute of Pharmaceutical and Biomedical Sciences
# The University of Strathclyde
#  Cathedral Street
# Glasgow
#  G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2018 The James Hutton Institute
# (c) The University of Strathclude 2019
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
"""Code to support pyani."""

import shutil

from logging import Logger
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Tuple,
)

import pandas as pd  # type: ignore

from Bio import SeqIO  # type: ignore

from . import pyani_config, download


class MatrixData(NamedTuple):

    """Convenience struct for matrix data returned by ORM."""

    name: str
    data: pd.DataFrame
    graphic_args: Dict


class Dependencies(NamedTuple):

    """Convenience struct for third-party dependency presence."""

    blast: Optional[str]
    legacy_blast: Optional[str]
    mummer: Optional[str]


# Class to hold ANI dataframe results
class ANIResults:

    """Holds ANI dataframe results."""

    def __init__(self, labels: List[str], mode: str) -> None:
        """Initialise with four empty, labelled dataframes.

        :param labels:
        :param mode:
        """
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

    def add_tot_length(
        self, qname: str, sname: str, value: float, sym: bool = True
    ) -> None:
        """Add a total length value to self.alignment_lengths.

        :param qname:
        :param sname:
        :param value:
        :param sym:
        """
        self.alignment_lengths.loc[qname, sname] = value
        if sym:
            self.alignment_lengths.loc[sname, qname] = value

    def add_sim_errors(
        self, qname: str, sname: str, value: float, sym: bool = True
    ) -> None:
        """Add a similarity error value to self.similarity_errors.

        :param qname:
        :param sname:
        :param value:
        :param sym:
        """
        self.similarity_errors.loc[qname, sname] = value
        if sym:
            self.similarity_errors.loc[sname, qname] = value

    def add_pid(self, qname: str, sname: str, value: float, sym: bool = True) -> None:
        """Add a percentage identity value to self.percentage_identity.

        :param qname:
        :param sname:
        :param value:
        :param sym:
        """
        self.percentage_identity.loc[qname, sname] = value
        if sym:
            self.percentage_identity.loc[sname, qname] = value

    def add_coverage(
        self, qname: str, sname: str, qcover: float, scover: Optional[float] = None
    ) -> None:
        """Add percentage coverage values to self.alignment_coverage.

        :param qname:
        :param sname:
        :param value:
        :param sym:
        """
        self.alignment_coverage.loc[qname, sname] = qcover
        if scover:
            self.alignment_coverage.loc[sname, qname] = scover

    @property
    def hadamard(self) -> float:
        """Return Hadamard matrix (identity * coverage)."""
        return self.percentage_identity * self.alignment_coverage

    @property
    def data(self) -> Iterator[Tuple[Any, str]]:
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


class BLASTfunctions(NamedTuple):

    """Convenience structure to hold BLAST functions."""

    db_func: Callable
    blastn_func: Callable


class BLASTexes(NamedTuple):

    """Convenience structure to hold BLAST executables."""

    format_exe: Path
    blast_exe: Path


# Class to hold/build BLAST commands
class BLASTcmds:

    """Class for construction of BLASTN and database formatting commands."""

    def __init__(
        self, funcs: BLASTfunctions, exes: BLASTexes, prefix: str, outdir: Path
    ) -> None:
        """Instantiate class.

        :param funcs:  BLASTfunctions, containing functions for this BLAST analysis
        :param exes:  BLASTexes, containing executables for this BLAST analysis
        :param prefix:  str, prefix for outputs from this BLAST analysis
        :param outdir:  Path to output directory for this BLAST analysis
        """
        self.funcs = funcs
        self.exes = exes
        self.prefix = prefix
        self.outdir = outdir

    def build_db_cmd(self, fname: Path) -> str:
        """Return database format/build command.

        :param fname:
        """
        return self.funcs.db_func(fname, self.outdir, self.exes.format_exe)[0]

    def get_db_name(self, fname: Path) -> str:
        """Return database filename.

        :param fname:
        """
        return self.funcs.db_func(fname, self.outdir, self.exes.format_exe)[1]

    def build_blast_cmd(self, fname: Path, dbname: Path):
        """Return BLASTN command.

        :param fname:  Path to query file
        :param dbname:  Path to database
        """
        return self.funcs.blastn_func(fname, dbname, self.outdir, self.exes.blast_exe)


# UTILITY FUNCTIONS
# =================

# Make a dictionary of assembly download info
def make_asm_dict(taxon_ids: Iterable[str], retries: int) -> Dict:
    """Return dict of assembly UIDs, keyed by each passed taxon ID.

    :param taxon_ids:  Iterable of NCBI taxonomy IDs
    :param retries:  Number of Entrez retry attempts
    """
    asm_dict = dict()

    for tid in taxon_ids:
        asm_uids = download.get_asm_uids(tid, retries)
        asm_dict[tid] = asm_uids.asm_ids

    return asm_dict


# Read sequence annotations in from file
def get_labels(filename: Path, logger: Logger = None) -> Dict:
    r"""Return dictionary of alternative sequence labels, or None.

    :param filename:  path to file containing tab-separated table of labels
    :param logger:  logging object

    Input files should be formatted as <hash>\t<key>\t<label>, one pair per line.
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
                    hash_md5, key, label = line.strip().split("\t")
                except ValueError:
                    if logger:
                        logger.warning("Problem with class file: %s", filename)
                        logger.warning("line %d: %s", count, line.strip())
                        logger.warning("(skipping line)")
                    continue
                else:
                    labeldict[key] = label
    return labeldict


# Return the total length of sequences in a passed FASTA file
def get_genome_length(filename: Path) -> int:
    """Return total length of all sequences in a FASTA file.

    :param filename:  path to FASTA file
    """
    with open(filename, "r") as ifh:
        return sum([len(record) for record in SeqIO.parse(ifh, "fasta")])


# Helper function to label results matrices from Run objects
def label_results_matrix(matrix: pd.DataFrame, labels: Dict) -> pd.DataFrame:
    """Return results matrix dataframe with labels.

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
def has_dependencies() -> Dependencies:
    """Return NamedTuple indicating if 3rd dependencies are available."""
    return Dependencies(
        shutil.which("blastn"), shutil.which("blastall"), shutil.which("nucmer")
    )
