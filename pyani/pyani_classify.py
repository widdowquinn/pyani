# Copyright 2016, The James Hutton Insitute
# Author: Leighton Pritchard
#
# This code is part of the pyani package, and is governed by its licence.
# Please see the LICENSE file that should have been included as part of
# this package.

"""Code to help classify input genomes according to ANI output."""

import csv

import networkx as nx
import pandas as pd


# Identify input files from input directory
def get_pyani_inputs(indir, logger=None):
    """Check files in the input directory, and return the names of the %id
    and %cov output files from the analysis.
    """
    


# Load pyani output into dataframes
def load_pyani_data(covfile, idfile, labelfile=None, logger=None):
    """Load coverage and identity matrices from pyani output, and an optional
    correponding labels.txt file for the sequences.

    Returns a tuple of (coverage data (DataFrame), ID data (DataFrame), and
    names for genomes in the data (dict))
    """
    # Load pyani outputs as DataFrames
    covdata = pd.DataFrame.from_csv(covfile, sep='\t')
    iddata = pd.DataFrame.from_csv(idfile, sep='\t')

    # Generate dictionary of genome labels. Use the labels.txt info if
    # available, and the names from the coverage data if not.
    if labelfile:
        with open(labels, mode='r') as fh:
            reader = csv.reader(fh, delimiter='\t')
            labeldict = {rows[0]:rows[1] for rows in reader
                         if rows[0] in covdata.columns}
    else:
        labeldict = {name:name for name in covdata.columns}

    return(covdata, iddata, labeldict)
