#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
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
# Copyright (c) 2017-2019 The James Hutton Institute
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
"""Provides the index subcommand for pyani."""

import os

from Bio import SeqIO

from pyani import download, pyani_files


def subcmd_index(args, logger):
    """Generate a file with the MD5 hash for each genome in an input directory.

    :param args:  Namespace, received command-line arguments
    :param logger:  logging object

    Identify the genome files in the input directory, and generate a single
    MD5 for each so that <genome>.fna produces <genome>.md5

    Genome files (FASTA) are identified from the file extension.
    """
    # Get list of FASTA files in the input directory
    logger.info("Scanning directory %s for FASTA files", args.indir)
    fpaths = pyani_files.get_fasta_paths(args.indir)
    logger.info("Found FASTA files:")
    logger.info([f"\t{fpath}\n" for fpath in fpaths])

    # Lists of class/label information
    classes = []
    labels = []

    # Create MD5 hash for each file, if needed
    for fpath in fpaths:
        hashfname = fpath.with_suffix(".md5")
        if os.path.isfile(hashfname):
            logger.info("%s already indexed (using existing hash)", fpath)
            with open(hashfname, "r") as ifh:
                datahash = ifh.readline().split()[0]
        else:
            # Write an .md5 hash file
            datahash = download.create_hash(fpath)
            logger.info("Writing hash to %s", hashfname)
            with open(hashfname, "w") as hfh:
                hfh.write(f"{datahash}\t{fpath}\n")

        # Parse the file and get the label/class information
        with open(fpath, "r") as sfh:
            label = list(SeqIO.parse(sfh, "fasta"))[0].description.split(" ", 1)[-1]
        labels.append("\t".join([datahash, fpath.stem, label]))
        classes.append("\t".join([datahash, fpath.stem, label]))

    # Write class and label files
    classfname = os.path.join(args.indir, args.classfname)
    logger.info("Writing classes file to %s", classfname)
    if os.path.exists(classfname):
        logger.warning("Class file %s exists, not overwriting", classfname)
    else:
        with open(classfname, "w") as ofh:
            ofh.write("\n".join(classes) + "\n")

    labelfname = os.path.join(args.indir, args.labelfname)
    logger.info("Writing labels file to %s", labelfname)
    if os.path.exists(labelfname):
        logger.warning("Labels file %s exists, not overwriting", labelfname)
    else:
        with open(labelfname, "w") as ofh:
            ofh.write("\n".join(labels) + "\n")
