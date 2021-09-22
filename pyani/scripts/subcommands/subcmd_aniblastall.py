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
"""Provides the aniblastall subcommand for pyani."""

import datetime
import json
import logging
import os
import subprocess

from argparse import Namespace
from logging import Logger
from itertools import permutations
from pathlib import Path
from typing import List, NamedTuple, Tuple, Dict

from Bio import SeqIO
from tqdm import tqdm

from pyani import (
    PyaniException,
    aniblastall,
    pyani_config,
    pyani_jobs,
    run_sge,
    run_multiprocessing as run_mp,
)
from pyani.pyani_files import collect_existing_output
from pyani.pyani_orm import (
    add_run,
    add_run_genomes,
    add_blastdb,
    Comparison,
    filter_existing_comparisons,
    get_session,
    PyaniORMException,
    update_comparison_matrices,
)
from pyani.pyani_tools import termcolor


# Convenience struct describing a pairwise comparison job for the SQLAlchemy
# implementation
class ComparisonJob(NamedTuple):

    """Pairwise comparison job for the SQLAlchemy implementation."""

    query: str
    subject: str
    blastcmd: str
    outfile: Path
    fragsize: int
    job: pyani_jobs.Job


def subcmd_aniblastall(args: Namespace) -> None:
    """Perform ANIblastall on all genome files in an input directory.

    :param args:  Namespace, command-line arguments
    :param logger:

    Finds ANI by the ANIblastall method, as described in ... some paper.

    All FASTA format files (selected by suffix) in the input directory are fragmented into (by default 1020nt) consecutive sections, and a BLASTALL database constructed from the whole genome input. The BLAST+ blastall tool is then used to query each set of fragments against each BLAST+ database, in turn.

    For each query, the BLAST+ .tab output is parsed to obtain alignment length, identity and similarity error count. Alignments below a threshold are not included in the calculation (this introduces systematic bias with respect to ANIm). The results are processed to calculate the ANI percentages, coverage, and similarity error.

    The calculated values are stored in the local SQLite3 database.
    """
    # Create logger
    logger = logging.getLogger(__name__)

    # Announce the analysis
    logger.info(termcolor("Running ANIblastall analysis", "red"))

    # Get BLASTALL version - this will be used in the database entreis
    blastall_version = aniblastall.get_version(args.blastall_exe)
    logger.info(termcolor(f"BLAST+ blastall version: {blastall_version}", "cyan"))

    # Use provided name, or make new one for this analysis
    start_time = datetime.datetime.now()
    name = args.name or "_".join(["ANIblastall", start_time.isoformat()])
    logger.info(termcolor(f"Analysis name: {name}", "cyan"))

    # Connect to existing database (which may be "clean" or have old analyses)
    logger.debug(f"Connecting to database {args.dbpath}")
    try:
        session = get_session(args.dbpath)
    except Exception:
        logger.error(
            f"Could not connect to database {args.dbpath} (exiting)", exc_info=True
        )
        raise SystemExit(1)

    # Add information about this run to the database
    logger.debug(f"Adding run info to database {args.dbpath}...")
    try:
        run = add_run(
            session,
            method="ANIblastall",
            cmdline=args.cmdline,
            date=start_time,
            status="started",
            name=name,
        )
    except PyaniORMException:
        logger.error("Could not add run to the database (exiting)", exc_info=True)
        raise SystemExit(1)
    logger.debug(f"\t...added run ID: {run} to the database")

    # Identify input files for comparison, and populate the database
    logger.debug(f"Adding files for {run} to database...")
    try:
        genome_ids = add_run_genomes(
            session, run, args.indir, args.classes, args.labels
        )
    except PyaniORMException:
        logger.error(
            f"Could not add genomes to database for run {run} (exiting)", exc_info=True
        )
        raise SystemExit(1)
    logger.debug(f"\t...added gnome IDs: {genome_ids}")

    # Get list of genomes for this anlaysis from the database
    logger.info("Compiling genomes for comparison")
    genomes = run.genomes.all()
    logger.debug(f"\tCollected {len(genomes)} genomes for this run")

    # Create output directories. We create the amin parent directory (args.outdir), but
    # also subdirectories for the BLAST databases.
    logger.debug(f"Creating output directory {args.outdir}")
    try:
        os.makedirs(args.outdir, exist_ok=True)
    except IOError:
        logger.error(
            f"Could not create output directory {args.outdir} (exiting)", exc_info=True
        )
        raise SystemError(1)
    fragdir = Path(str(args.outdir)) / "fragments"
    blastdbdir = Path(str(args.outdir)) / "blastalldbs"
    logger.debug("\t...creating subdirectories")
    os.makedirs(fragdir, exist_ok=True)
    os.makedirs(blastdbdir, exist_ok=True)

    # Create a new sequence fragment file and a new BLAST+ database for each input genome,
    # and add this data to the database as a row in BlastDB
    logger.info("Creating input sequence fragment files")
    fragfiles = {}
    fraglens = {}
    for genome in genomes:
        fragpath, fragsizes = fragment_fasta_file(
            Path(str(genome.path)), Path(str(fragdir)), args.fragsize
        )
        logger.info(f"fragsizes: {type(fragsizes)}")
        fragfiles.update({Path(genome.path).stem: fragpath})
        fraglens.update({Path(genome.path).stem: fragsizes})

        logger.info("Constructing formatdb command")
        dbcmd, blastdbpath = aniblastall.construct_formatdb_cmd(
            Path(genome.path), blastdbdir
        )

        logger.info("Running subprocess")
        subprocess.run(
            dbcmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

        logger.info("Adding blastdb")
        add_blastdb(
            session, genome, run, fragpath, blastdbpath, json.dumps(fragsizes), dbcmd
        )


def fragment_fasta_file(inpath: Path, outdir: Path, fragsize: int) -> Tuple[Path, str]:
    """Return path to fragmented sequence file and JSON of fragment lengths.

    :param inpath:  Path to genome file
    :param outdir:  Path to directory to hold fragmented files
    :param fragsize:  size of genome fragments

    Returns a tuple of ``(path, json)`` where ``path`` is the path to the fragment
    file and ``json`` is a JSON-ified dictionary of fragment lengths, keyed by
    fragment sequence ID.
    """
    # Generate fragments for the input sequence, looping over each contig/
    # chromosome in the input file and breaking into sections of length
    # fragsize
    sizedict = {}
    outseqs = []
    count = 0
    for seq in SeqIO.parse(inpath, "fasta"):
        idx = 0
        while idx < len(seq):
            count += 1
            newseq = seq[idx : idx + fragsize]
            newseq.id = f"frag{count:05d}"
            outseqs.append(newseq)
            sizedict[newseq.id] = len(newseq)
            idx += fragsize

    # Write fragments to output file
    fragpath = outdir / f"{inpath.stem}-fragments.fna"
    SeqIO.write(outseqs, fragpath, "fasta")
    return fragpath, sizedict
