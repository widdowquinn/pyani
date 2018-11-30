#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""subcmd_anib.py

Provides the anib subcommand for pyani

(c) The James Hutton Institute 2017-18

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD6 9LH,
Scotland,
UK

The MIT License

Copyright (c) 2017-18 The James Hutton Institute

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

import datetime

from pyani import anib, pyani_orm
from pyani.pyani_orm import Run, add_run_genomes
from pyani.pyani_tools import last_exception


def subcmd_anib(args, logger):
    """Perform ANIb on all genome files in an input directory.

    Finds ANI using the ANIb method as described in Richter et al. (2009)
    Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

    All input genome sequences are fragmented into 1000nt sections. These are
    compiled into BLAST+ nucleotide databases (one per genome), and for each
    pairwise genome comparison between genomes A and B, the fragments from
    genome A are queried against the database for genome B; then the fragments
    from genome B are queried against the database for genome A. All queries
    are performed with BLASTN+

    The BLASTN+ output files are processed for each comparison to calculate,
    for the 'best' matches to each query: total alignment length, average
    percentage identity, coverage, 'similarity errors' (unaligned bases).

    The calculated values are deposited in the SQLite3 database used for the
    current analysis.
    """
    # Announce the analysis
    logger.info("Running ANIb analysis")

    # Get current BLAST version
    blast_version = anib.get_version(args.blastn_exe)
    logger.info("Current BLASTN+ version: %s", blast_version)

    # Use the provided name or construct one for this analysis
    start_time = datetime.datetime.now()
    name = args.name or "{}_{}".format("ANIb", start_time.isoformat())
    logger.info("Analysis name: %s", name)

    # Add analysis info to the database
    # 1. Connect to the active database
    logger.info("Connecting to database: %s", args.dbpath)
    try:
        session = pyani_orm.get_session(args.dbpath)
    except Exception:
        logger.error("Could not connect to database %s (exiting)", args.dbpath)
        logger.error(last_exception())
        raise SystemExit(1)

    # 2. Add run info to database
    logger.info("Adding run information to database %s", args.dbpath)
    run = Run(
        method="ANIb",
        cmdline=args.cmdline,
        date=start_time,
        status="started",
        name=name,
    )
    try:
        session.add(run)
        session.commit()
        logger.info("Added run %s to the database", run)
    except Exception:
        logger.error("Could not add run to the database (exiting)")
        logger.error(last_exception())
        raise SystemExit(1)

    # 3. Identify input files for comparison and populate the database
    logger.info("Adding genomes for the run to the database")
    logger.info("\tInput directory: %s", args.indir)
    logger.info("\tClasses file: %s", args.classes)
    logger.info("\tLabels file: %s", args.labels)
    try:
        add_run_genomes(session, run, args.indir, args.classes, args.labels)
    except:
        logger.error("Could not add genomes to database (exiting)")
        logger.error(last_exception())
        raise SystemExit(1)
