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
import os

from itertools import permutations

from Bio import SeqIO
from tqdm import tqdm

from pyani import (
    anib,
    last_exception,
    pyani_config,
    run_multiprocessing as run_mp,
    run_sge,
)
from pyani.pyani_files import (
    collect_existing_output,
    get_fasta_and_hash_paths,
    read_hash_string,
)
from pyani.pyani_orm import (
    BlastDB,
    Comparison,
    Genome,
    Run,
    add_run_genomes,
    filter_existing_comparisons,
    get_comparison_dict,
    get_session,
    update_comparison_matrices,
)


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
    # 1. connect to the active database
    logger.info("Connecting to database: %s", args.dbpath)
    try:
        session = get_session(args.dbpath)
    except Exception:
        logger.error("Could not connect to database %s (exiting)", args.dbpath)
        logger.error(last_exception())
        raise SystemExit(1)

    # 2. add run info to database
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

    # 3. identify input files for comparison and populate the database
    logger.info("Adding genomes for the run to the database")
    logger.info("\tInput directory: %s", args.indir)
    logger.info("\tClasses file: %s", args.classes)
    logger.info("\tLabels file: %s", args.labels)
    try:
        add_run_genomes(session, run, args.indir, args.classes, args.labels)
    except Exception:
        logger.error("Could not add genomes to database (exiting)")
        logger.error(last_exception())
        raise SystemExit(1)

    # Construct command lines for BLAST nucleotide database construction and
    # BLASTN+ comparisons
    # 1. create output directory
    logger.info("Generating ANIb command-lines")
    blastdir = os.path.join(os.path.join(args.outdir, pyani_config.ALIGNDIR["ANIb"]))
    try:
        logger.info("Creating output directory %s", blastdir)
        os.makedirs(blastdir, exist_ok=True)
    except IOError:
        logger.error("Could not create output directory %s (exiting)", blastdir)
        logger.error(last_exception())
        raise SystemExit(1)

    # 2. generate (Genome, Genome) tuples for all comparisons
    logger.info("Compiling pairwise comparisons")
    genomes = run.genomes.all()
    logger.info("Collated %d genomes for this run", len(genomes))
    comparisons = list(permutations(tqdm(genomes, disable=args.disable_tqdm), 2))
    logger.info("Total pairwise comparisons to perform: %d", len(comparisons))

    # 3. check for existing comparisons with the same software package, version
    #    and settings; these will be removed from the current list of comparisons
    logger.info("Checking database for existing comparison files...")
    comparisons_to_run = filter_existing_comparisons(
        session, run, comparisons, "blastn", blast_version, args.fragsize, None
    )
    logger.info(
        "...after checking, %d comparisons remain to be run", len(comparisons_to_run)
    )

    # 4. if there aren't any comparisons left to run, update the result matrices and
    #    return
    if len(comparisons_to_run) == 0:
        logger.info(
            "All comparisons are present in the database (skipping comparisons)"
        )
        logger.info("Updating summary matrices with existing results")
        update_comparison_matrices(session, run)
        return

    # 5. if we're in recovery mode, we want to salvage output from an old run
    #    and may not need to run all jobs. If this is so, we prepare a list of
    #    existing output files
    if args.recovery:
        logger.warning("Entering recovery mode")
        logger.info(
            "\tIn this mode, existing comparison output from %s is reused", blastdir
        )
        existingfiles = collect_existing_output(blastdir, "blastn", args)
    else:
        existingfiles = None
    print(existingfiles)

    # 6. construct BLAST nucleotide databases with the appropriate fragment
    #    size for each genome, and add these to the database
    fragdir = os.path.join(args.outdir, "fragments")
    os.makedirs(fragdir, exist_ok=True)
    logger.info(
        "Fragmenting all FASTA files in %s to %dnt in %s",
        args.indir,
        args.fragsize,
        fragdir,
    )
    fragdata = anib.fragment_genomes(
        args.indir, args.fragsize, fragdir, args.makeblastdb_exe
    )

    # Build the BLAST+ nucleotide databases in-place, handing off commands
    # to the appropriate scheduler
    logger.info("Building databases for BLAST+ comparisons")
    joblist = anib.generate_makeblastdb_jobs(fragdata)
    run_jobs(joblist, args, logger)

    # Add the created BLAST+ databases to the pyani database, recording relative
    # path (so that the database can be bundled with results, if a suitable pyani
    # analysis structure is followed). This is transactional, so either all
    # databases are added (with reference to the current run) or none are.
    logger.info("Adding BLAST+ databases to pyani database")
    for fraginfo in fragdata:
        print(fraginfo)
        fragsource = (
            session.query(Genome).filter(Genome.genome_hash == fraginfo.hash).first()
        )
        print(fragsource)
        blastdb = BlastDB(
            genome=fragsource,
            run=run,
            path=fraginfo.dbpath,
            size=fraginfo.fragcount,
            dbcmd=fraginfo.formatcmd,
        )
        session.add(blastdb)
    session.commit()


def run_jobs(joblist, args, logger):
    """Pass ANIb jobs to a scheduler

    joblist           list of Jobs
    args              command-line arguments for the run
    logger            logging output
    """
    if args.scheduler == "multiprocessing":
        logger.info("Running jobs with multiprocessing")
        if not args.workers:
            logger.info("(using maximum number of worker threads)")
        else:
            logger.info("(using %d worker threads, if available)", args.workers)
        cumval = run_mp.run_dependency_graph(
            [_ for _ in joblist], workers=args.workers, logger=logger
        )
        if 0 < cumval:
            logger.error("At least one command failed, please investigate (exiting)")
            raise PyaniException("Multiprocessing run failed in ANIb")
        else:
            logger.info("Multiprocessing run completed without error")
    else:
        logger.info("Running jobs with SGE")
        logger.info("Setting jobarray group size to %d", args.sgegroupsize)
        run_sge.run_dependency_graph(
            [_.job for _ in joblist],
            logger=logger,
            jgprefix=args.jobprefix,
            sgegroupsize=args.sgegroupsize,
            sgeargs=args.sgeargs,
        )
