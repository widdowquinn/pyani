#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
#
# Author: Leighton Pritchard
# Contact: leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
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
"""Provides the anim subcommand for pyani."""

import datetime
import os

from collections import namedtuple
from itertools import combinations
from pathlib import Path

from tqdm import tqdm

from pyani import (
    PyaniException,
    anim,
    pyani_config,
    pyani_jobs,
    run_sge,
    run_multiprocessing as run_mp,
)
from pyani.pyani_files import collect_existing_output
from pyani.pyani_orm import (
    Comparison,
    PyaniORMException,
    add_run,
    add_run_genomes,
    filter_existing_comparisons,
    get_session,
    update_comparison_matrices,
)


# Convenience struct describing a pairwise comparison job for the SQLAlchemy
# implementation
ComparisonJob = namedtuple(
    "ComparisonJob", "query subject filtercmd nucmercmd outfile job"
)

# Convenience struct describing an analysis run
RunData = namedtuple("RunData", "method name date cmdline")

# Convenience struct for a single nucmer comparison result
ComparisonResult = namedtuple(
    "ComparisonResult", "qid sid aln_length sim_errs pid qlen slen qcov scov"
)

# Convenience struct for comparison program data/info
ProgData = namedtuple("ProgData", "program version")

# Convenience struct for comparison parameters
# Use default of zero for fragsize or else db queries will not work
# SQLite/Python nulls do not match up well
ProgParams = namedtuple("ProgParams", "fragsize maxmatch")


def subcmd_anim(args, logger):
    """Perform ANIm on all genome files in an input directory.

    :param args:  Namespace, command-line arguments
    :param logger:  logging object

    Finds ANI by the ANIm method, as described in Richter et al (2009)
    Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

    All FASTA format files (selected by suffix) in the input directory
    are compared against each other, pairwise, using NUCmer (whose path must
    be provided).

    For each pairwise comparison, the NUCmer .delta file output is parsed to
    obtain an alignment length and similarity error count for every unique
    region alignment between the two organisms, as represented by
    sequences in the FASTA files. These are processed to calculated aligned
    sequence lengths, average nucleotide identity (ANI) percentages, coverage
    (aligned percentage of whole genome - forward direction), and similarity
    error count for each pairwise comparison.

    The calculated values are deposited in the SQLite3 database being used for
    the analysis.

    For each pairwise comparison the NUCmer output is stored in the output
    directory for long enough to extract summary information, but for each run
    the output is gzip compressed. Once all runs are complete, the outputs
    for each comparison are concatenated into a single gzip archive.
    """
    # Announce the analysis
    logger.info("Running ANIm analysis")

    # Get current nucmer version
    nucmer_version = anim.get_version(args.nucmer_exe)
    logger.info(f"MUMMer nucmer version: {nucmer_version}")

    # Use the provided name or make one for the analysis
    start_time = datetime.datetime.now()
    name = args.name or "_".join(["ANIm", start_time.isoformat()])
    logger.info(f"Analysis name: {name}")

    # Get connection to existing database. This may or may not have data
    logger.info(f"Connecting to database {args.dbpath}")
    try:
        session = get_session(args.dbpath)
    except Exception:
        logger.error(
            f"Could not connect to database {args.dbpath} (exiting)", exc_info=True
        )
        raise SystemExit(1)

    # Add information about this run to the database
    logger.info(f"Adding run info to database {args.dbpath}...")
    try:
        run = add_run(
            session,
            method="ANIm",
            cmdline=args.cmdline,
            date=start_time,
            status="started",
            name=name,
        )
    except PyaniORMException:
        logger.error(
            f"Could not add run {run} to the database (exiting)", exc_info=True
        )
        raise SystemExit(1)
    logger.info(f"...added run ID: {run} to the database")

    # Identify input files for comparison, and populate the database
    logger.info(f"Adding genomes for run {run} to database...")
    try:
        genome_ids = add_run_genomes(
            session, run, args.indir, args.classes, args.labels
        )
    except PyaniORMException:
        logger.error(f"Could not add genomes to database for run {run} (exiting)")
        raise SystemExit(1)
    logger.info(f"\t...added genome IDs: {genome_ids}")

    # Generate commandlines for NUCmer analysis and output compression
    logger.info("Generating ANIm command-lines")
    deltadir = os.path.join(os.path.join(args.outdir, pyani_config.ALIGNDIR["ANIm"]))
    logger.info("NUCmer output will be written temporarily to %s", deltadir)

    # Create output directories
    logger.info(f"Creating output directory {deltadir}")
    try:
        os.makedirs(deltadir, exist_ok=True)
    except IOError:
        logger.error(
            f"Could not create output directory {deltadir} (exiting)", exc_info=True
        )
        raise SystemError(1)

    # Get list of genome IDs for this analysis from the database
    logger.info("Compiling genomes for comparison")
    genomes = run.genomes.all()
    logger.info(f"Collected {len(genomes)} genomes for this run")

    # Generate all pair combinations of genome IDs as a list of (Genome, Genome) tuples
    logger.info(
        "Compiling pairwise comparisons (this can take time for large datasets)..."
    )
    comparisons = list(combinations(tqdm(genomes, disable=args.disable_tqdm), 2))
    logger.info(f"\t...total parwise comparisons to be performed: {len(comparisons)}")

    # Check for existing comparisons; if one has been done (for the same
    # software package, version, and setting) we add the comparison to this run,
    # but remove it from the list of comparisons to be performed
    logger.info("Checking database for existing comparison data...")
    comparisons_to_run = filter_existing_comparisons(
        session, run, comparisons, "nucmer", nucmer_version, None, args.maxmatch
    )
    logger.info(
        f"\t...after check, still need to run {len(comparisons_to_run)} comparisons"
    )

    # If there are no comparisons to run, update the Run matrices and exit
    # from this function
    if not comparisons_to_run:
        logger.info("All comparison results present in database (skipping comparisons)")
        logger.info("Updating summary matrices with existing results")
        update_comparison_matrices(session, run)
        return

    # If we are in recovery mode, we are salvaging output from a previous
    # run, and do not necessarily need to rerun all the jobs. In this case,
    # we prepare a list of output files we want to recover from the results
    # in the output directory.
    if args.recovery:
        logger.warning("Entering recovery mode")
        logger.info(
            f"\tIn this mode, existing comparison output from {deltadir} is reused"
        )
        existingfiles = collect_existing_output(deltadir, "nucmer", args)
        logger.info(
            f"\tIdentified {len(existingfiles)} existing output files for reuse"
        )
    else:
        existingfiles = None
        logger.info(f"\tIdentified no existing output files")

    # Create list of NUCmer jobs for each comparison still to be performed
    logger.info("Creating NUCmer jobs for ANIm")
    joblist = generate_joblist(comparisons_to_run, existingfiles, args, logger)
    logger.info(f"Generated {len(joblist)} jobs, {len(comparisons_to_run)} comparisons")

    # Pass jobs to appropriate scheduler
    logger.info(f"Passing {len(joblist)} jobs to {args.scheduler}...")
    run_anim_jobs(joblist, args, logger)
    logger.info("...jobs complete")

    # Process output and add results to database
    # This requires us to drop out of threading/multiprocessing: Python's SQLite3
    # interface doesn't allow sharing connections and cursors
    logger.info("Adding comparison results to database...")
    update_comparison_results(joblist, run, session, nucmer_version, args, logger)
    update_comparison_matrices(session, run)
    logger.info("...database updated.")


def generate_joblist(comparisons, existingfiles, args, logger):
    """Return list of ComparisonJobs.

    :param comparisons:  list of (Genome, Genome) tuples
    :param existingfiles:  list of pre-existing nucmer output files
    :param args:  Namespace of command-line arguments for the run
    :param logger:  logging object
    """
    joblist = []  # will hold ComparisonJob structs
    for idx, (query, subject) in enumerate(
        tqdm(comparisons, disable=args.disable_tqdm)
    ):
        ncmd, dcmd = anim.construct_nucmer_cmdline(
            query.path,
            subject.path,
            args.outdir,
            args.nucmer_exe,
            args.filter_exe,
            args.maxmatch,
        )
        logger.debug("Commands to run:\n\t%s\n\t%s", ncmd, dcmd)
        outprefix = ncmd.split()[3]  # prefix for NUCmer output
        if args.nofilter:
            outfname = Path(outprefix + ".delta")
        else:
            outfname = Path(outprefix + ".filter")
        logger.debug("Expected output file for db: %s", outfname)

        # If we're in recovery mode, we don't want to repeat a computational
        # comparison that already exists, so we check whether the ultimate
        # output is in the set of existing files and, if not, we add the jobs
        # TODO: something faster than a list search (dict or set?)
        # The comparisons collections always gets updated, so that results are
        # added to the database whether they come from recovery mode or are run
        # in this call of the script.
        if args.recovery and os.path.split(outfname)[-1] in existingfiles:
            logger.debug("Recovering output from %s, not building job", outfname)
        else:
            logger.debug("Building job")
            # Build jobs
            njob = pyani_jobs.Job("%s_%06d-n" % (args.jobprefix, idx), ncmd)
            fjob = pyani_jobs.Job("%s_%06d-f" % (args.jobprefix, idx), dcmd)
            fjob.add_dependency(njob)
            joblist.append(ComparisonJob(query, subject, dcmd, ncmd, outfname, fjob))
    return joblist


def run_anim_jobs(joblist, args, logger):
    """Pass ANIm nucmer jobs to the scheduler.

    :param joblist:           list of ComparisonJob namedtuples
    :param args:              command-line arguments for the run
    :param logger:            logging output
    """
    if args.scheduler == "multiprocessing":
        logger.info("Running jobs with multiprocessing")
        if not args.workers:
            logger.info("(using maximum number of worker threads)")
        else:
            logger.info("(using %d worker threads, if available)", args.workers)
        cumval = run_mp.run_dependency_graph(
            [_.job for _ in joblist], workers=args.workers, logger=logger
        )
        if cumval > 0:
            logger.error(
                "At least one NUCmer comparison failed. "
                + "Please investigate (exiting)"
            )
            raise PyaniException("Multiprocessing run failed in ANIm")
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


def update_comparison_results(joblist, run, session, nucmer_version, args, logger):
    """Update the Comparison table with the completed result set.

    :param joblist:         list of ComparisonJob namedtuples
    :param run:             Run ORM object for the current ANIm run
    :param session:         active pyanidb session via ORM
    :param nucmer_version:  version of nucmer used for the comparison
    :param args:            command-line arguments for this run
    :param logger:          logging output

    The Comparison table stores individual comparison results, one per row.
    """
    # Add individual results to Comparison table
    for job in tqdm(joblist, disable=args.disable_tqdm):
        logger.debug("\t%s vs %s", job.query.description, job.subject.description)
        aln_length, sim_errs = anim.parse_delta(job.outfile)
        qcov = aln_length / job.query.length
        scov = aln_length / job.subject.length
        try:
            pid = 1 - sim_errs / aln_length
        except ZeroDivisionError:  # aln_length was zero (no alignment)
            pid = 0
        run.comparisons.append(
            Comparison(
                query=job.query,
                subject=job.subject,
                aln_length=aln_length,
                sim_errs=sim_errs,
                identity=pid,
                cov_query=qcov,
                cov_subject=scov,
                program="nucmer",
                version=nucmer_version,
                fragsize=None,
                maxmatch=args.maxmatch,
            )
        )

    # Populate db
    logger.info("Committing results to database")
    session.commit()
