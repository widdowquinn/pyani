#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""subcmd_anim.py

Provides the anim subcommand for pyani

(c) The James Hutton Institute 2017-2018

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

Copyright (c) 2017-2018 The James Hutton Institute

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
import sqlite3

from collections import namedtuple
from itertools import combinations

from tqdm import tqdm

from pyani import (
    anim,
    pyani_config,
    pyani_orm,
    pyani_files,
    pyani_jobs,
    pyani_tools,
    run_sge,
    run_multiprocessing as run_mp,
)
from pyani.pyani_orm import Run, Genome, Label, Class
from pyani.pyani_tools import last_exception


# Named tuple describing a pairwise comparison
Comparison = namedtuple("Comparison", "query_id subject_id cmdline outfile")


def subcmd_anim(args, logger):
    """Perform ANIm on all genome files in an input directory.

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

    # Use the provided name or make one for the analysis
    start_time = datetime.datetime.now()
    if args.name is None:
        name = "_".join(["ANIm", start_time.isoformat()])
    else:
        name = args.name

    # Add info for this analysis to the database
    # First, get a connection to an existing database (which may or may not have data)
    logger.info("Connecting to database %s", args.dbpath)
    try:
        session = pyani_orm.get_session(args.dbpath)
    except Exception:
        logger.error("Could not connect to database %s (exiting)", args.dbpath)
        raise SystemExit(1)

    # Add this run to the database
    logger.info("Adding run information to database %s", args.dbpath)
    run = Run(
        method="ANIm",
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

    # Identify input files for comparison, and populate the database
    logger.info("Adding this run's genomes to database")
    add_run_genomes(session, run, args, logger)
    logger.info("Input genome/hash files added to database")

    # Generate commandlines for NUCmer analysis and output compression
    logger.info("Generating ANIm command-lines")
    deltadir = os.path.join(os.path.join(args.outdir, pyani_config.ALIGNDIR["ANIm"]))
    logger.info("NUCmer output will be written temporarily to %s", deltadir)

    # Create output directories
    logger.info("Creating output directory %s", deltadir)
    try:
        os.makedirs(deltadir, exist_ok=True)
    except IOError:
        logger.error("Could not create output directory (exiting)")
        logger.error(last_exception())
        raise SystemError(1)

    # Get list of genome IDs for this analysis from the database
    logger.info("Compiling genomes for comparison")
    genomes = run.genomes.all()
    logger.info("Collated {} genomes for this run".format(len(genomes)))

    # Generate all pair combinations of genome IDs
    logger.info(
        "Compiling pairwise comparisons (this can take time for large datasets)"
    )
    comparisons = list(combinations(tqdm(genomes, disable=args.disable_tqdm), 2))
    logger.debug("Total pairwise comparisons {}".format(len(comparisons)))

    # Check for existing comparisons; if one has been done (for the same
    # software package, version, and setting) we remove it from the list
    # of comparisons to be performed, but we add a new entry to the
    # runs_comparisons table.
    # TODO: turn this into a generator or some such?
    # nucmer_version = anim.get_version(args.nucmer_exe)

    import sys

    sys.exit(0)

    # # Existing entries for the comparison:run link table
    # new_link_ids = [
    #     (genome1, genome2)
    #     for (genome1, genome2) in comparisons
    #     if pyani_orm.get_comparison(
    #         args.dbpath, qid, sid, "nucmer", nucmer_version, maxmatch=args.maxmatch
    #     )
    #     is not None
    # ]
    # logger.info(
    #     "Existing comparisons to be associated with new run:\n\t%s", new_link_ids
    # )
    # if len(new_link_ids) > 0:
    #     for (qid, sid) in tqdm(new_link_ids, disable=args.disable_tqdm):
    #         pyani_orm.add_comparison_link(
    #             args.dbpath,
    #             run_id,
    #             qid,
    #             sid,
    #             "nucmer",
    #             nucmer_version,
    #             maxmatch=args.maxmatch,
    #         )

    # # If we are in recovery mode, we are salvaging output from a previous
    # # run, and do not necessarily need to rerun all the jobs. In this case,
    # # we prepare a list of output files we want to recover from the results
    # # in the output directory.
    # if args.recovery:
    #     logger.warning("Entering recovery mode")
    #     logger.info(
    #         "\tIn this mode, existing comparison output from %s is reused", deltadir
    #     )
    #     # Obtain collection of expected output files already present in directory
    #     if args.nofilter:
    #         suffix = ".delta"
    #     else:
    #         suffix = ".filter"
    #     existingfiles = [
    #         fname
    #         for fname in os.listdir(deltadir)
    #         if os.path.splitext(fname)[-1] == suffix
    #     ]
    #     logger.info("Identified %d existing output files", len(existingfiles))
    # else:
    #     existingfiles = None

    # # New comparisons to be run for this analysis
    # # TODO: Can we parallelise this as a function called with multiprocessing?
    # logger.info("Excluding comparisons present in database")
    # comparison_ids = [
    #     (qid, sid)
    #     for (qid, sid) in comparison_ids
    #     if pyani_orm.get_comparison(
    #         args.dbpath, qid, sid, "nucmer", nucmer_version, maxmatch=args.maxmatch
    #     )
    #     is None
    # ]
    # logger.debug("Comparisons still to be performed:\n\t%s", comparison_ids)
    # logger.info("Total comparisons to be conducted: %d", len(comparison_ids))

    # if not len(comparison_ids):
    #     logger.info(
    #         "All comparison results already present in database "
    #         + "(skipping comparisons)"
    #     )
    # else:
    #     # Create list of NUCmer jobs for each comparison still to be
    #     # performed
    #     logger.info("Creating NUCmer jobs for ANIm")
    #     joblist, comparisons = generate_joblist(
    #         comparison_ids, existingfiles, args, logger
    #     )
    #     logger.info("Generated %s jobs, %d comparisons", len(joblist), len(comparisons))

    #     # Pass commands to the appropriate scheduler
    #     logger.info("Passing %d jobs to scheduler", len(joblist))
    #     if args.scheduler == "multiprocessing":
    #         logger.info("Running jobs with multiprocessing")
    #         if not args.workers:
    #             logger.info("(using maximum number of worker threads)")
    #         else:
    #             logger.info("(using %d worker threads, if available)", args.workers)
    #         cumval = run_mp.run_dependency_graph(
    #             joblist, workers=args.workers, logger=logger
    #         )
    #         if 0 < cumval:
    #             logger.error(
    #                 "At least one NUCmer comparison failed. "
    #                 + "Please investigate (exiting)"
    #             )
    #             raise pyani_tools.PyaniException(
    #                 "Multiprocessing run " + "failed in ANIm"
    #             )
    #         else:
    #             logger.info("Multiprocessing run completed without error")
    #     else:
    #         logger.info("Running jobs with SGE")
    #         logger.info("Setting jobarray group size to %d", args.sgegroupsize)
    #         run_sge.run_dependency_graph(
    #             joblist,
    #             logger=logger,
    #             jgprefix=args.jobprefix,
    #             sgegroupsize=args.sgegroupsize,
    #             sgeargs=args.sgeargs,
    #         )

    #     # Process output and add results to database
    #     # We have to drop out of threading/multiprocessing to do this: Python's
    #     # SQLite3 interface doesn't allow sharing connections and cursors
    #     logger.info("Adding comparison results to database")
    #     for comp in tqdm(comparisons, disable=args.disable_tqdm):
    #         aln_length, sim_errs = anim.parse_delta(comp.outfile)
    #         qlen = pyani_orm.get_genome_length(args.dbpath, comp.query_id)
    #         slen = pyani_orm.get_genome_length(args.dbpath, comp.subject_id)
    #         qcov = aln_length / qlen
    #         scov = aln_length / slen
    #         pid = 1 - sim_errs / aln_length
    #         comp_id = pyani_orm.add_comparison(
    #             args.dbpath,
    #             comp.query_id,
    #             comp.subject_id,
    #             aln_length,
    #             sim_errs,
    #             pid,
    #             qcov,
    #             scov,
    #             "nucmer",
    #             nucmer_version,
    #             maxmatch=args.maxmatch,
    #         )
    #         link_id = pyani_orm.add_comparison_link(
    #             args.dbpath,
    #             run_id,
    #             comp.query_id,
    #             comp.subject_id,
    #             "nucmer",
    #             nucmer_version,
    #             maxmatch=args.maxmatch,
    #         )
    #         logger.debug(
    #             "Added ID %s vs %s, as comparison %s (link: %s)",
    #             comp.query_id,
    #             comp.subject_id,
    #             comp_id,
    #             link_id,
    #         )


def load_classes_labels(path):
    """Returns a dictionary of genome classes or labels keyed by hash

    The expected format of the classes and labels files is:

    <HASH>\t<FILESTEM>\t<CLASS>|<LABEL>,

    where <HASH> is the MD5 hash of the genome data (this is not checked);
    <FILESTEM> is the path to the genome file (this is intended to be a
    record for humans to audit, it's not needed for the database interaction;
    and <CLASS>|<LABEL> is the class or label associated with that genome.
    """
    datadict = {}
    with open(path, "r") as ifh:
        for line in ifh.readlines():
            genomehash, _, data = line.strip().split("\t")
            datadict[genomehash] = data
    return datadict


def add_run_genomes(session, run, args, logger):
    """Add genomes to a run in the pyani database

    session    live SQLAlchemy session to pyani database
    run        a pyani_orm.Run object describing a pyani run
    args       parsed command-line arguments
    logger     logging object
    """
    infiles = pyani_files.get_fasta_and_hash_paths(args.indir)
    # Get labels and classes, keyed by hash
    if args.classes:
        classes = load_classes_labels(args.classes)
    if args.labels:
        labels = load_classes_labels(args.labels)
    # Get hash string and sequence description for each FASTA/hash pair,
    # and add info to the current database
    for fastafile, hashfile in infiles:
        # Get genome data
        logger.info("Processing genome hash: %s", hashfile)
        inhash, _ = pyani_files.read_hash_string(hashfile)
        logger.info("Processing genome sequence: %s", fastafile)
        indesc = pyani_files.read_fasta_description(fastafile)
        abspath = os.path.abspath(fastafile)
        genome_len = pyani_tools.get_genome_length(abspath)
        outstr = [
            "FASTA file:\t%s" % abspath,
            "description:\t%s" % indesc,
            "hash file:\t%s" % hashfile,
            "MD5 hash:\t%s" % inhash,
            "Total length:\t%d" % genome_len,
        ]
        logger.info("\t" + "\n\t".join(outstr))

        # If it's not already there, add the passed genome to the database and
        # associate with the passed run
        logger.info("Checking if genome is in database...")
        genome = session.query(Genome).filter(Genome.genome_hash == inhash).first()
        if not isinstance(genome, Genome):
            # No existing genome with this hash
            logger.info("Adding genome to database...")
            try:
                genome = Genome(
                    genome_hash=inhash,
                    path=abspath,
                    length=genome_len,
                    description=indesc,
                )
                session.add(genome)
                session.commit()
                logger.info("Added genome %s to database", genome)
            except Exception:
                logger.error("Could not add genome to database (exiting)")
                logger.error(last_exception())
                raise SystemExit(1)
        try:
            logger.info("Connecting run with genome")
            genome.runs.append(run)
        except Exception:
            logger.error("Could not link genome to current run (exiting)")
            logger.error(last_exception())
            raise SystemExit(1)

        # If there is an associated class or label for this genome, add it
        if inhash in classes:
            try:
                gclass = Class(genome=genome, run=run, genome_class=classes[inhash])
                session.add(gclass)
                session.commit()
            except Exception:
                logger.error("Could not add genome class to database (exiting)")
                logger.error(last_exception())
                raise SystemExit(1)
        if inhash in labels:
            try:
                glabel = Label(genome=genome, run=run, genome_label=labels[inhash])
                session.add(glabel)
                session.commit()
            except Exception:
                logger.error("Could not add genome label to database (exiting)")
                logger.error(last_exception())
                raise SystemExit(1)


def generate_joblist(comparison_ids, existingfiles, args, logger):
    """Returns tuple of ANIm jobs, and comparisons."""
    joblist, comparisons = [], []
    for idx, (qid, sid) in enumerate(tqdm(comparison_ids, disable=args.disable_tqdm)):
        qpath = pyani_orm.get_genome_path(args.dbpath, qid)
        spath = pyani_orm.get_genome_path(args.dbpath, sid)
        ncmd, dcmd = anim.construct_nucmer_cmdline(
            qpath, spath, args.outdir, args.nucmer_exe, args.filter_exe, args.maxmatch
        )
        logger.debug("Commands to run:\n\t%s\n\t%s", ncmd, dcmd)
        outprefix = ncmd.split()[3]  # prefix for NUCmer output
        if args.nofilter:
            outfname = outprefix + ".delta"
        else:
            outfname = outprefix + ".filter"
        logger.debug("Expected output file for db: %s", outfname)

        # If we're in recovery mode, we don't want to repeat a computational
        # comparison that already exists, so we check whether the ultimate
        # output is in the set of existing files and, if not, we add the jobs
        # TODO: something faster than a list search (dict or set?)
        # The comparisons collections always gets updated, so that results are
        # added to the database whether they come from recovery mode or are run
        # in this call of the script.
        comparisons.append(Comparison(qid, sid, dcmd, outfname))
        if args.recovery and os.path.split(outfname)[-1] in existingfiles:
            logger.debug("Recovering output from %s, not building job", outfname)
        else:
            logger.debug("Building job")
            # Build jobs
            njob = pyani_jobs.Job("%s_%06d-n" % (args.jobprefix, idx), ncmd)
            fjob = pyani_jobs.Job("%s_%06d-f" % (args.jobprefix, idx), dcmd)
            fjob.add_dependency(njob)
            joblist.append(fjob)
    return (joblist, comparisons)
