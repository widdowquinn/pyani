#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
# (c) University of Strathclyde 2019-2020
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# 161 Cathedral Street,
# Glasgow,
# G4 0RE
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2017-2019 The James Hutton Institute
# Copyright (c) 2019-2020 University of Strathclyde
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
"""Provides the anib subcommand for pyani."""

import datetime
import json
import logging
import os
import subprocess

from argparse import Namespace
from itertools import permutations
from pathlib import Path
from typing import List, NamedTuple, Tuple, Dict

from Bio import SeqIO
from tqdm import tqdm

from pyani import (
    PyaniException,
    anib,
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


def subcmd_anib(args: Namespace) -> None:
    """Perform ANIb on all genome files in an input directory.

    :param args:  Namespace, command-line arguments

    Finds ANI by the ANIb method, as described in Goris J, Konstantinidis KT,
    Klappenbach JA, Coenye T, Vandamme P, et al. (2007) DNA-DNA hybridization
    values and their relationship to whole-genome sequence similarities.
    Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0.

    All FASTA format files (selected by suffix) in the input directory are
    fragmented into (by default 1020nt) consecutive sections, and a BLAST+
    database constructed from the whole genome input. The BLAST+ blastn tool
    is then used to query each set of fragments against each BLAST+ database,
    in turn.

    For each query, the BLAST+ .tab output is parsed to obtain alignment length,
    identity and similarity error count. Alignments below a threshold are not
    included in the calculation (this introduces systematic bias with respect to
    ANIm). The results are processed to calculate the ANI percentages, coverage,
    and similarity error.

    The calculated values are stored in the local SQLite3 database.
    """
    # Create logger
    logger = logging.getLogger(__name__)

    # Announce the analysis
    logger.info(termcolor("Running ANIb analysis", "red"))

    # Get BLAST+ version - this will be used in the database entries
    blastn_version = anib.get_version(args.blastn_exe)
    logger.info(termcolor("BLAST+ blastn version: %s", "cyan"), blastn_version)

    # Use provided name, or make new one for this analysis
    start_time = datetime.datetime.now()
    name = args.name or "_".join(["ANIb", start_time.isoformat()])
    logger.info(termcolor("Analysis name: %s", "cyan"), name)

    # Connect to existing database (which may be "clean" or have old analyses)
    logger.debug("Connecting to database %s", args.dbpath)
    try:
        session = get_session(args.dbpath)
    except Exception:
        logger.error(
            "Could not connect to database %s (exiting)", args.dbpath, exc_info=True
        )
        raise SystemExit(1)

    # Add information about this run to the database
    logger.debug("Adding run info to database %s...", args.dbpath)
    try:
        run = add_run(
            session,
            method="ANIb",
            cmdline=args.cmdline,
            date=start_time,
            status="started",
            name=name,
        )
    except PyaniORMException:
        logger.error("Could not add run to the database (exiting)", exc_info=True)
        raise SystemExit(1)
    logger.debug("\t...added run ID: %s to the database", run)

    # Identify input files for comparison, and populate the database
    logger.debug("Adding files for %s to database...", run)
    try:
        genome_ids = add_run_genomes(
            session, run, args.indir, args.classes, args.labels
        )
    except PyaniORMException:
        logger.error(
            "Could not add genomes to database for run %s (exiting)", run, exc_info=True
        )
        raise SystemExit(1)
    logger.debug("\t...added genome IDs: %s", genome_ids)

    # Get list of genomes for this analysis from the database
    logger.info("Compiling genomes for comparison")
    genomes = run.genomes.all()
    logger.debug("\tCollected %s genomes for this run", len(genomes))

    # Create output directories. We create the main parent directory (args.outdir), but
    # also subdirectories for the BLAST databases,
    logger.debug("Creating output directory %s", args.outdir)
    try:
        os.makedirs(args.outdir, exist_ok=True)
    except IOError:
        logger.error(
            "Could not create output directory %s (exiting)", args.outdir, exc_info=True
        )
        raise SystemError(1)
    fragdir = Path(str(args.outdir)) / "fragments"
    blastdbdir = Path(str(args.outdir)) / "blastdbs"
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
        fragfiles.update({genome: fragpath})
        fraglens.update({genome: fragsizes})

        dbcmd, blastdbpath = anib.construct_makeblastdb_cmd(
            Path(genome.path), blastdbdir  # fragpath, blastdbdir
        )  # args.outdir)

        subprocess.run(
            dbcmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

        add_blastdb(
            session, genome, run, fragpath, blastdbpath, json.dumps(fragsizes), dbcmd
        )

    # Generate all pair permutations of genome IDs as a list of (Genome, Genome) tuples
    logger.info(
        "Compiling pairwise comparisons (this can take time for large datasets)..."
    )
    comparisons = list(permutations(tqdm(genomes, disable=args.disable_tqdm), 2))
    logger.info("\t...total pairwise comparisons to be performed: %s", len(comparisons))

    # Check for existing comparisons; if one has already been done (for the same
    # software package, version, and setting) we add the comparison to this run,
    # but remove it from the list of comparisons to be performed
    logger.info("Checking database for existing comparison data...")
    comparisons_to_run = filter_existing_comparisons(
        session, run, comparisons, "blastn", blastn_version, args.fragsize, False
    )
    logger.info(
        "\t...after check, still need to run %s comparisons", len(comparisons_to_run)
    )

    # If there are no comparisons to run, update the Run matrices and exit
    # from this function
    if not comparisons_to_run:
        logger.info(
            termcolor(
                "All comparison results present in database (skipping comparisons)",
                "magenta",
            )
        )
        logger.info("Updating summary matrices with existing results")
        update_comparison_matrices(session, run)
        return

    # If we are in recovery mode, we are salvaging output from a previous
    # run, and do not necessarily need to rerun all the jobs. In this case,
    # we prepare a list of output files we want to recover from the results
    # in the output directory.
    # ¶ Should this use output files, or pull from the database?
    if args.recovery:
        logger.warning("Entering recovery mode...")
        logger.debug(
            "\tIn this mode, existing comparison output from %s is reused", args.outdir
        )
        existing_files = collect_existing_output(args.outdir, "blastn", args)
        if existing_files:
            logger.debug(
                "\tIdentified %s existing output files for reuse, %s (etc)",
                len(existing_files),
                existing_files[0],
            )
        else:
            logger.debug("\tIdentified no existing output files")
    else:
        existing_files = list()
        logger.debug("\tAssuming no pre-existing output files")

    ## Split the input genome files into contiguous fragments of the specified size,
    ## as described in Goris et al. We create a new directory to hold sequence
    ## fragments, away from the main genomes
    # logger.info("Splitting input genome files into %snt fragments...", args.fragsize)
    # fragdir = Path(args.outdir) / "fragments"
    # os.makedirs(fragdir, exist_ok=True)
    # fragfiles, fraglens = anib.fragment_fasta_files(
    #     [Path(str(_.path)) for _ in genomes],
    #     Path(args.outdir) / "fragments",
    #     args.fragsize,
    # )
    # logger.debug("...wrote %s fragment files to %s", len(fragfiles), fragdir)

    # Create list of BLASTN jobs for each comparison still to be performed
    logger.info("Creating blastn jobs for ANIb...")
    # This method considered, but doesn't get fraglens
    # fragfiles =  [file for file in fragdir.iterdir()]
    joblist = generate_joblist(
        comparisons_to_run, existing_files, fragfiles.values(), fraglens.values(), args
    )
    logger.debug("...created %s blastn jobs", len(joblist))

    # Pass jobs to appropriate scheduler
    logger.debug("Passing %s jobs to %s...", len(joblist), args.scheduler)
    run_anib_jobs(joblist, args)
    logger.info("...jobs complete.")

    # Process output and add results to database
    # This requires us to drop out of threading/multiprocessing: Python's SQLite3
    # interface doesn't allow sharing connections and cursors
    logger.info("Adding comparison results to database...")
    update_comparison_results(joblist, run, session, blastn_version, fraglens, args)
    update_comparison_matrices(session, run)
    logger.info("...database updated.")


def generate_joblist(
    comparisons: List,
    existing_files: List,
    fragfiles: List,
    fragsizes: List,
    args: Namespace,
) -> List[ComparisonJob]:
    """Return list of ComparisonJobs.

    :param comparisons:  list of (Genome, Genome) tuples for which comparisons are needed
    :param existing_files:  list of pre-existing BLASTN+ outputs
    :param fragfiles:  list of files containing genome fragments
    :param fragsizes:  list of fragment lengths
    :param args:  Namespace, command-line arguments
    """
    logger = logging.getLogger(__name__)

    existing_files = set(existing_files)  # Path objects hashable

    joblist = []  # will hold ComparisonJob structs
    for idx, (query, subject) in enumerate(
        tqdm(comparisons, disable=args.disable_tqdm)
    ):
        qprefix, qsuffix = (
            args.outdir / "fragments" / Path(query.path).stem,
            Path(query.path).suffix,
        )
        qfrags = Path(f"{qprefix}-fragments{qsuffix}")

        blastcmd = anib.generate_blastn_commands(
            qfrags, Path(subject.path), args.outdir, args.blastn_exe
        )
        logger.debug("Commands to run:\n\t%s\n", blastcmd)
        outprefix = blastcmd.split()[2][:-10]  # prefix for blastn output
        outfname = Path(outprefix + ".blast_tab")
        logger.debug("Expected output file for db: %s", outfname)

        # If we are in recovery mode, we are salvaging output from a previous
        # run, and do not necessarily need to rerun all the jobs. In this case,
        # we prepare a list of output files we want to recover from the results
        # in the output directory.

        # The comparisons collection always gets updated, so that results are
        # added to the database whether they come from recovery mode or are run
        # in this call of the script.
        if args.recovery and outfname in existing_files:
            logger.debug("Recovering output from %s, not submitting job", outfname)
            # Need to track the expected output, but set the job itself to None.
            joblist.append(
                ComparisonJob(query, subject, blastcmd, outfname, args.fragsize, None)
            )
        else:
            logger.debug("Building job")
            # Build job
            blastjob = pyani_jobs.Job("%s_%06d-blast" % (args.jobprefix, idx), blastcmd)
            joblist.append(
                ComparisonJob(
                    query, subject, blastcmd, outfname, args.fragsize, blastjob
                )
            )
    return joblist
    # raise NotImplementedError


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


def run_anib_jobs(joblist: List[ComparisonJob], args: Namespace) -> None:
    """Pass ANIb blastn jobs to the scheduler.

    :param joblist:           list of ComparisonJob namedtuples
    :param args:              command-line arguments for the run
    """
    logger = logging.getLogger(__name__)
    logger.debug("Scheduler: %s", args.scheduler)

    # Entries with None seen in recovery mode:
    jobs = [_.job for _ in joblist if _.job]

    if args.scheduler == "multiprocessing":
        logger.info("Running jobs with multiprocessing")
        if not args.workers:
            logger.debug("(using maximum number of worker threads)")
        else:
            logger.debug("(using %d worker threads, if available)", args.workers)
        cumval = run_mp.run_dependency_graph(jobs, workers=args.workers)
        if cumval > 0:
            logger.error(
                "At least one blastn comparison failed. Please investigate (exiting)"
            )
            raise PyaniException("Multiprocessing run failed in ANIb")
        logger.info("Multiprocessing run completed without error")
    elif args.scheduler.lower() == "sge":
        logger.info("Running jobs with SGE")
        logger.debug("Setting jobarray group size to %d", args.sgegroupsize)
        logger.debug("Joblist contains %d jobs", len(joblist))
        run_sge.run_dependency_graph(
            jobs,
            jgprefix=args.jobprefix,
            sgegroupsize=args.sgegroupsize,
            sgeargs=args.sgeargs,
        )
    else:
        logger.error(termcolor("Scheduler %s not recognised", "red"), args.scheduler)
        raise SystemError(1)


#
#
def update_comparison_results(
    joblist: List[ComparisonJob],
    run,
    session,
    blastn_version: str,
    fraglens: Dict,
    args: Namespace,
) -> None:
    """Update the Comparision table with the completed result set.

    :param joblist:         list of ComparisonJob namedtuples
    :param run:             Run ORM object for the current ANIb run
    :param session:         active pyanidb session via ORM
    :param blastn_version:  version of blastn used for the comparison
    :param fraglens:     dictionary of fragment lengths for each genome
    :param args:            command-line arguments for this run

    The Comparison table stores individual comparison results, one per row.
    """
    logger = logging.getLogger(__name__)

    # Add individual results to Comparison table
    for job in tqdm(joblist, disable=args.disable_tqdm):
        logger.debug("\t%s vs %s", job.query.description, job.subject.description)
        aln_length, sim_errs, ani_pid = anib.parse_blast_tab(job.outfile, fraglens)

        qcov = aln_length / job.query.length
        scov = aln_length / job.subject.length
        run.comparisons.append(
            Comparison(
                query=job.query,
                subject=job.subject,
                aln_length=int(aln_length),
                sim_errs=int(sim_errs),
                identity=ani_pid,
                cov_query=qcov,
                cov_subject=scov,
                program="blastn",
                version=blastn_version,
                fragsize=job.fragsize,
                maxmatch=False,
            )
        )
    #
    # Populate db
    logger.debug("Committing results to database")
    session.commit()
