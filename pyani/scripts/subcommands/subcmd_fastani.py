# -*- coding: utf-8 -*-
# (c) The University of Strathclyde 2021–Present
# Author: Bailey Harrington
#
# Contact: bailey.harrington@strath.ac.uk
#
# Bailey Harrington,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G4 0RE,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2021–Present University of Strathclyde
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

import datetime
import json  # is this needed? [LP: for this, probably not]
import logging
import os

from argparse import (
    Namespace,
)  # unsure what this does [LP: makes it easier to mock CLI stuff; also mypy hints]
from itertools import (
    permutations,
    combinations,
)  # may be specific to anib [LP: yes; ANIb is asymmetrical, but fastANI will take care of what permutations does if you pass a file list; I'd expect you'd only need this when passing jobs out via SLURM/SGE]
from pathlib import Path, PosixPath
from typing import List, Tuple, NamedTuple
from Bio import SeqIO  # unsure what this does [LP: bioinformatics file format IO]
from tqdm import tqdm  # unsure what this does [LP: progress bars]

from pyani import (
    fastani,
    pyani_jobs,
    run_sge,
    pyani_config,
    run_multiprocessing as run_mp,
)
from pyani.pyani_files import collect_existing_output
from pyani.pyani_orm import (
    Comparison,
    PyaniException,
    PyaniORMException,
    add_run,
    add_run_genomes,
    filter_existing_comparisons,
    get_session,
    update_comparison_matrices,
)
from pyani.pyani_tools import termcolor


# Convenience struct describing a pairwise comparison job for the SQLAlchemy
# implementation
class ComparisonJob(NamedTuple):

    """Pairwise comparison job for the SQLAlchemy implementation"""

    query: str
    ref: str
    fastcmd: str
    outfile: Path
    fragLen: int
    kmerSize: int
    minFraction: float
    job: pyani_jobs.Job


# Convenience struct describing an analysis run
class RunData(NamedTuple):

    """Convenience struct describing an analysis run."""

    method: str
    name: str
    date: datetime.datetime
    cmdline: str


class ComparisonResult(NamedTuple):

    """Convenience struct for a single fastani comparison result."""

    qid: float
    rid: float
    aln_length: int
    sim_errs: int
    pid: float
    qlen: int
    rlen: int
    qcov: float
    rcov: float


def subcmd_fastani(args: Namespace) -> None:
    """Perform fastANI on all genome files in an input directory.

    :param args:  Namespace, command-line arguments

    Finds ANI by the fastANI method, as described in Jain et al (2018)
    Nature Communications 9, 5114. doi:10.1038/s41467-018-07641-9.

    All FASTA format files (selected by suffix) in the input directory
    are compared against each other, pairwise, using fastANI (whose path
    must be provided).

    For each pairwise comparison, the fastANI .fastani file output is parsed
    to obtain an alignment length and similarity error countfor the two
    organisms, as represented by sequences in the FASTA files. These are
    processed to calculate aligned sequence lengths, average nucleotide
    identity (ANI) percentages, coverage (aligned percentage of whole
    genome), and similarity error count for each pairwise comparison.

    The calculated values are deposited in the SQLite3 database being used
    for the analysis.

    For each pairwise comparison, the fastANI output is stored in the output
    directory for long enough to extract summary information, but for each
    run the output is gzip compressed. Once all runs are complete, the
    outputs for each comparison are concatenated into a single gzip archive.

    """
    logger = logging.getLogger(__name__)

    # announce that we're starting
    logger.info(termcolor("Running fastANI analysis", "red"))

    # Get current fastani version
    logger.info(termcolor("fastANI executable: %s", args.fastani_exe))
    fastani_version = fastani.get_version(args.fastani_exe)
    logger.info(termcolor("FastANI version: %s", "cyan"), fastani_version)

    # Use provided name, or make new one for this analysis
    start_time = datetime.datetime.now()
    name = args.name or "_".join(["FastANI", start_time.isoformat()])
    logger.info(termcolor("Analysis name: %s", "cyan"), name)

    # Connect to existing database (which may be "clean", or have old analyses)
    logger.debug("Connecting to database %s", args.dbpath)
    try:
        session = get_session(args.dbpath)
        session.echo = True
    except Exception:  # is there a less generic option?
        logger.error(
            "Could not connect to database %s; exiting", args.dbpath, exc_info=True
        )
        raise SystemExit(1)

    # Add information about this run to the database
    logger.debug("Adding run info to database %s...", args.dbpath)
    try:
        run, run_id = add_run(
            session,
            method="FastANI",
            cmdline=args.cmdline,
            date=start_time,
            status="started",
            name=name,
        )
    except PyaniORMException:
        logger.error("Could not add run to the database; (exiting)", exc_info=True)
        raise SystemExit(1)
    logger.debug(
        "\t...added run ID: %s to the database", run_id
    )  # this should use the run_id

    # Identify input files for comparison, and populate the database
    logger.debug("Adding genomes for run %s to database...", run_id)
    try:
        genome_ids = add_run_genomes(
            session, run, args.indir, args.classes, args.labels
        )
    except PyaniORMException:
        logger.error(
            "Could not add genomes to database for run %d; (exiting)",
            run_id,
            exc_info=True,
        )  # this differs from subcmd_anim.py
    logger.debug("\t...added genome IDs: %s", genome_ids)

    # Generate command-liens for fastANI analysis
    logger.info("Generating fastANI command-lines")
    fastanidir = args.outdir / pyani_config.ALIGNDIR["fastANI"]
    logger.debug("fastANI output will be written temporarily to %s", fastanidir)

    # Create output directories. We create the main parent directory (args.outdir), but
    # also subdirectories for the _________________
    logger.debug("Creating output directory %s", fastanidir)
    try:
        fastanidir.mkdir(exist_ok=True, parents=True)
    except IOError:
        logger.error(
            "Could not create output directory %s; (exiting)",
            fastanidir,
            exc_info=True,
        )
        raise SystemError(1)

    # Get list of genomes for this analysis from the database
    logger.info("Compiling genomes for comparison")
    genomes = run.genomes.all()
    logger.debug("\tCollected %d genomes for this run", len(genomes))

    # Generate all pair permutations of genome IDs as a list of (Genome, Geneme) tuples
    logger.info(
        "Compiling pairwise comparisons (this can take time for large datasets)..."
    )
    comparisons = list(permutations(tqdm(genomes, disable=args.disable_tqdm), 2))
    logger.info("\t...total pairwise comparisons to be performed: %d", len(comparisons))

    # Check for existing comparisons; if one has already been done (for the same
    # software package, version, and setting) we add the comparison to this run,
    # but remove it from the list of comparisons to be performed
    logger.info("Checking database for existing comparison data...")
    comparisons_to_run = filter_existing_comparisons(
        session,
        run,
        comparisons,
        "fastANI",
        fastani_version,
        fragsize=args.fragLen,  # fragsize
        maxmatch=False,  # maxmatch
        kmersize=args.kmerSize,
        minmatch=args.minFraction,
    )
    logger.info(
        "\t...after check, still need to run %d comparisons", len(comparisons_to_run)
    )  # this is indicating 0 comparisons — need to diagnose

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
    if args.recovery:
        logger.warning("Entering recovery mode...")
        logger.debug(
            "\tIn this mode, existing comparison output from %s is reused", args.outdir
        )
        existingfiles = collect_existing_output(args.outdir, "fastani", args)
        logger.debug(
            "\tIdentified %d existing output files for reuse", len(existingfiles)
        )
    else:
        existingfiles = None  # in anim this was an empty list; in anib None
        logger.debug("\tIdentified no existing output files")

    # Create list of FastANI jobs for each comparison still to be performed
    logger.info("Creating fastani jobs for fastANI...")
    job = generate_joblist(comparisons_to_run, existingfiles, args)
    logger.debug("...created %d fastani jobs", len(job))

    # Pass jobs to appropriate scheduler
    logger.debug("Passing %d jobs to %s...", len(job), args.scheduler)
    run_fastani_jobs(job, args)
    logger.info("...jobs complete")

    # Process output and add results to database
    # This requires us to drop out of threading/multiprocessing: Python's SQLite3
    # interface doesn't allow sharing connecitons and cursors
    logger.info("Adding comparison results to database...")
    update_comparison_results(job, run, session, fastani_version, args)
    update_comparison_matrices(session, run)
    logger.info("...database updated.")


def generate_joblist(
    comparisons: List,  # in ANIm: List[Tuple]
    existingfiles: List,  # in ANIm: List[Path]
    args: Namespace,
) -> List[ComparisonJob]:
    """Return list of ComparisonJobs

    :param comparisons: list of (Genome, Genome) tuples for which comparisons are needed
    :param existing files: list of pre-existing FastANI outputs
    :param args: Namespace, command-line arguments
    """
    logger = logging.getLogger(__name__)

    joblist = []  # will hold ComparisonJob structs
    for idx, (query, ref) in enumerate(tqdm(comparisons, disable=args.disable_tqdm)):
        # ¶ need to make a thing that produces this
        fastcmd = fastani.construct_fastani_cmdline(
            query.path,
            ref.path,
            args.outdir,
            args.fastani_exe,
            args.fragLen,
            args.kmerSize,
            args.minFraction,
        )
        logger.debug("Commands to run:\n\t%s", fastcmd)
        outfile = fastcmd.split()[6]  # ¶ should this be hard-coded???
        # ¶ There is likely no need to test for different output files
        # ¶ The only exception might be something with --matrix
        outfname = Path(outfile)  # ¶ unsure about this suffix
        logger.debug("Expected output file for db: %s", outfname)

        # If we're in recovery mode, we don't want to repeat a computational
        # comparison that already exists, so we check whether the ultimate
        # output is in the set of existing files and, if not, we add the jobs
        # TODO: something faster than a list search (dict or set?)
        # The comparisons collections always gets updated, so that results are
        # added to the database whether they come from recovery mode or are run
        # in this call of the script.
        if args.recovery and outfname.name in existingfiles:
            logger.debug("Recovering output from %s, not building job", outfname)
            # Need to track the expected output, but set the job itself to None:
            joblist.append(
                ComparisonJob(
                    query,
                    ref,
                    fastcmd,
                    outfname,
                    args.fragLen,
                    args.kmerSize,
                    args.minFraction,
                    None,
                )
            )
        else:
            logger.debug("Building job")
            # Build jobs
            fastjob = pyani_jobs.Job("f{args.jobprefix}_{idx:06d}", fastcmd)
            joblist.append(
                ComparisonJob(
                    query,
                    ref,
                    fastcmd,
                    outfname,
                    args.fragLen,
                    args.kmerSize,
                    args.minFraction,
                    fastjob,
                )
            )
    return joblist


def run_fastani_jobs(joblist: List[ComparisonJob], args: Namespace) -> None:
    """Pass fastANI jobs to the scheduler.

    :param joblist:         list of ComparisonJob namedtuples
    :param args:            command-line arguments for the run
    """
    logger = logging.getLogger(__name__)
    logger.debug("Scheduler: %s", args.scheduler)

    # ¶ Possibly tie this to the --threads option of fastANI?
    if args.scheduler == "multiprocessing":
        logger.info("Running jobs with multiprocessing")
        if not args.workers:
            logger.debug("(using maximum number of worker threads)")
        else:
            logger.debug("(using %d worker threads, if available)", args.workers)
        cumval = run_mp.run_dependency_graph(
            [_.job for _ in joblist], workers=args.workers
        )
        # ¶ Determine if this is needed
        if cumval > 0:
            logger.error(
                "At least one fastANI comparison failed. Please investigate (exiting)"
            )
            raise PyaniException("Multiprocessing run failed in fastANI")
        logger.info("Multiprocessing run completed without error")
    elif args.scheduler.lower() == "sge":
        logger.info("Running jobs with SGE")
        logger.debug("Setting jobarray group size to %d", args.sgegroupsize)
        logger.debug("Joblist contains %d jobs", len(joblist))
        run_sge.run_dependency_graph(
            [_.job for _ in joblist],
            jgprefix=args.jobprefix,
            sgegroupsize=args.sgegroupsize,
            sgeargs=args.sgeargs,
        )
    # elif args.scheduler.lower() == "slurm":
    #     logger.info("Running jobs with Slurm")
    #     logger.debug("Setting jobarray group size to...")
    #     logger.debug("Joblist contains  %d jobs", len(joblist))
    else:
        logger.error(termcolor("Scheduler %s not recognised", "red"), args.scheduler)
        raise SystemError(1)


def update_comparison_results(
    joblist: List[ComparisonJob], run, session, fastani_version: str, args: Namespace
) -> None:
    """Update the Comparison table with the completed result set.

    :param joblist:          list of ComparisonJob namedtuples
    :param run:              Run ORM object for the current ANIm run
    :param session:          active pyanidb session via ORM
    :param fastani_version:  version of fastANI used for the comparison
    :param args:             command-line arguments for this run

    The Comparison table stores individual comparison results, one per row.
    """
    logger = logging.getLogger(__name__)

    # Add individual results to Comparison table
    for job in tqdm(joblist, disable=args.disable_tqdm):
        logger.debug("\t%s vs %s", job.query.description, job.ref.description)
        # ¶ fastANI allows many runs to share the same outfile; might have to alter this
        # ¶ May also need to use something other than this Comparison object
        # ¶ Or add new columns?
        try:
            contents = fastani.parse_fastani_file(job.outfile)
        except fastani.PyaniFastANIException:
            contents = fastani.ComparisonResult(job.query, job.ref, 0, 0, 0)
            print(contents)
        # if len(contents) > 1:
        #     raise ValueError(
        #         "fastANI output file %s has more than one line", job.outfile
        #     )
        logger.debug("Parsed fastANI file contents: %s", contents)
        query, ref, ani, matches, num_frags = contents

        ani = float(ani)  # should be in the range 0–1
        aln_length = matches * job.fragLen  # should be in bp units
        sim_errs = (
            int(num_frags) * job.fragLen - matches * job.fragLen
        )  # should be in bp units
        qcov = (
            float(matches) * job.fragLen / job.query.length
        )  # should be inthe range 0–1
        # try:
        #   pid = (1 - sim_errs) / int(aln_length)
        # except ZeroDivisionError:  # aln_length was zero (no alignment)
        #   pid = 0
        run.comparisons.append(
            Comparison(
                query=job.query,
                subject=job.ref,
                aln_length=int(aln_length),
                sim_errs=int(sim_errs),
                identity=ani,
                cov_query=qcov,
                cov_subject=None,
                program="fastANI",
                version=fastani_version,
                fragsize=job.fragLen,
                maxmatch=False,
                kmersize=job.kmerSize,
                minmatch=job.minFraction,
            )
        )
    # Populate db
    logger.debug("Committing results to database")
    session.commit()
