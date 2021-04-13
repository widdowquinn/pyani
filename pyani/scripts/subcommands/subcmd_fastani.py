# insert copyright info

import datetime
import json  # is this needed?
import logging
import os

from argparse import Namespace  # unsure what this does
from itertools import permutations  # may be specific to anib
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO   # unsure what this does
from tqdm import tqdm   # unsure what this does

from pyani import fastani
from pyani.pyani_files import collect_existing_output
from pyani.pyani_orm import (
    PyaniORMException,
    add_run,
    add_run_genomes,
    filter_existing_comparisons,
    get_session,
    update_comparison_matrices,
)
from pyani.pyani_tools import termcolor

def subcmd_fastani(args: Namespace) -> None:
    """
    """
    logger = logging.getLogger(__name__)
	
    # announce that we're starting
    logger.info(
        termcolor("running FastANI analysis", "red")
    )
	
    # Do we need to get the version of anything? FastANI?
    # fastani_version = fastani.get_version(args.fastani_exe)
    # logger.info(termcolor("FastANI version: %s" % fastani_version, "cyan")
    # the anib version of this is odd 
	
    # Use provided name, or make new one for this analysis
    start_time = datetime.datetime.now()
    name = args.name or "_".join(["FastANI", start_time.isoformat()])
    logger.info(termcolor("Analysis anme: %s" % name, "cyan"))
	
    # Connect to existing database (which may be "clean", or have old analyses)
    logger.debug("Connecting to database %s" % args.dbpath)
    try:
        session = get_session(arge.dbpath)
    except Exception:   # is there a less generic option?
        logger.error(
            "Could not connect to database %s; exiting" %s args.dbpath),
            exc_info = true
        )
        raise SystemExit(1)
	
    # Add information about this run to the database
    logger.debug("Adding run info to database %s..." % args.dbpath)
    try:
        run = add_run(
            session,
            method = "FastANI",
            cmdline = args.cmdline,
            data = start_time,
            status = "started",
            name = name,
        )
    except PyaniORMException:
        logger.error(
            "Could not add run to the database; (exiting)",
            exc_info = True
        )
        raise SystemExit(1)
    logger.debug("\t...added run ID: %s to the database" % run)
    try:
        genome_ids = ad_run_genomes(
            session,
            run,
            args.indir,
            args.classes,
            args.labels
        )
    except PyaniORMException:
        logger.error(
            "Could not add genomes to database for run %s; (exiting)" % run),
            exc_info = True
        )
    logger.debug("\t...added genome IDs: %s" % genome_ids)
    
    # Get list of genomes for this analysis from the database
    logger.info("Compiling genomes for comparison")
    genomes = run.genomes.all()
    logger.debug("\tCollected %s genomes for theis run" % len(genomes))
    
    # Create output directories. We create the main parent directory (args.outdir), but
    # also subdirectories for the _________________
    logger.debug("Creating output directory %s" % args.outdir)
    try:
        os.makedires(args.outdir, exist_ok = True)
    except IOError:
        logger.error(
            "Could not create output directory %s; (exiting)" % args.outdir,
            exc_info = True
        )
        raise SystemError(1)
    
    # Generate all pair permutations of genome IDs as a list of (Genome, Geneme) tuples
    logger.info(
        "Compiling pairwise comparisons (this can take time for large datasets)..."
    )
    comparisons = list(permutations(tqdm(genomes, disable = args.disable_tqdm), 2))
    logger.info("\t...total pairwise comparisons to be performed: %d" % len(comparisons))
    
    # Check for existing comparisons; if one has already been done (for the same
    # software package, version, and setting) we add the comparison to this run,
    # but remove it from the list of comparisons to be performed
    logger.info("Checking database for existing comparison data...")
    comparisons_to_run = filter_existing_comparisons(
        session,
        run,
        comparisons,
        "fastani",
        fastani_version,
        None,   # fragsize
        None    # maxmatch
    )
    logger.info(
        "\t...after check, still need to run %d comparisons" % len(comparisons)
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
        update_coparison_matrices(session, run)
        return
    
    # If we are in recovery mode, we are salvaging output from a previous
    # run, and do not necessarily need to rerun all the jobs. In this case,
    # we prepare a list of output files we want to recover from the results
    # in the output directory.
    if args.recovery:
        logger.warning("Entering recovery mode...")
        logger.debug(
            "\tIn this mode, existing comparison output from %s is reused" % args.outdir
        )
        existingfiles = collect_existing_output(args.outdir, "fastani", args)
        logger.debug(
            "\tIdentified %s existing output files for reuse" % len(existingfiles)
        )
    else:
        existingfiles = None    # in anim this was an empty list; in anib None
        logger.debug("\tIdentified no existing output files")

    # Create list of FastANI jobs for each comparison still to be performed
    logger.info("Creating fastani jobs for fastANI...")  # is this line necessary/accurate?
    joblist = generate_joblist(
        comparisons_to_run, existingfiles, args
    )
    logger.dbug("...created %d fastani jobs" % len(joblist))
    # do we need lines here, like in ANIm, for the scheduler and adding results
    # to the database?
    raise NotImplementedError # needed?

def generate_joblist(
    comparisons: List,      # in ANIm: List[Tuple]
    existingfiles: List,    # in ANIm: List[Path]
    args: Namespace,
) -> NotImplementedError:   # in ANIm: List[ComparisonJob]
    """Return list of ComparisonJobs
    
    :param comparisons: list of (Genome, Genome) tuples for which comparisons are needed
    :param existing files: list of pre-existing FastANI outputs
    :param args: Namespace, command-line arguments
    """
    # logger = logging.getLogger(__name__)  # in ANIm: uncommented
    raise NotImplementedError

# is this needed for fastANI?
def update_comparison_results():
    pass
