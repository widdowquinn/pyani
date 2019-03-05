#!/usr/bin/env python3
#
# average_nucleotide_identity.py
#
# This script calculates Average Nucleotide Identity (ANI) according to one of
# a number of alternative methods described in, e.g.
#
# Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for the
# prokaryotic species definition. Proc Natl Acad Sci USA 106: 19126-19131.
# doi:10.1073/pnas.0906412106. (ANI1020, ANIm, ANIb)
#
# Goris J, Konstantinidis KT, Klappenbach JA, Coenye T, Vandamme P, et al.
# (2007) DNA-DNA hybridization values and their relationship to whole-genome
# sequence similarities. Int J Syst Evol Micr 57: 81-91.
# doi:10.1099/ijs.0.64483-0.
#
# ANI is proposed to be the appropriate in silico substitute for DNA-DNA
# hybridisation (DDH), and so useful for delineating species boundaries. A
# typical percentage threshold for species boundary in the literature is 95%
# ANI (e.g. Richter et al. 2009).
#
# All ANI methods follow the basic algorithm:
# - Align the genome of organism 1 against that of organism 2, and identify
#   the matching regions
# - Calculate the percentage nucleotide identity of the matching regions, as
#   an average for all matching regions
# Methods differ on: (1) what alignment algorithm is used, and the choice of
# parameters (this affects the aligned region boundaries); (2) what the input
# is for alignment (typically either fragments of fixed size, or the most
# complete assembly available); (3) whether a reciprocal comparison is
# necessary or desirable.
#
# ANIm: uses MUMmer (NUCmer) to align the input sequences.
# ANIb: uses BLASTN to align 1000nt fragments of the input sequences
# TETRA: calculates tetranucleotide frequencies of each input sequence
#
# This script takes as main input a directory containing a set of
# correctly-formatted FASTA multiple sequence files. All sequences for a
# single organism should be contained in only one sequence file. The names of
# these files are used for identification, so it would be advisable to name
# them sensibly.
#
# Output is written to a named directory. The output files differ depending on
# the chosen ANI method.
#
# ANIm: MUMmer/NUCmer .delta files, describing the sequence
#       alignment; tab-separated format plain text tables describing total
#       alignment lengths, and total alignment percentage identity
#
# ANIb: FASTA sequences describing 1000nt fragments of each input sequence;
#       BLAST nucleotide databases - one for each set of fragments; and BLASTN
#       output files (tab-separated tabular format plain text) - one for each
#       pairwise comparison of input sequences. There are potentially a lot of
#       intermediate files.
#
# TETRA: Tab-separated text file describing the Z-scores for each
#        tetranucleotide in each input sequence.
#
# In addition, all methods produce a table of output percentage identity (ANIm
# and ANIb) or correlation (TETRA), between each sequence.
#
# If graphical output is chosen, the output directory will also contain PDF
# files representing the similarity between sequences as a heatmap with
# row and column dendrograms.
#
# DEPENDENCIES
# ============
#
# o Biopython (http://www.biopython.org)
#
# o BLAST+ executable in the $PATH, or available on the command line (ANIb)
#       (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
#
# o MUMmer executables in the $PATH, or available on the command line (ANIm)
#       (http://mummer.sourceforge.net/)
#
# For graphical output
# --------------------
#
# o R with shared libraries installed on the system, for graphical output
#       (http://cran.r-project.org/)
#
# o Rpy2 (http://rpy.sourceforge.net/rpy2.html)
#
#
# USAGE
# =====
#
# calculate_ani.py [options]
#
# Options:
#   -h, --help            show this help message and exit
#   -o OUTDIRNAME, --outdir=OUTDIRNAME
#                         Output directory
#   -i INDIRNAME, --indir=INDIRNAME
#                         Input directory name
#   -v, --verbose         Give verbose output
#   -f, --force           Force file overwriting
#   -s, --fragsize        Sequence fragment size for ANIb
#   --skip_nucmer         Skip NUCmer runs, for testing (e.g. if output already
#                         present)
#   --skip_blast          Skip BLAST runs, for testing (e.g. if output already
#                         present)
#   --noclobber           Don't nuke existing files
#   -g, --graphics        Generate heatmap of ANI
#   -m METHOD, --method=METHOD
#                         ANI method
#   --maxmatch            Override MUMmer settings and allow all matches in
#                         NUCmer
#   --nucmer_exe=NUCMER_EXE
#                         Path to NUCmer executable
#   --blast_exe=BLAST_EXE
#                         Path to BLASTN+ executable
#   --makeblastdb_exe=MAKEBLASTDB_EXE
#                         Path to BLAST+ makeblastdb executable
#
# (c) The James Hutton Institute 2013-2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
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
# Copyright (c) 2010-2019 The James Hutton Institute
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

import json
import logging
import logging.handlers
import os
import pandas as pd
import random
import shutil
import sys
import tarfile
import time
import traceback

from argparse import ArgumentParser

from pyani import (
    anib,
    anim,
    tetra,
    pyani_config,
    pyani_files,
    pyani_graphics,
    pyani_tools,
)
from pyani import run_multiprocessing as run_mp
from pyani import run_sge
from pyani.pyani_config import params_mpl, ALIGNDIR, FRAGSIZE, TETRA_FILESTEMS
from pyani import __version__ as VERSION


# Process command-line arguments
def parse_cmdline():
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="average_nucleotide_identity.py")
    parser.add_argument(
        "--version", action="version", version="%(prog)s: pyani " + VERSION
    )
    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdirname",
        action="store",
        default=None,
        required=True,
        help="Output directory (required)",
    )
    parser.add_argument(
        "-i",
        "--indir",
        dest="indirname",
        action="store",
        default=None,
        required=True,
        help="Input directory name (required)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Give verbose output",
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force file overwriting",
    )
    parser.add_argument(
        "-s",
        "--fragsize",
        dest="fragsize",
        action="store",
        default=FRAGSIZE,
        type=int,
        help="Sequence fragment size for ANIb " "(default %i)" % FRAGSIZE,
    )
    parser.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        action="store",
        default=None,
        help="Logfile location",
    )
    parser.add_argument(
        "--skip_nucmer",
        dest="skip_nucmer",
        action="store_true",
        default=False,
        help="Skip NUCmer runs, for testing " + "(e.g. if output already present)",
    )
    parser.add_argument(
        "--skip_blastn",
        dest="skip_blastn",
        action="store_true",
        default=False,
        help="Skip BLASTN runs, for testing " + "(e.g. if output already present)",
    )
    parser.add_argument(
        "--noclobber",
        dest="noclobber",
        action="store_true",
        default=False,
        help="Don't nuke existing files",
    )
    parser.add_argument(
        "--nocompress",
        dest="nocompress",
        action="store_true",
        default=False,
        help="Don't compress/delete the comparison output",
    )
    parser.add_argument(
        "-g",
        "--graphics",
        dest="graphics",
        action="store_true",
        default=False,
        help="Generate heatmap of ANI",
    )
    parser.add_argument(
        "--gformat",
        dest="gformat",
        action="store",
        default="pdf,png,eps",
        help="Graphics output format(s) [pdf|png|jpg|svg] "
        "(default pdf,png,eps meaning three file formats)",
    )
    parser.add_argument(
        "--gmethod",
        dest="gmethod",
        action="store",
        default="mpl",
        choices=["mpl", "seaborn"],
        help="Graphics output method (default mpl)",
    )
    parser.add_argument(
        "--labels",
        dest="labels",
        action="store",
        default=None,
        help="Path to file containing sequence labels",
    )
    parser.add_argument(
        "--classes",
        dest="classes",
        action="store",
        default=None,
        help="Path to file containing sequence classes",
    )
    parser.add_argument(
        "-m",
        "--method",
        dest="method",
        action="store",
        default="ANIm",
        choices=["ANIm", "ANIb", "ANIblastall", "TETRA"],
        help="ANI method (default ANIm)",
    )
    parser.add_argument(
        "--scheduler",
        dest="scheduler",
        action="store",
        default="multiprocessing",
        choices=["multiprocessing", "SGE"],
        help="Job scheduler (default multiprocessing, i.e. locally)",
    )
    parser.add_argument(
        "--workers",
        dest="workers",
        action="store",
        default=None,
        type=int,
        help="Number of worker processes for multiprocessing "
        "(default zero, meaning use all available cores)",
    )
    parser.add_argument(
        "--SGEgroupsize",
        dest="sgegroupsize",
        action="store",
        default=10000,
        type=int,
        help="Number of jobs to place in an SGE array group " "(default 10000)",
    )
    parser.add_argument(
        "--SGEargs",
        dest="sgeargs",
        action="store",
        default=None,
        type=str,
        help="Additional arguments for qsub",
    )
    parser.add_argument(
        "--maxmatch",
        dest="maxmatch",
        action="store_true",
        default=False,
        help="Override MUMmer to allow all NUCmer matches",
    )
    parser.add_argument(
        "--nucmer_exe",
        dest="nucmer_exe",
        action="store",
        default=pyani_config.NUCMER_DEFAULT,
        help="Path to NUCmer executable",
    )
    parser.add_argument(
        "--filter_exe",
        dest="filter_exe",
        action="store",
        default=pyani_config.FILTER_DEFAULT,
        help="Path to delta-filter executable",
    )
    parser.add_argument(
        "--blastn_exe",
        dest="blastn_exe",
        action="store",
        default=pyani_config.BLASTN_DEFAULT,
        help="Path to BLASTN+ executable",
    )
    parser.add_argument(
        "--makeblastdb_exe",
        dest="makeblastdb_exe",
        action="store",
        default=pyani_config.MAKEBLASTDB_DEFAULT,
        help="Path to BLAST+ makeblastdb executable",
    )
    parser.add_argument(
        "--blastall_exe",
        dest="blastall_exe",
        action="store",
        default=pyani_config.BLASTALL_DEFAULT,
        help="Path to BLASTALL executable",
    )
    parser.add_argument(
        "--formatdb_exe",
        dest="formatdb_exe",
        action="store",
        default=pyani_config.FORMATDB_DEFAULT,
        help="Path to BLAST formatdb executable",
    )
    parser.add_argument(
        "--write_excel",
        dest="write_excel",
        action="store_true",
        default=False,
        help="Write Excel format output tables",
    )
    parser.add_argument(
        "--rerender",
        dest="rerender",
        action="store_true",
        default=False,
        help="Rerender graphics output without recalculation",
    )
    parser.add_argument(
        "--subsample",
        dest="subsample",
        action="store",
        default=None,
        help="Subsample a percentage [0-1] or specific "
        + "number (1-n) of input sequences",
    )
    parser.add_argument(
        "--seed",
        dest="seed",
        action="store",
        default=None,
        help="Set random seed for reproducible subsampling.",
    )
    parser.add_argument(
        "--jobprefix",
        dest="jobprefix",
        action="store",
        default="ANI",
        help="Prefix for SGE jobs (default ANI).",
    )
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))


# Create output directory if it doesn't exist
def make_outdir():
    """Make the output directory, if required.

    This is a little involved.  If the output directory already exists,
    we take the safe option by default, and stop with an error.  We can,
    however, choose to force the program to go on, in which case we can
    either clobber the existing directory, or not.  The options turn out
    as the following, if the directory exists:

    DEFAULT: stop and report the collision
    FORCE: continue, and remove the existing output directory
    NOCLOBBER+FORCE: continue, but do not remove the existing output
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error(
                "Output directory %s would overwrite existing " + "files (exiting)",
                args.outdirname,
            )
            sys.exit(1)
        elif args.noclobber:
            logger.warning(
                "NOCLOBBER: not actually deleting directory %s", args.outdirname
            )
        else:
            logger.info(
                "Removing directory %s and everything below it", args.outdirname
            )
            shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s", args.outdirname)
    try:
        os.makedirs(args.outdirname)  # We make the directory recursively
        # Depending on the choice of method, a subdirectory will be made for
        # alignment output files
        if args.method != "TETRA":
            os.makedirs(os.path.join(args.outdirname, ALIGNDIR[args.method]))
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)


# Compress output directory and delete it
def compress_delete_outdir(outdir):
    """Compress the contents of the passed directory to .tar.gz and delete."""
    # Compress output in .tar.gz file and remove raw output
    tarfn = outdir + ".tar.gz"
    logger.info("\tCompressing output from %s to %s", outdir, tarfn)
    with tarfile.open(tarfn, "w:gz") as fh:
        fh.add(outdir)
    logger.info("\tRemoving output directory %s", outdir)
    shutil.rmtree(outdir)


# Calculate ANIm for input
def calculate_anim(infiles, org_lengths):
    """Returns ANIm result dataframes for files in input directory.

    - infiles - paths to each input file
    - org_lengths - dictionary of input sequence lengths, keyed by sequence

    Finds ANI by the ANIm method, as described in Richter et al (2009)
    Proc Natl Acad Sci USA 106: 19126-19131 doi:10.1073/pnas.0906412106.

    All FASTA format files (selected by suffix) in the input directory
    are compared against each other, pairwise, using NUCmer (which must
    be in the path). NUCmer output is stored in the output directory.

    The NUCmer .delta file output is parsed to obtain an alignment length
    and similarity error count for every unique region alignment between
    the two organisms, as represented by the sequences in the FASTA files.

    These are processed to give matrices of aligned sequence lengths,
    average nucleotide identity (ANI) percentages, coverage (aligned
    percentage of whole genome), and similarity error cound for each pairwise
    comparison.
    """
    logger.info("Running ANIm")
    logger.info("Generating NUCmer command-lines")
    deltadir = os.path.join(args.outdirname, ALIGNDIR["ANIm"])
    logger.info("Writing nucmer output to %s", deltadir)
    # Schedule NUCmer runs
    if not args.skip_nucmer:
        joblist = anim.generate_nucmer_jobs(
            infiles,
            args.outdirname,
            nucmer_exe=args.nucmer_exe,
            filter_exe=args.filter_exe,
            maxmatch=args.maxmatch,
            jobprefix=args.jobprefix,
        )
        if args.scheduler == "multiprocessing":
            logger.info("Running jobs with multiprocessing")
            if args.workers is None:
                logger.info("(using maximum number of available " + "worker threads)")
            else:
                logger.info("(using %d worker threads, if available)", args.workers)
            cumval = run_mp.run_dependency_graph(
                joblist, workers=args.workers, logger=logger
            )
            logger.info("Cumulative return value: %d", cumval)
            if 0 < cumval:
                logger.warning(
                    "At least one NUCmer comparison failed. " + "ANIm may fail."
                )
            else:
                logger.info("All multiprocessing jobs complete.")
        else:
            logger.info("Running jobs with SGE")
            logger.info("Jobarray group size set to %d", args.sgegroupsize)
            run_sge.run_dependency_graph(
                joblist,
                logger=logger,
                jgprefix=args.jobprefix,
                sgegroupsize=args.sgegroupsize,
                sgeargs=args.sgeargs,
            )
    else:
        logger.warning("Skipping NUCmer run (as instructed)!")

    # Process resulting .delta files
    logger.info("Processing NUCmer .delta files.")
    results = anim.process_deltadir(deltadir, org_lengths, logger=logger)
    if results.zero_error:  # zero percentage identity error
        if not args.skip_nucmer and args.scheduler == "multiprocessing":
            if 0 < cumval:
                logger.error(
                    "This has possibly been a NUCmer run failure, "
                    + "please investigate"
                )
                logger.error(last_exception())
                sys.exit(1)
            else:
                logger.error(
                    "This is possibly due to a NUCmer comparison "
                    + "being too distant for use. Please consider "
                    + "using the --maxmatch option."
                )
                logger.error(
                    "This is alternatively due to NUCmer run "
                    + "failure, analysis will continue, but please "
                    + "investigate."
                )
    if not args.nocompress:
        logger.info("Compressing/deleting %s", deltadir)
        compress_delete_outdir(deltadir)

    # Return processed data from .delta files
    return results


# Calculate TETRA for input
def calculate_tetra(infiles):
    """Calculate TETRA for files in input directory.

    - infiles - paths to each input file
    - org_lengths - dictionary of input sequence lengths, keyed by sequence

    Calculates TETRA correlation scores, as described in:

    Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for
    the prokaryotic species definition. Proc Natl Acad Sci USA 106:
    19126-19131. doi:10.1073/pnas.0906412106.

    and

    Teeling et al. (2004) Application of tetranucleotide frequencies for the
    assignment of genomic fragments. Env. Microbiol. 6(9): 938-947.
    doi:10.1111/j.1462-2920.2004.00624.x
    """
    logger.info("Running TETRA.")
    # First, find Z-scores
    logger.info("Calculating TETRA Z-scores for each sequence.")
    tetra_zscores = {}
    for filename in infiles:
        logger.info("Calculating TETRA Z-scores for %s", filename)
        org = os.path.splitext(os.path.split(filename)[-1])[0]
        tetra_zscores[org] = tetra.calculate_tetra_zscore(filename)
    # Then calculate Pearson correlation between Z-scores for each sequence
    logger.info("Calculating TETRA correlation scores.")
    tetra_correlations = tetra.calculate_correlations(tetra_zscores)
    return tetra_correlations


# Calculate ANIb for input
def unified_anib(infiles, org_lengths):
    """Calculate ANIb for files in input directory.

    - infiles - paths to each input file
    - org_lengths - dictionary of input sequence lengths, keyed by sequence

    Calculates ANI by the ANIb method, as described in Goris et al. (2007)
    Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0. There are
    some minor differences depending on whether BLAST+ or legacy BLAST
    (BLASTALL) methods are used.

    All FASTA format files (selected by suffix) in the input directory are
    used to construct BLAST databases, placed in the output directory.
    Each file's contents are also split into sequence fragments of length
    options.fragsize, and the multiple FASTA file that results written to
    the output directory. These are BLASTNed, pairwise, against the
    databases.

    The BLAST output is interrogated for all fragment matches that cover
    at least 70% of the query sequence, with at least 30% nucleotide
    identity over the full length of the query sequence. This is an odd
    choice and doesn't correspond to the twilight zone limit as implied by
    Goris et al. We persist with their definition, however.  Only these
    qualifying matches contribute to the total aligned length, and total
    aligned sequence identity used to calculate ANI.

    The results are processed to give matrices of aligned sequence length
    (aln_lengths.tab), similarity error counts (sim_errors.tab), ANIs
    (perc_ids.tab), and minimum aligned percentage (perc_aln.tab) of
    each genome, for each pairwise comparison. These are written to the
    output directory in plain text tab-separated format.
    """
    logger.info("Running %s", args.method)
    blastdir = os.path.join(args.outdirname, ALIGNDIR[args.method])
    logger.info("Writing BLAST output to %s", blastdir)
    # Build BLAST databases and run pairwise BLASTN
    if not args.skip_blastn:
        # Make sequence fragments
        logger.info("Fragmenting input files, and writing to %s", args.outdirname)
        # Fraglengths does not get reused with BLASTN
        fragfiles, fraglengths = anib.fragment_fasta_files(
            infiles, blastdir, args.fragsize
        )
        # Export fragment lengths as JSON, in case we re-run with --skip_blastn
        with open(os.path.join(blastdir, "fraglengths.json"), "w") as outfile:
            json.dump(fraglengths, outfile)

        # Which executables are we using?
        # if args.method == "ANIblastall":
        #    format_exe = args.formatdb_exe
        #    blast_exe = args.blastall_exe
        # else:
        #    format_exe = args.makeblastdb_exe
        #    blast_exe = args.blastn_exe

        # Run BLAST database-building and executables from a jobgraph
        logger.info("Creating job dependency graph")
        jobgraph = anib.make_job_graph(
            infiles, fragfiles, anib.make_blastcmd_builder(args.method, blastdir)
        )
        # jobgraph = anib.make_job_graph(infiles, fragfiles, blastdir,
        #                               format_exe, blast_exe, args.method,
        #                               jobprefix=args.jobprefix)
        if args.scheduler == "multiprocessing":
            logger.info("Running jobs with multiprocessing")
            logger.info("Running job dependency graph")
            if args.workers is None:
                logger.info("(using maximum number of available " + "worker threads)")
            else:
                logger.info("(using %d worker threads, if available)", args.workers)
            cumval = run_mp.run_dependency_graph(
                jobgraph, workers=args.workers, logger=logger
            )
            if 0 < cumval:
                logger.warning(
                    "At least one BLAST run failed. " + "%s may fail.", args.method
                )
            else:
                logger.info("All multiprocessing jobs complete.")
        else:
            run_sge.run_dependency_graph(jobgraph, logger=logger)
            logger.info("Running jobs with SGE")
    else:
        # Import fragment lengths from JSON
        if args.method == "ANIblastall":
            with open(os.path.join(blastdir, "fraglengths.json"), "rU") as infile:
                fraglengths = json.load(infile)
        else:
            fraglengths = None
        logger.warning("Skipping BLASTN runs (as instructed)!")

    # Process pairwise BLASTN output
    logger.info("Processing pairwise %s BLAST output.", args.method)
    try:
        data = anib.process_blast(
            blastdir, org_lengths, fraglengths=fraglengths, mode=args.method
        )
    except ZeroDivisionError:
        logger.error("One or more BLAST output files has a problem.")
        if not args.skip_blastn:
            if 0 < cumval:
                logger.error(
                    "This is possibly due to BLASTN run failure, "
                    + "please investigate"
                )
            else:
                logger.error(
                    "This is possibly due to a BLASTN comparison "
                    + "being too distant for use."
                )
        logger.error(last_exception())
    if not args.nocompress:
        logger.info("Compressing/deleting %s", blastdir)
        compress_delete_outdir(blastdir)

    # Return processed BLAST data
    return data


# Write ANIb/ANIm/TETRA output
def write(results):
    """Write ANIb/ANIm/TETRA results to output directory.

    - results - results object from analysis

    Each dataframe is written to an Excel-format file (if args.write_excel is
    True), and plain text tab-separated file in the output directory. The
    order of result output must be reflected in the order of filestems.
    """
    logger.info("Writing %s results to %s", args.method, args.outdirname)
    if args.method == "TETRA":
        out_excel = os.path.join(args.outdirname, TETRA_FILESTEMS[0]) + ".xlsx"
        out_csv = os.path.join(args.outdirname, TETRA_FILESTEMS[0]) + ".tab"
        if args.write_excel:
            results.to_excel(out_excel, index=True)
        results.to_csv(out_csv, index=True, sep="\t")

    else:
        for dfr, filestem in results.data:
            out_excel = os.path.join(args.outdirname, filestem) + ".xlsx"
            out_csv = os.path.join(args.outdirname, filestem) + ".tab"
            logger.info("\t%s", filestem)
            if args.write_excel:
                dfr.to_excel(out_excel, index=True)
            dfr.to_csv(out_csv, index=True, sep="\t")


# Draw ANIb/ANIm/TETRA output
def draw(filestems, gformat):
    """Draw ANIb/ANIm/TETRA results

    - filestems - filestems for output files
    - gformat - the format for output graphics
    """
    # Draw heatmaps
    for filestem in filestems:
        fullstem = os.path.join(args.outdirname, filestem)
        outfilename = fullstem + ".%s" % gformat
        infilename = fullstem + ".tab"
        df = pd.read_csv(infilename, index_col=0, sep="\t")
        logger.info("Writing heatmap to %s", outfilename)
        params = pyani_graphics.Params(
            params_mpl(df)[filestem],
            pyani_tools.get_labels(args.labels),
            pyani_tools.get_labels(args.classes),
        )
        if args.gmethod == "mpl":
            pyani_graphics.heatmap_mpl(
                df, outfilename=outfilename, title=filestem, params=params
            )
        elif args.gmethod == "seaborn":
            pyani_graphics.heatmap_seaborn(
                df, outfilename=outfilename, title=filestem, params=params
            )


# Subsample the input files
def subsample_input(infiles):
    """Returns a random subsample of the input files.

    - infiles: a list of input files for analysis
    """
    logger.info("--subsample: %s", args.subsample)
    try:
        samplesize = float(args.subsample)
    except TypeError:  # Not a number
        logger.error(
            "--subsample must be int or float, got %s (exiting)", type(args.subsample)
        )
        sys.exit(1)
    if samplesize <= 0:  # Not a positive value
        logger.error("--subsample must be positive value, got %s", str(args.subsample))
        sys.exit(1)
    if int(samplesize) > 1:
        logger.info("Sample size integer > 1: %d", samplesize)
        k = min(int(samplesize), len(infiles))
    else:
        logger.info("Sample size proportion in (0, 1]: %.3f", samplesize)
        k = int(min(samplesize, 1.0) * len(infiles))
    logger.info("Randomly subsampling %d sequences for analysis", k)
    if args.seed:
        logger.info("Setting random seed with: %s", args.seed)
        random.seed(args.seed)
    else:
        logger.warning("Subsampling without specified random seed!")
        logger.warning("Subsampling may NOT be easily reproducible!")
    return random.sample(infiles, k)


# Run as script
if __name__ == "__main__":

    # Parse command-line
    args = parse_cmdline()

    # Set up logging
    logger = logging.getLogger("average_nucleotide_identity.py: %s" % time.asctime())
    t0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter("%(levelname)s: %(message)s")
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, "w")
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging", args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info("pyani version: %s", VERSION)
    logger.info(args)
    logger.info("command-line: %s", " ".join(sys.argv))

    # Have we got an input and output directory? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s", args.indirname)
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    if args.rerender:  # Rerendering, we want to overwrite graphics
        args.force, args.noclobber = True, True
    make_outdir()
    logger.info("Output directory: %s", args.outdirname)

    # Check for the presence of space characters in any of the input filenames
    # or output directory. If we have any, abort here and now.
    filenames = [args.outdirname] + os.listdir(args.indirname)
    for fname in filenames:
        if " " in os.path.abspath(fname):
            logger.error("File or directory '%s' contains whitespace", fname)
            logger.error("This will cause issues with MUMmer and BLAST")
            logger.error("(exiting)")
            sys.exit(1)

    if args.labels and not os.path.isfile(args.labels):
        logger.error("Missing labels file: %s", args.labels)
        sys.exit(1)
    if args.classes and not os.path.isfile(args.classes):
        logger.error("Missing classes file: %s", args.classes)
        sys.exit(1)

    # Have we got a valid method choice?
    # Dictionary below defines analysis function, and output presentation
    # functions/settings, dependent on selected method.
    methods = {
        "ANIm": (calculate_anim, pyani_config.ANIM_FILESTEMS),
        "ANIb": (unified_anib, pyani_config.ANIB_FILESTEMS),
        "TETRA": (calculate_tetra, pyani_config.TETRA_FILESTEMS),
        "ANIblastall": (unified_anib, pyani_config.ANIBLASTALL_FILESTEMS),
    }
    if args.method not in methods:
        logger.error("ANI method %s not recognised (exiting)", args.method)
        logger.error("Valid methods are: %s", list(methods.keys()))
        sys.exit(1)
    logger.info("Using ANI method: %s", args.method)

    # Skip calculations (or not) depending on rerender option
    if args.rerender:
        logger.warning("--rerender option used")
        logger.warning("Producing graphics with no new recalculations")
    else:
        # Have we got a valid scheduler choice?
        schedulers = ["multiprocessing", "SGE"]
        if args.scheduler not in schedulers:
            logger.error("scheduler %s not recognised (exiting)", args.scheduler)
            logger.error("Valid schedulers are: %s", "; ".join(schedulers))
            sys.exit(1)
        logger.info("Using scheduler method: %s", args.scheduler)

        # Get input files
        logger.info("Identifying FASTA files in %s", args.indirname)
        infiles = pyani_files.get_fasta_files(args.indirname)
        logger.info("Input files:\n\t%s", "\n\t".join(infiles))

        # Are we subsampling? If so, make the selection here
        if args.subsample:
            infiles = subsample_input(infiles)
            logger.info("Sampled input files:\n\t%s", "\n\t".join(infiles))

        # Get lengths of input sequences
        logger.info("Processing input sequence lengths")
        org_lengths = pyani_files.get_sequence_lengths(infiles)
        logger.info(
            "Sequence lengths:\n"
            + os.linesep.join(
                ["\t%s: %d" % (k, v) for k, v in list(org_lengths.items())]
            )
        )

        # Run appropriate method on the contents of the input directory,
        # and write out corresponding results.
        logger.info("Carrying out %s analysis", args.method)
        if args.method == "TETRA":
            results = methods[args.method][0](infiles)
        else:
            results = methods[args.method][0](infiles, org_lengths)
        write(results)

    # Do we want graphical output?
    if args.graphics or args.rerender:
        logger.info("Rendering output graphics")
        logger.info("Formats requested: %s", args.gformat)
        for gfmt in args.gformat.split(","):
            logger.info("Graphics format: %s", gfmt)
            logger.info("Graphics method: %s", args.gmethod)
            draw(methods[args.method][1], gfmt)

    # Report that we've finished
    logger.info("Done: %s.", time.asctime())
    logger.info("Time taken: %.2fs", (time.time() - t0))
