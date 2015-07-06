#!/usr/bin/env python
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
# (c) The James Hutton Institute 2013-2015
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
# Copyright (c) 2010-2014 The James Hutton Institute
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
import shutil
import sys
import traceback

from argparse import ArgumentParser

from pyani import anib, anim, tetra, pyani_config, pyani_files, pyani_graphics
from pyani.run_multiprocessing import multiprocessing_run
from pyani.pyani_config import params_mpl


# Process command-line arguments
def parse_cmdline(args):
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="average_nucleotide_identity.py")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default=None,
                        help="Output directory")
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default=None,
                        help="Input directory name")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-f", "--force", dest="force",
                        action="store_true", default=False,
                        help="Force file overwriting")
    parser.add_argument("-s", "--fragsize", const="fragsize",
                        action="store_const", default=pyani_config.FRAGSIZE,
                        help="Sequence fragment size for ANIb")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("--skip_nucmer", dest="skip_nucmer",
                        action="store_true", default=False,
                        help="Skip NUCmer runs, for testing " +
                        "(e.g. if output already present)")
    parser.add_argument("--skip_blastn", dest="skip_blastn",
                        action="store_true", default=False,
                        help="Skip BLASTN runs, for testing " +
                        "(e.g. if output already present)")
    parser.add_argument("--noclobber", dest="noclobber",
                        action="store_true", default=False,
                        help="Don't nuke existing files")
    parser.add_argument("-g", "--graphics", dest="graphics",
                        action="store_true", default=False,
                        help="Generate heatmap of ANI")
    parser.add_argument("--gformat", dest="gformat",
                        action="store", default="pdf",
                        help="Graphics output format [pdf|png|jpg|svg]")
    # parser.add_argument("--gmethod", dest="gmethod",
    #                     action="store", default="mpl",
    #                     help="Graphics output method [mpl|R]")
    parser.add_argument("--labels", dest="labels",
                        action="store", default=None,
                        help="Path to file containing sequence labels")
    parser.add_argument("--classes", dest="classes",
                        action="store", default=None,
                        help="Path to file containing sequence classes")
    parser.add_argument("-m", "--method", dest="method",
                        action="store", default="ANIm",
                        help="ANI method [ANIm|ANIb|ANIblastall|TETRA]")
    parser.add_argument("--scheduler", dest="scheduler",
                        action="store", default="multiprocessing",
                        help="Job scheduler [multiprocessing|SGE]")
    parser.add_argument("--maxmatch", dest="maxmatch",
                        action="store_true", default=False,
                        help="Override MUMmer to allow all NUCmer matches")
    parser.add_argument("--nucmer_exe", dest="nucmer_exe",
                        action="store", default=pyani_config.NUCMER_DEFAULT,
                        help="Path to NUCmer executable")
    parser.add_argument("--blastn_exe", dest="blastn_exe",
                        action="store", default=pyani_config.BLASTN_DEFAULT,
                        help="Path to BLASTN+ executable")
    parser.add_argument("--makeblastdb_exe", dest="makeblastdb_exe",
                        action="store",
                        default=pyani_config.MAKEBLASTDB_DEFAULT,
                        help="Path to BLAST+ makeblastdb executable")
    parser.add_argument("--blastall_exe", dest="blastall_exe",
                        action="store", default=pyani_config.BLASTALL_DEFAULT,
                        help="Path to BLASTALL executable")
    parser.add_argument("--formatdb_exe", dest="formatdb_exe",
                        action="store",
                        default=pyani_config.FORMATDB_DEFAULT,
                        help="Path to BLAST formatdb executable")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """ Returns last exception as a string, or use in logging.
    """
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


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
            logger.error("Output directory %s would " % args.outdirname +
                         "overwrite existing files (exiting)")
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything below it" %
                        args.outdirname)
            if args.noclobber:
                logger.warning("NOCLOBBER: not actually deleting directory")
            else:
                shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s" % args.outdirname)
    try:
        os.makedirs(args.outdirname)   # We make the directory recursively
    except OSError:
        # This gets thrown if the directory exists. If we've forced overwrite/
        # delete and we're not clobbering, we let things slide
        if args.noclobber and args.force:
            logger.info("NOCLOBBER+FORCE: not creating directory")
        else:
            logger.error(last_exception)
            sys.exit(1)


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
    # Schedule NUCmer runs
    if not args.skip_nucmer:
        cmdlist = anim.generate_nucmer_commands(infiles, args.outdirname,
                                                nucmer_exe=args.nucmer_exe,
                                                maxmatch=args.maxmatch)
        logger.info("NUCmer commands:\n" + os.linesep.join(cmdlist))
        if args.scheduler == 'multiprocessing':
            logger.info("Running jobs with multiprocessing")
            cumval = multiprocessing_run(cmdlist, verbose=args.verbose)
            logger.info("Cumulative return value: %d" % cumval)
            if 0 < cumval:
                logger.warning("At least one NUCmer comparison failed. " +
                               "ANIm may fail.")
            else:
                logger.info("All multiprocessing jobs complete.")
        else:
            logger.info("Running jobs with SGE")
            raise NotImplementedError
    else:
        logger.warning("Skipping NUCmer run (as instructed)!")

    # Process resulting .delta files
    logger.info("Processing NUCmer .delta files.")
    try:
        data = anim.process_deltadir(args.outdirname, org_lengths)
    except ZeroDivisionError:
        logger.error("One or more NUCmer output files has a problem.")
        if not args.skip_nucmer:
            if 0 < cumval:
                logger.error("This is possibly due to NUCmer run failure, " +
                             "please investigate")
            else:
                logger.error("This is possibly due to a NUCmer comparison " +
                             "being too distant for use. Please consider " +
                             "using the --maxmatch option.")
        logger.error(last_exception())
    return data


# Calculate TETRA for input
def calculate_tetra(infiles, org_lengths):
    """Calculate TETRA for files in input directory.

    - infiles - paths to each input file
    - org_lengths - dictionary of input sequence lengths, keyed by sequence

    Calculates TETRA correlation scores, as described in:

    Richter M, Rossello-Mora R (2009) Shifting the genomic gold standard for the
    prokaryotic species definition. Proc Natl Acad Sci USA 106: 19126-19131.
    doi:10.1073/pnas.0906412106.

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
        logger.info("Calculating TETRA Z-scores for %s" % filename)
        org = os.path.splitext(os.path.split(filename)[-1])[0]
        tetra_zscores[org] = tetra.calculate_tetra_zscore(filename)
    # Then calculate Pearson correlation between Z-scores for each sequence
    logger.info("Calculating TETRA correlation scores.")
    tetra_correlations = tetra.calculate_correlations(tetra_zscores)
    return (tetra_correlations, )


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
    logger.info("Running %s" % args.method)
    # Build BLAST databases and run pairwise BLASTN
    if not args.skip_blastn:
        # Make sequence fragments
        logger.info("Fragmenting input files, and writing to %s" %
                    args.outdirname)
        # Fraglengths does not get reused with BLASTN
        fragfiles, fraglengths = anib.fragment_FASTA_files(infiles,
                                                           args.outdirname,
                                                           args.fragsize)
        # Export fragment lengths as JSON, in case we re-run BLASTALL with
        # --skip_blastn
        if args.method == "ANIblastall":
            with open(os.path.join(args.outdirname,
                                   'fraglengths.json'), 'w') as outfile:
                json.dump(fraglengths, outfile)

        # Which executables are we using?
        if args.method == "ANIblastall":
            blastdb_exe = args.formatdb_exe
            blastn_exe = args.blastall_exe
        else:
            blastdb_exe = args.makeblastdb_exe
            blastn_exe = args.blastn_exe

        # Build BLASTN databases
        logger.info("Constructing %s BLAST databases" % args.method)
        cmdlist = anib.generate_blastdb_commands(infiles, args.outdirname,
                                                 blastdb_exe=blastdb_exe,
                                                 mode=args.method)
        logger.info("Generated commands:\n%s" % '\n'.join(cmdlist))
        if args.scheduler == 'multiprocessing':
            logger.info("Running jobs with multiprocessing")
            cumval = multiprocessing_run(cmdlist, verbose=args.verbose)
            if 0 < cumval:
                logger.warning("At least one makeblastdb run failed. " +
                               "%s may fail." % args.method)
            else:
                logger.info("All multiprocessing jobs complete.")
        else:
            logger.info("Running jobs with SGE")
            raise NotImplementedError

        # Run pairwise BLASTN
        logger.info("Running %s BLASTN jobs" % args.method)
        cmdlist = anib.generate_blastn_commands(fragfiles, args.outdirname,
                                                blastn_exe, mode=args.method)
        logger.info("Generated commands:\n%s" % '\n'.join(cmdlist))
        if args.scheduler == 'multiprocessing':
            logger.info("Running jobs with multiprocessing")
            cumval = multiprocessing_run(cmdlist, verbose=args.verbose)
            logger.info("Cumulative return value: %d" % cumval)
            if 0 < cumval:
                logger.warning("At least one BLASTN comparison failed. " +
                               "%s may fail." % args.method)
            else:
                logger.info("All multiprocessing jobs complete.")
        else:
            logger.info("Running jobs with SGE")
            raise NotImplementedError
    else:
        # Import fragment lengths from JSON
        if args.method == "ANIblastall":
            with open(os.path.join(args.outdirname, 'fraglengths.json'),
                      'rU') as infile:
                fraglengths = json.load(infile)
        else:
            fraglengths = None
        logger.warning("Skipping BLASTN runs (as instructed)!")

    # Process pairwise BLASTN output
    logger.info("Processing pairwise %s BLAST output." % args.method)
    try:
        data = anib.process_blast(args.outdirname, org_lengths,
                                  fraglengths=fraglengths, mode=args.method)
    except ZeroDivisionError:
        logger.error("One or more BLAST output files has a problem.")
        if not args.skip_blastn:
            if 0 < cumval:
                logger.error("This is possibly due to BLASTN run failure, " +
                             "please investigate")
            else:
                logger.error("This is possibly due to a BLASTN comparison " +
                             "being too distant for use.")
        logger.error(last_exception())
    return data


# Read sequence annotations in from file
def get_labels(filename):
    """Returns a dictionary of alternative sequence labels, or None

    - filename - path to file containing tab-separated table of labels

    Input files should be formatted as <key>\t<label>, one pair per line.
    """
    try:
        labeldict = {}
        with open(filename, 'rU') as fh:
            for line in fh.readlines():
                key, label = line.strip().split('\t')
                labeldict[key] = label
        return labeldict
    except:
        None


# Write ANIb/ANIm/TETRA output
def write(results, filestems):
    """Write ANIb/ANIm/TETRA results to output directory.

    - results - tuple of dataframes from analysis

    Each dataframe is written to an Excel-format file, and plain text
    tab-separated file in the output directory. The order of result output
    must be reflected in the order of filestems.
    """
    logger.info("Writing %s results to %s" % (args.method, args.outdirname))
    for df, filestem in zip(results, filestems):
        logger.info("\t%s" % filestem)
        df.to_excel(os.path.join(args.outdirname, filestem) + '.xlsx',
                    index=True)
        df.to_csv(os.path.join(args.outdirname, filestem) + '.tab',
                  index=True, sep="\t")


# Draw ANIb/ANIm/TETRA output
def draw(results, filestems):
    """Draw ANIb/ANIm/TETRA results

    - results - tuple of dataframes from ANIb analysis
    """
    # Draw heatmaps
    for df, filestem in zip(results, filestems):
        fullstem = os.path.join(args.outdirname, filestem)
        outfilename = fullstem + '.%s' % args.gformat
        infilename = fullstem + '.tab'
        logger.info("Writing heatmap to %s" % outfilename)
        # if args.gmethod == "mpl":
        pyani_graphics.heatmap_mpl(df, outfilename=outfilename,
                                   title=filestem,
                                   cmap=params_mpl(df)[filestem][0],
                                   vmin=params_mpl(df)[filestem][1],
                                   vmax=params_mpl(df)[filestem][2],
                                   labels=get_labels(args.labels),
                                   classes=get_labels(args.classes))
        # elif args.gmethod == "R":
        #     rstr = pyani_graphics.heatmap_r(infilename, outfilename,
        #                                     gformat=args.gformat.lower(),
        #                                     title=filestem,
        #                                     cmap=params_r(df)[filestem][0],
        #                                     vmin=params_r(df)[filestem][1],
        #                                     vmax=params_r(df)[filestem][2],
        #                                     labels=get_labels(args.labels),
        #                                     classes=get_labels(args.classes))
        #     logger.info("Executed R code:\n%s" % rstr)


# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # Set up logging
    logger = logging.getLogger('average_nucleotide_identity.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(args)

    # Have we got an input and output directory? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % args.indirname)
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s" % args.outdirname)

    # Have we got a valid method choice?
    # Dictionary below defines analysis function, and output presentation
    # functions/settings, dependent on selected method.
    methods = {"ANIm": (calculate_anim, pyani_config.ANIM_FILESTEMS),
               "ANIb": (unified_anib, pyani_config.ANIB_FILESTEMS),
               "TETRA": (calculate_tetra, pyani_config.TETRA_FILESTEMS),
               "ANIblastall": (unified_anib,
                               pyani_config.ANIBLASTALL_FILESTEMS)}
    if args.method not in methods:
        logger.error("ANI method %s not recognised (exiting)" % args.method)
        logger.error("Valid methods are: %s" % methods.keys())
        sys.exit(1)
    logger.info("Using ANI method: %s" % args.method)

    # Have we got a valid scheduler choice?
    schedulers = ["multiprocessing", "SGE"]
    if args.scheduler not in schedulers:
        logger.error("scheduler %s not recognised (exiting)" % args.scheduler)
        logger.error("Valid schedulers are: %s" % '; '.join(schedulers))
        sys.exit(1)
    logger.info("Using scheduler method: %s" % args.scheduler)

    # Get input files
    logger.info("Identifying FASTA files in %s" % args.indirname)
    infiles = pyani_files.get_fasta_files(args.indirname)
    logger.info("Input files:\n\t%s" % '\n\t'.join(infiles))

    # Get lengths of input sequences
    logger.info("Processing input sequence lengths")
    org_lengths = pyani_files.get_sequence_lengths(infiles)
    logger.info("Sequence lengths:\n" +
                os.linesep.join(["\t%s: %d" % (k, v) for
                                 k, v in org_lengths.items()]))

    # Run appropriate method on the contents of the input directory,
    # and write out corresponding results.
    logger.info("Carrying out %s analysis" % args.method)
    results = methods[args.method][0](infiles, org_lengths)
    write(results, methods[args.method][1])

    # Do we want graphical output?
    if args.graphics:
        logger.info("Rendering output graphics")
        logger.info("Graphics format: %s" % args.gformat)
        # logger.info("Graphics method: %s" % args.gmethod)
        draw(results, methods[args.method][1])

    # Report that we've finished
    logger.info("Done.")
