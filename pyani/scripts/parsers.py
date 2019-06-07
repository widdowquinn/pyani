# -*- coding: utf-8 -*-
"""parsers.py

Provides command-line/subcommand parsers for the pyani.py script.

- createdb:        create local database for pyani analyses
- download:        download genome files from NCBI
- index:           index local genome files
- anib:            perform ANIb analysis
- aniblastall:     perform ANIb analysis with legacy BLAST
- anim:            perform ANIm analysis
- classify:        classify input sequences on the basis of ANI output
- report:          report analysis tables from database
- plot:            report analysis as images

(c) The James Hutton Institute 2016-2019
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

Leighton Pritchard,
Information and Computing Sciences,
James Hutton Institute,
Errol Road,
Invergowrie,
Dundee,
DD2 5DA,
Scotland,
UK

The MIT License

Copyright (c) 2016-2019 The James Hutton Institute

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

import sys

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from . import subcommands
from .. import pyani_config


# Common parser for all subcommands
def build_parser_common():
    """Return the common argument parser for all script subcommands.

    Common arguments are:

    -l, --logfile (specify logfile output path)
    -v, --verbose (produce verbose output on STDOUT)
    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "-l",
        "--logfile",
        dest="logfile",
        action="store",
        default=None,
        help="logfile location",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="report verbose progress to log",
    )
    parser.add_argument(
        "--disable_tqdm",
        action="store_true",
        dest="disable_tqdm",
        default=False,
        help="Turn off tqdm progress bar",
    )
    return parser


# Common parser for all subcommands
def build_parser_scheduler():
    """Return the common argument parser for job scheduling.

    Common arguments are:

    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "--scheduler",
        dest="scheduler",
        action="store",
        default="multiprocessing",
        choices=["multiprocessing", "SGE"],
        help="Job scheduler (default multiprocessing, " + "i.e. locally)",
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
        "--jobprefix",
        dest="jobprefix",
        action="store",
        default="PYANI",
        help="Prefix for SGE jobs (default PYANI).",
    )
    return parser


# Subcommand parsers
def build_parser_download(subps, parents=None):
    """Return a command-line parser for the download subcommand.

    The download subcommand takes specific arguments:

    -t, --taxon (NCBI taxonomy IDs - comma-separated list, or one ID)
    --email     (email for providing to Entrez services)
    --retries   (number of Entrez retry attempts to make)
    --batchsize (number of Entrez records to download in a batch)
    --timeout   (how long to wait for Entrez query timeout)
    -f, --force (allow existing directory overwrite)
    --noclobber (don't replace existing files)
    --labels    (path to write labels file)
    --classes   (path to write classes file)
    """
    parser = subps.add_parser(
        "download", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional argument: output directory
    parser.add_argument(
        action="store", dest="outdir", default=None, help="output directory"
    )
    # Required arguments for NCBI query
    parser.add_argument(
        "-t",
        "--taxon",
        dest="taxon",
        action="store",
        default=None,
        help="NCBI taxonomy IDs (required, " + "comma-separated list)",
        required=True,
    )
    parser.add_argument(
        "--email",
        dest="email",
        action="store",
        default=None,
        help="Email associated with NCBI queries (required)",
        required=True,
    )
    # Arguments controlling connection to NCBI for download
    parser.add_argument(
        "--retries",
        dest="retries",
        action="store",
        default=20,
        help="Number of Entrez retry attempts per request",
    )
    parser.add_argument(
        "--batchsize",
        dest="batchsize",
        action="store",
        default=10000,
        help="Entrez record return batch size",
    )
    parser.add_argument(
        "--timeout",
        dest="timeout",
        action="store",
        default=10,
        help="Timeout for URL connection (s)",
    )
    # Arguments controlling local filehandling
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Allow download to existing directory",
    )
    parser.add_argument(
        "--noclobber",
        dest="noclobber",
        action="store_true",
        default=False,
        help="Don't replace existing files",
    )
    # Names for output files
    parser.add_argument(
        "--labels",
        dest="labelfname",
        action="store",
        default="labels.txt",
        help="Filename for labels file",
    )
    parser.add_argument(
        "--classes",
        dest="classfname",
        action="store",
        default="classes.txt",
        help="Filename for classes file",
    )
    # Output for Kraken
    parser.add_argument(
        "--kraken",
        dest="kraken",
        action="store_true",
        default=False,
        help="Modify downloaded sequence ID for Kraken",
    )
    # Dry-run: do everything except download
    parser.add_argument(
        "--dry-run",
        dest="dryrun",
        action="store_true",
        default=False,
        help="Dry run only, do not download or overwrite.",
    )
    parser.set_defaults(func=subcommands.subcmd_download)


def build_parser_classify(subps, parents=None):
    """Return a command-line parser for the classify subcommand.

    The classify subcommand takes specific arguments:

    --cov_min       (minimum coverage threshold for an edge)
    --id_min        (minimum identity threshold for an edge)
    --resolution    (number of identity thresholds to test)
    """
    parser = subps.add_parser(
        "classify", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional arguments: output directory and run ID
    parser.add_argument(
        action="store", dest="outdir", default=None, help="output directory"
    )
    parser.add_argument(
        action="store", dest="run_id", default=None, help="run ID to classify"
    )
    # Other optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    # Parameters for classification: minimum %coverage, %identity,
    # and the resolution of thresholds to test
    parser.add_argument(
        "--cov_min",
        dest="cov_min",
        action="store",
        type=float,
        default=0.5,
        help="minimum %%coverage for an edge",
    )
    parser.add_argument(
        "--id_min",
        dest="id_min",
        action="store",
        type=float,
        default=0.8,
        help="minimum %%identity for an edge",
    )
    parser.add_argument(
        "--resolution",
        dest="resolution",
        action="store",
        type=int,
        default=1500,
        help="number of identity thresholds to test",
    )
    parser.set_defaults(func=subcommands.subcmd_classify)


def build_parser_plot(subps, parents=None):
    """Return a command-line parser for the plot subcommand.

    The plot subcommand takes specific arguments:

    --method        (graphics method to use)
    """
    parser = subps.add_parser(
        "plot", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional arguments: output directory and run ID
    parser.add_argument(
        action="store", dest="outdir", default=None, help="output directory"
    )
    parser.add_argument(
        action="store", dest="run_id", default=None, help="run ID to plot"
    )
    # Other optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    # Graphics methods and formats
    parser.add_argument(
        "--formats",
        dest="formats",
        action="store",
        default="png",
        help="graphics output format (seaborn/mpl)",
    )
    parser.add_argument(
        "--method",
        dest="method",
        action="store",
        default="seaborn",
        help="graphics method to use for plotting",
        choices=["seaborn", "mpl", "plotly"],
    )
    parser.set_defaults(func=subcommands.subcmd_plot)


def build_parser_index(subps, parents=None):
    """Return a command-line parser for the index subcommand.

    The index subcommand takes a single positional argument:

    indir           (directory containing input genome sequence files)
    """
    parser = subps.add_parser(
        "index", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional argument: input directory
    parser.add_argument(
        action="store", dest="indir", default=None, help="input directory"
    )
    # Names for output files
    parser.add_argument(
        "--labels",
        dest="labelfname",
        action="store",
        default="labels.txt",
        help="Filename for labels file",
    )
    parser.add_argument(
        "--classes",
        dest="classfname",
        action="store",
        default="classes.txt",
        help="Filename for classes file",
    )
    parser.set_defaults(func=subcommands.subcmd_index)


def build_parser_createdb(subps, parents=None):
    """Return a command-line parser for the createdb subcommand."""
    parser = subps.add_parser(
        "createdb", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Path to database (default: .pyani/pyanidb)
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        dest="force",
        default=False,
        help="force creation of new empty database",
    )
    parser.set_defaults(func=subcommands.subcmd_createdb)


def build_parser_db(subps, parents=None):
    """Return a command-line parser for the db subcommand."""
    parser = subps.add_parser(
        "db", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        action="store",
        type=int,
        dest="run_id",
        default=None,
        help="work with this runID",
    )
    # Path to database (default: .pyani/pyanidb)
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    parser.add_argument(
        "--relabel",
        action="store",
        dest="relabelfname",
        default=None,
        help="relabel genomes with hashes in passed file",
    )
    parser.add_argument(
        "--reclass",
        action="store",
        dest="reclassfname",
        default=None,
        help="reclass genomes with hashes in passed file",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        dest="force",
        default=False,
        help="force database operation",
    )
    parser.set_defaults(func=subcommands.subcmd_db)


# Common parser for all subcommands
def build_parser_run_common():
    """Return the common argument parser for analysis run subcommands.

    Common arguments are:

    --name     (human-readable name for the run)
    --labels   (genome labels for this run)
    --classes  (genome classes for this run)
    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "--name",
        dest="name",
        action="store",
        default=None,
        help="Provide a name to help identify the analysis run",
    )
    parser.add_argument(
        "--classes",
        dest="classes",
        action="store",
        default=None,
        help="Path to genome classes for the run",
    )
    parser.add_argument(
        "--labels",
        dest="labels",
        action="store",
        default=None,
        help="Path to genome labels for the run",
    )
    parser.add_argument(
        "--recovery",
        dest="recovery",
        action="store_true",
        default=False,
        help="Recover comparison results for this run",
    )
    return parser


def build_parser_anim(subps, parents=None):
    """Return a command-line parser for the anim subcommand."""
    parser = subps.add_parser(
        "anim", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional arguments: input and output directories
    parser.add_argument(
        action="store", dest="indir", default=None, help="input genome directory"
    )
    parser.add_argument(
        action="store",
        dest="outdir",
        default=None,
        help="output analysis results directory",
    )
    # Optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    parser.add_argument(
        "--nucmer_exe",
        dest="nucmer_exe",
        action="store",
        default=pyani_config.NUCMER_DEFAULT,
        help="path to NUCmer executable",
    )
    parser.add_argument(
        "--filter_exe",
        dest="filter_exe",
        action="store",
        default=pyani_config.FILTER_DEFAULT,
        help="path to delta-filter executable",
    )
    parser.add_argument(
        "--maxmatch",
        dest="maxmatch",
        action="store_true",
        default=False,
        help="override MUMmer to allow all NUCmer matches",
    )
    parser.add_argument(
        "--nofilter",
        dest="nofilter",
        action="store_true",
        default=False,
        help="do not use delta-filter for 1:1 NUCmer matches",
    )
    parser.set_defaults(func=subcommands.subcmd_anim)


def build_parser_anib(subps, parents=None):
    """Return a command-line parser for the anib subcommand."""
    parser = subps.add_parser(
        "anib", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.set_defaults(func=subcommands.subcmd_anib)


def build_parser_aniblastall(subps, parents=None):
    """Return a command-line parser for the aniblastall subcommand."""
    parser = subps.add_parser(
        "aniblastall", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.set_defaults(func=subcommands.subcmd_aniblastall)


def build_parser_report(subps, parents=None):
    """Return a command-line parser for the report subcommand."""
    parser = subps.add_parser(
        "report", parents=parents, formatter_class=ArgumentDefaultsHelpFormatter
    )
    # Required positional argument: output directory
    parser.add_argument(
        action="store",
        dest="outdir",
        default=None,
        help="output analysis results directory",
    )
    # Optional arguments
    parser.add_argument(
        "--dbpath",
        action="store",
        dest="dbpath",
        default=".pyani/pyanidb",
        help="path to pyani database",
    )
    parser.add_argument(
        "--runs",
        action="store_true",
        dest="show_runs",
        default=False,
        help="Report table of analysis runs in database",
    )
    parser.add_argument(
        "--genomes",
        action="store_true",
        dest="show_genomes",
        default=False,
        help="Report table of genomes in database",
    )
    parser.add_argument(
        "--runs_genomes",
        action="store_true",
        dest="show_runs_genomes",
        default=False,
        help="Report table of all genomes for each run in " + "database",
    )
    parser.add_argument(
        "--genomes_runs",
        action="store_true",
        dest="show_genomes_runs",
        default=False,
        help="Report table of all runs in which each genome "
        + "in the database participates",
    )
    parser.add_argument(
        "--run_results",
        action="store",
        dest="run_results",
        default=False,
        help="Report table of results for a pyani run",
    )
    parser.add_argument(
        "--run_matrices",
        action="store",
        dest="run_matrices",
        default=False,
        help="Report matrices of results for a pyani run",
    )
    parser.add_argument(
        "--formats",
        dest="formats",
        action="store",
        default=None,
        help="Output formats (in addition to .tab)",
    )
    parser.set_defaults(func=subcommands.subcmd_report)


# Process command-line
def parse_cmdline(argv=None):
    """Parse command-line arguments for script.

    The script offers a single main parser, with subcommands for the actions:

    download    - download all available NCBI assemblies below the passed
                  taxonomy ID
    index       - index genome sequence files in a subdirectory, for analysis
    createdb    - generate SQLite database for data and analysis results
    anim        - conduct ANIm analysis
    anib        - conduct ANIb analysis
    aniblastall - conduct ANIblastall analysis
    report      - generate output describing analyses, genomes, and results
    plot        - generate graphical output describing results
    classify    - produce graph-based classification of genomes on the basis
                  of ANI analysis
    """
    # Main parent parser
    parser_main = ArgumentParser(
        prog="pyani.py", formatter_class=ArgumentDefaultsHelpFormatter
    )
    subparsers = parser_main.add_subparsers(
        title="subcommands", description="valid subcommands", help="additional help"
    )

    # Common parsers
    parser_common = build_parser_common()
    parser_scheduler = build_parser_scheduler()
    parser_run_common = build_parser_run_common()

    # Add subcommand parsers
    build_parser_download(subparsers, parents=[parser_common])
    build_parser_index(subparsers, parents=[parser_common])
    build_parser_createdb(subparsers, parents=[parser_common])
    build_parser_db(subparsers, parents=[parser_common])
    build_parser_anim(
        subparsers, parents=[parser_common, parser_scheduler, parser_run_common]
    )
    build_parser_anib(
        subparsers, parents=[parser_common, parser_scheduler, parser_run_common]
    )
    build_parser_aniblastall(
        subparsers, parents=[parser_common, parser_scheduler, parser_run_common]
    )
    build_parser_report(subparsers, parents=[parser_common])
    build_parser_plot(subparsers, parents=[parser_common])
    build_parser_classify(subparsers, parents=[parser_common])

    # Parse arguments
    if argv is None:
        argv = sys.argv[1:]
    return parser_main.parse_args(argv)
