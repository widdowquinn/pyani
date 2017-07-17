# -*- coding: utf-8 -*-
"""parsers.py

This module provides command-line (and subcommand) parsers for the pyani.py
script.

- download:        download genome files from NCBI
- classify:        classify input sequences on the basis of ANI output

(c) The James Hutton Institute 2016-2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

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

Copyright (c) 2016-2017 The James Hutton Institute

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

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from . import subcommands


# Common parser for all subcommands
def build_parser_common():
    """Returns the common argument parser for all script subcommands

    Common arguments are:

    -l, --logfile (specify logfile output path)
    -v, --verbose (produce verbose output on STDOUT)
    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument('-l', '--logfile', dest='logfile',
                        action='store', default=None,
                        help='logfile location')
    parser.add_argument('-v', '--verbose', action='store_true',
                        dest='verbose', default=False,
                        help='report verbose progress to log')    
    return parser


# Subcommand parsers
def build_parser_download(subparsers, parents=None):
    """Returns a command-line parser for the download subcommand

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
    parser = subparsers.add_parser('download', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    # Required positional argument: output directory
    parser.add_argument(action='store',
                        dest='outdir', default=None,
                        help='output directory')
    # Required arguments for NCBI query
    parser.add_argument("-t", "--taxon", dest="taxon",
                        action="store", default=None,
                        help="NCBI taxonomy IDs (required, " +
                        "comma-separated list)",
                        required=True)
    parser.add_argument("--email", dest="email",
                        action="store", default=None,
                        help="Email associated with NCBI queries (required)",
                        required=True)
    # Arguments controlling connection to NCBI for download
    parser.add_argument("--retries", dest="retries",
                        action="store", default=20,
                        help="Number of Entrez retry attempts per " +\
                        "request")
    parser.add_argument("--batchsize", dest="batchsize",
                        action="store", default=10000,
                        help="Entrez record return batch size")
    parser.add_argument("--timeout", dest="timeout",
                        action="store", default=10,
                        help="Timeout for URL connection (s)")
    # Arguments controlling local filehandling
    parser.add_argument("-f", "--force", dest="force",
                        action="store_true", default=False,
                        help="Allow download to existing directory")
    parser.add_argument("--noclobber", dest="noclobber",
                        action="store_true", default=False,
                        help="Don't replace existing files")
    # Names for output files
    parser.add_argument("--labels", dest="labelfname",
                        action="store", default="labels.txt",
                        help="Filename for labels file")
    parser.add_argument("--classes", dest="classfname",
                        action="store", default="classes.txt",
                        help="Filename for classes file")
    parser.set_defaults(func=subcommands.subcmd_download)


def build_parser_classify(subparsers, parents=None):
    """Returns a command-line parser for the classify subcommand

    The classify subcommand takes specific arguments:

    -o, --outdir    (directory to write classification output)
    --labels        (path to input file labels)
    --cov_min       (minimum coverage threshold for an edge)
    --id_min        (minimum identity threshold for an edge)
    --resolution    (number of identity thresholds to test)
    """
    parser = subparsers.add_parser('classify', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    # Required positional argument: input directory
    parser.add_argument(action='store',
                        dest='indir', default=None,
                        help='input directory')
    # Output directory, defaults to input directory indir
    parser.add_argument('-o', '--outdir', action='store',
                        dest='outdir', default=None,
                        help='output directory')
    # Label file, defaults to indir/labels.txt
    parser.add_argument('--labels', dest='labelfile',
                        action='store', default=None,
                        help='file with labels for input genomes')
    # Parameters for classification: minimum %coverage, %identity,
    # and the resolution of thresholds to test
    parser.add_argument('--cov_min', dest='cov_min',
                        action='store', type=float, default=0.5,
                        help='minimum %%coverage for an edge')
    parser.add_argument('--id_min', dest='id_min',
                        action='store', type=float, default=0.8,
                        help='minimum %%identity for an edge')
    parser.add_argument('--resolution', dest='resolution',
                        action='store', type=int, default=1500,
                        help='number of identity thresholds to test')
    parser.set_defaults(func=subcommands.subcmd_classify)
    

def build_parser_index(subparsers, parents=None):
    """Returns a command-line parser for the index subcommand

    The index subcommand takes a single positional argument:

    indir           (directory containing input genome sequence files)
    """
    parser = subparsers.add_parser('index', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    # Required positional argument: input directory
    parser.add_argument(action='store',
                        dest='indir', default=None,
                        help='input directory')
    parser.set_defaults(func=subcommands.subcmd_index)


def build_parser_createdb(subparsers, parents=None):
    """Returns a command-line parser for the createdb subcommand
    """
    parser = subparsers.add_parser('createdb', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    # Path to database (default: .pyani/pyanidb)
    parser.add_argument("--dbpath", action='store',
                        dest='dbpath', default='.pyani/pyanidb',
                        help='path to pyani database')
    parser.add_argument("-f", "--force", action='store_true',
                        dest='force', default=False,
                        help='force creation of new empty database')
    parser.set_defaults(func=subcommands.subcmd_createdb)


def build_parser_anim(subparsers, parents=None):
    """Returns a command-line parser for the anim subcommand
    """
    parser = subparsers.add_parser('anim', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    parser.set_defaults(func=subcommands.subcmd_anim)


def build_parser_anib(subparsers, parents=None):
    """Returns a command-line parser for the anib subcommand
    """
    parser = subparsers.add_parser('anib', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    parser.set_defaults(func=subcommands.subcmd_anib)


def build_parser_aniblastall(subparsers, parents=None):
    """Returns a command-line parser for the aniblastall subcommand
    """
    parser = subparsers.add_parser('aniblastall', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    parser.set_defaults(func=subcommands.subcmd_aniblastall)


def build_parser_render(subparsers, parents=None):
    """Returns a command-line parser for the render subcommand
    """
    parser = subparsers.add_parser('render', parents=parents,
                                   formatter_class=\
                                   ArgumentDefaultsHelpFormatter)
    parser.set_defaults(func=subcommands.subcmd_render)
    
    
# Process command-line
def parse_cmdline():
    """Parse command-line arguments for script.

    The script offers a single main parser, with subcommands for the actions:

    download - download all available NCBI assemblies below the passed
               taxonomy ID
    classify - produce a graph-based classification of each input genome
    """
    # Main parent parser
    parser_main = ArgumentParser(prog="pyani.py",
                                 formatter_class=ArgumentDefaultsHelpFormatter)
    subparsers = parser_main.add_subparsers(title="subcommands",
                                            description="valid subcommands",
                                            help="additional help")

    # Common parser, included with all the subcommand parsers
    parser_common = build_parser_common()

    # Add subcommand parsers
    build_parser_download(subparsers, parents=[parser_common])
    build_parser_index(subparsers, parents=[parser_common])
    build_parser_createdb(subparsers, parents=[parser_common])
    build_parser_anim(subparsers, parents=[parser_common])
    build_parser_anib(subparsers, parents=[parser_common])
    build_parser_aniblastall(subparsers, parents=[parser_common])
    build_parser_render(subparsers, parents=[parser_common])
    build_parser_classify(subparsers, parents=[parser_common])
    
    # Parse arguments
    return parser_main.parse_args()
