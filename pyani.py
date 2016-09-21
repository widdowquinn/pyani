#!/usr/bin/env python3
#
# pyani.py
#
# This script uses the pyani module to produce ANI analyses and classifications
# of prokaryotic genome sequences (draft or complete).
#
# (c) The James Hutton Institute 2016
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
# Copyright (c) 2016 The James Hutton Institute
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

import logging
import os
import shutil
import sys
import time
import traceback

from argparse import ArgumentParser

# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Process command-line
def parse_cmdline(args):
    """Parse command-line arguments for script.
    The script offers a single main parser, with subcommands for the actions:
    
    classify - produce a graph-based classification of each input genome
    """
    # Main parent parser
    parser_main = ArgumentParser(prog="pyani.py")
    subparsers = parser_main.add_subparsers(title="subcommands",
                                            description="valid subcommands",
                                            help="additional help")

    # A 'common' parser, with shared options for all subcommands
    parser_common = ArgumentParser(add_help=False)
    parser_common.add_argument('-l', '--logfile', dest='logfile',
                               action='store', default=None,
                               help='logfile location')
    parser_common.add_argument('-v', '--verbose', action='store_true',
                               dest='verbose', default=False,
                               help='report verbose progress to log')
    parser_common.add_argument('-i', '--indir', action='store',
                               dest='indir', default=None,
                               help='input directory')
    parser_common.add_argument('-o', '--outdir', action='store',
                               dest='outdir', default=None,
                               help='output directory')

    # Subcommand parsers
    parser_classify = subparsers.add_parser('classify',
                                            parents=[parser_common])

    # Genome classification options
    parser_classify.add_argument('--labels', dest='labelfile',
                                 action='store', default=None,
                                 help='file with labels for input genomes')
    parser_classify.add_argument('--cov_threshold', dest='cov_threshold',
                                 action='store', type=float, default=0.5,
                                 help='minimum %%coverage for an edge')
    parser_classify.add_argument('--id_threshold', dest='id_threshold',
                                 action='store', type=float, default=0.5,
                                 help='minimum %%identity for an edge')
    parser_classify.add_argument('--resolution', dest='resolution',
                                 action='store', type=int, default=1500,
                                 help='number of identity thresholds to test')
    
    # Parse arguments
    return parser_main.parse_args()


# Classify input genomes on basis of ANI coverage and identity output
def subcmd_classify(args):
    """Take a directory of pyani output files, and generate a series of
    classifications of the input data.
    """
    # Process the input directory and load data


    
###
# Run as script
if __name__ == '__main__':

    # Parse command-line
    assert len(sys.argv) > 1, 'pyani.py requires a valid subcommand'
    args = parse_cmdline(sys.argv)
    subcmd = sys.argv[1]

    # Set up logging
    logger = logging.getLogger('pyani.py: %s' % time.asctime())
    t0 = time.time()
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
            logger.error('Could not open %s for logging' %
                         args.logfile)
            logger.error(last_exception())
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info('Processed arguments: %s' % args)
    logger.info('command-line: %s' % ' '.join(sys.argv))

    # Define subcommand main functions, and distribute
    subcmds = {'classify': subcmd_classify}
    try:
        subcmds[subcmd](args)
    except KeyError:
        logger.error("Subcommand {0} not recognised".format(subcmd))
        raise NotImplementedError
    else:
        logger.error("Could not execute subcommand {0}".format(subcmd))
        logger.error(last_exception())
        sys.exit(1)

    # Exit cleanly (POSIX)
    sys.exit(0)
