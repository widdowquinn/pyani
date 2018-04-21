#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Implements the pyani script for classifying prokaryotic genomes.

(c) The James Hutton Institute 2017
Author: Leighton Pritchard

Contact: leighton.pritchard@hutton.ac.uk
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

Copyright (c) 2017 The James Hutton Institute
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

import logging
import sys
import time

from .parsers import parse_cmdline
from .. import __version__
from ..pyani_tools import last_exception


# Main function
def run_main(namespace=None):
    """Run main process for pyani.py script."""
    # If we need to (i.e. a namespace isn't passed), parse the command-line
    if namespace is None:
        args = parse_cmdline()
    else:
        args = namespace

    # Catch execution with no arguments
    if len(sys.argv) == 1:
        sys.stderr.write("pyani version: {0}\n".format(__version__))
        return 0

    # Set up logging
    logger = logging.getLogger('labbook.py: %s' % time.asctime())
    time0 = time.time()
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
        except OSError:
            logger.error('Could not open %s for logging', args.logfile)
            logger.error(last_exception())
            sys.exit(1)
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    args.cmdline = ' '.join(sys.argv)
    logger.info('Processed arguments: %s', args)
    logger.info('command-line: %s', args.cmdline)

    # Run the subcommand
    returnval = args.func(args, logger)
    logger.info('Completed. Time taken: %.3f', (time.time() - time0))
    return returnval
