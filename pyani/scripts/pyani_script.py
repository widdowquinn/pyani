#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
# (c) University of Strathclyde 2019
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@strath.ac.uk
#
# Leighton Pritchard,
# Strathclyde Institute for Pharmacy and Biomedical Sciences,
# Cathedral Street,
# Glasgow,
# G1 1XQ
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2016-2019 The James Hutton Institute
# Copyright (c) 2019 University of Strathclyde
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
"""Implements the pyani script for classifying prokaryotic genomes."""

import sys
import time

from logging import Logger
from typing import List, Optional

from .logger import build_logger
from .parsers import parse_cmdline
from .. import __version__


# Main function
def run_main(argv: Optional[List[str]] = None, logger: Optional[Logger] = None) -> int:
    """Run main process for pyani.py script.

    :param argv:
    :param logger:
    """
    # If we need to (i.e. a namespace isn't passed), parse the command-line
    if argv is None:
        args = parse_cmdline()
    else:
        args = parse_cmdline(argv)

    # Catch execution with no arguments
    if len(sys.argv) == 1:
        sys.stderr.write("pyani version: {0}\n".format(__version__))
        return 0

    # Set up logging
    time0 = time.time()
    if logger is None:
        logger = build_logger("pyani", args)

    # Run the subcommand
    returnval = args.func(args, logger)
    logger.info("Completed. Time taken: %.3f", (time.time() - time0))
    return returnval
