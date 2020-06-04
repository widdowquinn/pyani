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
"""Module providing support for package-level logging."""

import logging
import logging.config
import re
import sys

from argparse import Namespace
from pathlib import Path
from typing import Optional


class NoColorFormatter(logging.Formatter):

    """Log formatter that strips terminal colour escape codes from the log message."""

    ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

    def format(self, record):
        """Return logger message with terminal escapes removed."""
        return "[%s] [%s]: %s" % (
            record.levelname,
            record.name,
            re.sub(self.ANSI_RE, "", record.msg % record.args),
        )


def config_logger(args: Optional[Namespace] = None) -> None:
    """Configure package-level logging.

    :param args: CLI namespace; logfile is used to create a logfile,
        verbose and debug control logging level.

    We configure a logger at package level, from which the module will
    inherit. If CLI args are provided, these are used to define output
    streams, and logging level.
    """
    # Default logger for this module
    logger = logging.getLogger(__package__)
    logger.setLevel(logging.DEBUG)

    # Create and add STDERR handler
    errformatter = logging.Formatter("[%(levelname)s] [%(name)s]: %(message)s")
    errhandler = logging.StreamHandler(sys.stderr)
    if args is not None and args.verbose:
        errhandler.setLevel(logging.INFO)
    elif args is not None and args.debug:
        errhandler.setLevel(logging.DEBUG)
    else:
        errhandler.setLevel(logging.WARNING)
    errhandler.setFormatter(errformatter)
    logger.addHandler(errhandler)

    # If args.logfile is provided, add a FileHandler for logfile
    if args is not None and args.logfile is not None:
        logdir = args.logfile.parents[0]
        # Check that output directory exists and, if not, create it
        try:
            if not logdir == Path.cwd():
                logdir.mkdir(exist_ok=True)
        except OSError:
            logger.error(
                "Could not create log directory %s (exiting)", logdir, exc_info=True
            )
            raise SystemExit(1)

        # Create logfile handler
        logformatter = NoColorFormatter()
        loghandler = logging.FileHandler(args.logfile, mode="w", encoding="utf8")
        if args.debug:
            loghandler.setLevel(logging.DEBUG)
        else:
            loghandler.setLevel(logging.INFO)
        loghandler.setFormatter(logformatter)
        logger.addHandler(loghandler)
