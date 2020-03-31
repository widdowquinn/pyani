#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2017-2019
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
# Copyright (c) 2017-2019 The James Hutton Institute
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
"""Module providing a script logger.

This wraps the builtin logging module classes
"""

import logging
import sys
import time

from argparse import Namespace
from pathlib import Path


def build_logger(name: str, args: Namespace) -> logging.Logger:
    """Return a logger for this script.

    :param name:  str, name for logger
    :param args:  Namespace, command-line arguments

    Instantiates a logger for the script, and adds basic info.
    """
    logger = logging.getLogger("{}: {}".format(name, time.asctime))
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter("%(levelname)s: %(message)s")
    err_handler.setFormatter(err_formatter)

    # Verbose output?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # If a logfile was specified, use it
    if args.logfile is not None:
        logdir = args.logfile.parents[0]
        try:
            if not logdir == Path.cwd():
                logdir.mkdir(exist_ok=True)
            logstream = args.logfile.open("w")
        except OSError:
            logger.error("Could not open %s for logging", args.logfile, exc_info=True)
            raise SystemExit(1)
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)

    # Report arguments
    args.cmdline = " ".join(sys.argv)
    logger.info("Processed arguments: %s", args)
    logger.info("command-line: %s", args.cmdline)

    # Add citation information
    logger.info("CITATION INFO")
    pyani_citation = [
        "If you use pyani in your published work, please cite the",
        "following publication:",
        "\tPritchard, L., Glover, R. H., Humphris, S., Elphinstone, J. G., ",
        "\t& Toth, I.K. (2016) 'Genomics and taxonomy in diagnostics for ",
        "\tfood security: soft-rotting enterobacterial plant pathogens.'",
        "\tAnalytical Methods, 8(1), 12–24. http://doi.org/10.1039/C5AY02550H",
    ]
    for line in pyani_citation:
        logger.info(line)

    # Add dependency citations
    dep_citations = [
        "The authors of pyani gratefully acknowledge its dependence on",
        "the following bioinformatics software:",
        "\tMUMmer3: S. Kurtz, A. Phillippy, A.L. Delcher, M. Smoot, M. Shumway,",
        "\tC. Antonescu, and S.L. Salzberg (2004), 'Versatile and open software",
        "\tfor comparing large genomes' Genome Biology 5:R12",
        "\tBLAST+: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J.,",
        "\tBealer K., & Madden T.L. (2008) 'BLAST+: architecture and applications.'",
        "\tBMC Bioinformatics 10:421.",
        "\tBLAST: Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J.,",
        "\tZhang, Z., Miller, W. & Lipman, D.J. (1997) 'Gapped BLAST and PSI-BLAST:",
        "\ta new generation of protein database search programs.' Nucleic Acids Res.",
        "\t25:3389-3402\n",
    ]
    for line in dep_citations:
        logger.info(line)

    return logger
