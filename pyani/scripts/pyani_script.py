#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2016-2019
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
# Copyright (c) 2016-2019 The James Hutton Institute
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
"""Implements the pyani script for classifying prokaryotic genomes."""

import logging
import sys
import time
import subprocess

from typing import List, Optional

# from . import pyani_config

from pyani.logger import config_logger
from pyani.pyani_tools import termcolor

from .parsers import parse_cmdline
from .. import __version__

from pyani.anib import get_version as get_blast_version
from pyani.aniblastall import get_version as get_blastall_version
from pyani.anim import get_version as get_nucmer_version


CITATION_INFO = [
    termcolor(
        "If you use pyani in your work, please cite the following publication:",
        "green",
    ),
    termcolor(
        "\tPritchard, L., Glover, R. H., Humphris, S., Elphinstone, J. G.,",
        "yellow",
    ),
    termcolor(
        "\t& Toth, I.K. (2016) 'Genomics and taxonomy in diagnostics for", "yellow"
    ),
    termcolor(
        "\tfood security: soft-rotting enterobacterial plant pathogens.'", "yellow"
    ),
    termcolor(
        "\tAnalytical Methods, 8(1), 12–24. http://doi.org/10.1039/C5AY02550H",
        "yellow",
    ),
]


# Main function
def run_main(argv: Optional[List[str]] = None) -> int:
    """Run main process for pyani.py script.

    :param argv:
    """
    # If we need to (i.e. a namespace isn't passed), parse the command-line
    if argv is None:
        args = parse_cmdline()
    else:
        args = parse_cmdline(argv)
    # Catch execution with no arguments
    if len(sys.argv) == 1:
        sys.stderr.write(f"pyani version: {__version__}\n\n")
        return 0
    # Catch requests for citation and version information
    if sys.argv[1].startswith("-"):
        if args.citation:
            sys.stderr.write(f"pyani version: {__version__}\n\n")
            for line in CITATION_INFO:
                sys.stderr.write(f"{line}\n")
            return 0
        if args.version:
            sys.stderr.write(f"pyani version: {__version__}\n\n")
            return 0
        return 0

    # Set up logging
    time0 = time.time()
    logger = logging.getLogger(__name__)
    config_logger(args)

    # Boilerplate for log
    logger.info("Processed arguments: %s", args)
    args.cmdline = " ".join(sys.argv)
    logger.info("command-line: %s", args.cmdline)
    add_log_headers()

    # Run the subcommand
    returnval = args.func(args)
    logger.info(
        termcolor("Completed. Time taken: %.3f", bold=True), (time.time() - time0)
    )

    # Was citation info requested?
    if args.citation:
        sys.stdout.write("\n".join(CITATION_INFO) + "\n\n")

    return returnval


def add_log_headers():
    """Add headers to log output."""
    logger = logging.getLogger(__name__)

    # Add citation information to log
    logger.info(termcolor("CITATION INFO", bold=True))
    for line in CITATION_INFO:
        logger.info(line)

    # Add dependency citations
    logger.info(termcolor("DEPENDENCIES", bold=True))
    dep_citations = [
        "The authors of pyani gratefully acknowledge its dependence on",
        "the following bioinformatics software:",
        f"\t{termcolor('MUMmer3', 'cyan')}: S. Kurtz, A. Phillippy, A.L. Delcher, M. Smoot, M. Shumway,",
        "\tC. Antonescu, and S.L. Salzberg (2004), 'Versatile and open software",
        "\tfor comparing large genomes' Genome Biology 5:R12",
        f"\t{termcolor('BLAST+', 'cyan')}: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J.,",
        "\tBealer K., & Madden T.L. (2008) 'BLAST+: architecture and applications.'",
        "\tBMC Bioinformatics 10:421.",
        f"\t{termcolor('BLAST', 'cyan')}: Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J.,",
        "\tZhang, Z., Miller, W. & Lipman, D.J. (1997) 'Gapped BLAST and PSI-BLAST:",
        "\ta new generation of protein database search programs.' Nucleic Acids Res.",
        "\t25:3389-3402",
        f"\t{termcolor('Biopython', 'cyan')}: Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A,",
        "\tFriedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL",
        "\t(2009) Biopython: freely available Python tools for computational",
        "\tmolecular biology and bioinformatics. Bioinformatics, 25, 1422-1423",
    ]
    for line in dep_citations:
        logger.info(line)
