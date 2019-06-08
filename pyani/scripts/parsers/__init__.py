# -*- coding: utf-8 -*-
"""Module providing parser definitions

(c) The James Hutton Institute 2017-2019

Author: Leighton Pritchard
Contact: leighton.pritchard@hutton.ac.uk
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

Copyright (c) 2017-2019 The James Hutton Institute
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

from pyani.scripts.parsers import (
    download_parser,
    index_parser,
    createdb_parser,
    db_parser,
    anim_parser,
    anib_parser,
    aniblastall_parser,
    report_parser,
    plot_parser,
    classify_parser,
    scheduling_parser,
    common_parser,
    run_common_parser,
)


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
    parser_common = common_parser.build()
    parser_scheduler = scheduling_parser.build()
    parser_run_common = run_common_parser.build()

    # Add subcommand parsers
    download_parser.build(subparsers, parents=[parser_common])
    index_parser.build(subparsers, parents=[parser_common])
    createdb_parser.build(subparsers, parents=[parser_common])
    db_parser.build(subparsers, parents=[parser_common])
    anim_parser.build(
        subparsers, parents=[parser_common, parser_scheduler, parser_run_common]
    )
    anib_parser.build(
        subparsers, parents=[parser_common, parser_scheduler, parser_run_common]
    )
    aniblastall_parser.build(
        subparsers, parents=[parser_common, parser_scheduler, parser_run_common]
    )
    report_parser.build(subparsers, parents=[parser_common])
    plot_parser.build(subparsers, parents=[parser_common])
    classify_parser.build(subparsers, parents=[parser_common])

    # Parse arguments
    if argv is None:
        argv = sys.argv[1:]
    return parser_main.parse_args(argv)
