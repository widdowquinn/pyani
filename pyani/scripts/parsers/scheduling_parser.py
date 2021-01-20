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
# Cathedral Street,
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
"""Provides parser for arguments common to job scheduling."""

from argparse import ArgumentParser


def build() -> ArgumentParser:
    """Return the common argument parser for job scheduling.

    :param subps:  collection of subparsers in main parser
    :param parents:  parsers from which arguments are inherited

    Common arguments are:

    """
    parser = ArgumentParser(add_help=False)
    parser.add_argument(
        "--scheduler",
        dest="scheduler",
        action="store",
        default="multiprocessing",
        choices=["multiprocessing", "SGE", "SLURM"],
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
        #"--SGEgroupsize",
        "--groupsize",
        dest="sgegroupsize",
        action="store",
        default=10000,
        type=int,
        help="Number of jobs to place in an SGE array group " "(default 10000)",
    )
    parser.add_argument(
        #"--SGEargs",
        "--hpcargs",
        dest="sgeargs",
        action="store",
        default=None,
        type=str,
        help="Additional arguments for qsub/sbatch",
    )
    parser.add_argument(
        "--jobprefix",
        dest="jobprefix",
        action="store",
        default="PYANI",
        help="Prefix for SGE jobs (default PYANI).",
    )
    return parser
