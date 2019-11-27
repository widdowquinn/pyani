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
"""Wrapper for MUMmer 3.23 delta-filter script.

It is required in order to catch STDOUT ahead of the SGE job
runner so that pyani can run on SGE/OGE scheduling systems.

The wrapper does not modify the output of delta-filter, and
is called in exactly the same way, passing arguments through
directly. The only departure from this is that the first
argument is the path to delta-filter, and the final
argument denotes the output filtered delta file path, so that
redirection (which no longer works with SGE/OGE) is not
necessary.

For example, the delta-filter command

delta-filter [options] <delta file> > <filtered delta file>

becomes

delta_filter_wrapper.py delta-filter [options] <delta file> <filtered delta file>

This wrapper is not very robust, but will be improved in later
versions of pyani.
"""

import shlex
import subprocess
import sys


def run_main() -> int:
    """Run main process for delta_filter_wrapper.py."""
    # Parse command-line
    df_exe = sys.argv[1]
    args = sys.argv[2:-1]
    outfname = sys.argv[-1]

    sys.stdout.write("script called with %s" % sys.argv)

    # Run delta-filter, routing output to the named file
    cmd = [shlex.quote(df_exe)] + [shlex.quote(_) for _ in args]
    with open(outfname, "w") as ofh:
        ofh.write(
            subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
                shell=False,
            ).stdout.decode("utf-8")
        )

    # Exit
    return 0
