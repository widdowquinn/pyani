# -*- coding: utf-8 -*-
# (c) The James Hutton Institute 2019
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
"""Code for handling BLAST output files."""


def parse_blasttab(fhandle):
    """Return the passed BLAST tab output file as a list of lists.

    :param fhandle:  filehandle, contains BLAST output file

    This is used when testing for conserved BLAST output, as the
    exact format of the BLAST result can depend on the software version.
    For instance, the locally-installed version may be BLASTN+ 2.6.0,
    which reports match identity to 3sf, and the version in CI may be
    BLASTN+ 2.2.28, which reports to 2sf.

    Returning a list of lines, parsed into the appropriate data type,
    allows for direct comparison of line content independent of formatting.
    """
    retval = []
    for line in fhandle.readlines():
        splitline = line.split("\t")
        data = splitline[:2]  # First two columns are strings
        data += [float(_) for _ in splitline[2:]]  # The rest are numeric
    return retval
