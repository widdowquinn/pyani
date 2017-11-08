#!/usr/bin/env python

"""test_anib.py

Test anib.py module.

These tests are intended to be run from the repository root using:

nosetests -v

print() statements will be caught by nosetests unless there is an
error. They can also be recovered with the -s option.

(c) The James Hutton Institute 2017
Author: Leighton Pritchard

Contact:
leighton.pritchard@hutton.ac.uk

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

import os
import unittest

from nose.tools import (assert_equal,)

from pyani import (anib,)

class TestBLASTCmdline(unittest.TestCase):

    """Class defining tests of BLAST command-line generation."""

    def setUp(self):
        """Set parameters for tests."""
        self.indir = os.path.join('tests', 'test_input', 'anib')
        self.outdir = os.path.join('tests', 'test_output', 'anib')
        self.seqdir = os.path.join('tests', 'test_input', 'sequences')
        self.fmtdboutdir = os.path.join(self.outdir, 'formatdb')
        self.fmtdbcmd = ' '.join(["formatdb -p F -i",
                                  "tests/test_output/anib/formatdb/NC_002696.fna",
                                  "-t NC_002696"])
        self.makeblastdbdir = os.path.join(self.outdir, 'makeblastdb')
        self.makeblastdbcmd = ' '.join(["makeblastdb -dbtype nucl -in",
                                        "tests/test_input/sequences/NC_002696.fna",
                                        "-title NC_002696 -out",
                                        "tests/test_output/anib/makeblastdb/NC_002696.fna"])
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(self.fmtdboutdir, exist_ok=True)
        os.makedirs(self.makeblastdbdir, exist_ok=True)

    def test_formatdb_generation(self):
        """generate formatdb command-line."""
        cmd = anib.construct_formatdb_cmd(os.path.join(self.seqdir,
                                                       "NC_002696.fna"),
                                          self.fmtdboutdir)
        assert_equal(cmd[0], self.fmtdbcmd)  # correct command
        assert(os.path.isfile(cmd[1]))       # creates new file

    def test_makeblastdb_generation(self):
        """generate makeblastdb command-line."""
        cmd = anib.construct_makeblastdb_cmd(os.path.join(self.seqdir,
                                                          "NC_002696.fna"),
                                             self.makeblastdbdir)
        assert_equal(cmd[0], self.makeblastdbcmd)  # correct command
        assert(os.path.isfile(cmd[1]))             # creates new file
