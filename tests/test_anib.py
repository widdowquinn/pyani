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
                                        os.path.join('tests', 'test_output',
                                                     'anib', 'makeblastdb',
                                                     'NC_002696.fna')])
        self.blastdbfnames = [os.path.join(self.seqdir, fname) for fname in
                              ('NC_002696.fna', 'NC_010338.fna')]
        self.blastdbtgt = [(' '.join(['makeblastdb -dbtype nucl -in',
                                      'tests/test_input/sequences/NC_002696.fna',
                                      '-title NC_002696 -out',
                                      'tests/test_output/anib/NC_002696.fna']),
                            'tests/test_output/anib/NC_002696.fna'),
                           (' '.join(['makeblastdb -dbtype nucl -in',
                                      'tests/test_input/sequences/NC_010338.fna',
                                      '-title NC_010338 -out',
                                      'tests/test_output/anib/NC_010338.fna']),
                            'tests/test_output/anib/NC_010338.fna')]
        self.blastdbtgtlegacy = [(' '.join(['formatdb -p F -i',
                                            'tests/test_output/anib/NC_002696.fna',
                                            '-t NC_002696']),
                                  'tests/test_output/anib/NC_002696.fna'),
                                 (' '.join(['formatdb -p F -i',
                                            'tests/test_output/anib/NC_010338.fna',
                                            '-t NC_010338']),
                                  'tests/test_output/anib/NC_010338.fna')]
        self.blastncmd = ' '.join(["blastn -out",
                                   os.path.join("tests", "test_output", "anib",
                                                "NC_002696_vs_NC_010338.blast_tab"),
                                   "-query tests/test_input/sequences/NC_002696.fna",
                                   "-db tests/test_input/sequences/NC_010338.fna",
                                   "-xdrop_gap_final 150 -dust no -evalue 1e-15",
                                   "-max_target_seqs 1 -outfmt '6 qseqid sseqid",
                                   "length mismatch pident nident qlen slen qstart",
                                   "qend sstart send positive ppos gaps' -task blastn"])
        self.blastallcmd = ' '.join(['blastall -p blastn -o',
                                     os.path.join('tests', 'test_output', 'anib',
                                                  'NC_002696_vs_NC_010338.blast_tab'),
                                     '-i tests/test_input/sequences/NC_002696.fna',
                                     '-d tests/test_input/sequences/NC_010338.fna',
                                     '-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8'])
        self.blastntgt = [' '.join(["blastn -out",
                                    os.path.join("tests", "test_output", "anib",
                                                 "NC_002696_vs_NC_010338.blast_tab"),
                                    "-query tests/test_input/sequences/NC_002696.fna",
                                    "-db tests/test_input/sequences/NC_010338.fna",
                                    "-xdrop_gap_final 150 -dust no -evalue 1e-15",
                                    "-max_target_seqs 1 -outfmt '6 qseqid sseqid",
                                    "length mismatch pident nident qlen slen qstart",
                                    "qend sstart send positive ppos gaps' -task blastn"]),
                          ' '.join(["blastn -out",
                                    os.path.join("tests", "test_output", "anib",
                                                 "NC_010338_vs_NC_002696.blast_tab"),
                                    "-query tests/test_input/sequences/NC_010338.fna",
                                    "-db tests/test_input/sequences/NC_002696.fna",
                                    "-xdrop_gap_final 150 -dust no -evalue 1e-15",
                                    "-max_target_seqs 1 -outfmt '6 qseqid sseqid length",
                                    "mismatch pident nident qlen slen qstart qend",
                                    "sstart send positive ppos gaps' -task blastn"])]
        self.blastalltgt = [' '.join(['blastall -p blastn -o',
                                      os.path.join('tests', 'test_output', 'anib',
                                                   'NC_002696_vs_NC_010338.blast_tab'),
                                      '-i tests/test_input/sequences/NC_002696.fna',
                                      '-d tests/test_input/sequences/NC_010338.fna',
                                      '-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8']),
                            ' '.join(['blastall -p blastn -o',
                                      os.path.join('tests', 'test_output', 'anib',
                                                   'NC_010338_vs_NC_002696.blast_tab'),
                                      '-i tests/test_input/sequences/NC_010338.fna',
                                      '-d tests/test_input/sequences/NC_002696.fna',
                                      '-X 150 -q -1 -F F -e 1e-15 -b 1 -v 1 -m 8'])]
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

    def test_blastdb_commands(self):
        """generate both BLAST+ and legacy BLAST db commands."""
        # BLAST+
        cmds = anib.generate_blastdb_commands(self.blastdbfnames, self.outdir,
                                              mode="ANIb")
        assert_equal(cmds, self.blastdbtgt)
        # legacy
        cmds = anib.generate_blastdb_commands(self.blastdbfnames, self.outdir,
                                              mode="ANIblastall")
        assert_equal(cmds, self.blastdbtgtlegacy)

    def test_blastn_generation(self):
        """generate BLASTN+ command-line."""
        cmd = anib.construct_blastn_cmdline(self.blastdbfnames[0],
                                            self.blastdbfnames[1],
                                            self.outdir)
        assert_equal(cmd, self.blastncmd)

    def test_blastall_generation(self):
        """generate legacy BLASTN command-line."""
        cmd = anib.construct_blastall_cmdline(self.blastdbfnames[0],
                                              self.blastdbfnames[1],
                                              self.outdir)
        print(cmd)
        assert_equal(cmd, self.blastallcmd)

    def test_blastn_commands(self):
        """Generate both BLASTN+ and legacy BLASTN commands."""
        # BLAST+
        cmds = anib.generate_blastn_commands(self.blastdbfnames, self.outdir,
                                             mode="ANIb")
        assert_equal(cmds, self.blastntgt)
        cmds = anib.generate_blastn_commands(self.blastdbfnames, self.outdir,
                                             mode="ANIblastall")
        assert_equal(cmds, self.blastalltgt)

